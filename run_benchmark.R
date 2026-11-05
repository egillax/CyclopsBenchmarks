# Cyclops 4.0 benchmark: legacy / default / super on HADES workloads
# ==================================================================
#
# Reads the data objects made by extract_cdm.R and times only the Cyclops calls, each fit in a fresh
# R process (callr) with one Cyclops library on the path:
#
#   ps    propensity model as CohortMethod::createPs() fits it (LR, 10-fold auto CV)
#   plp   LASSO as PatientLevelPrediction fits it (LR, 75% training split, 3-fold auto CV)
#   cox   Cox outcome model on the PS cohort: fixed-variance fit (1 and 8 threads) and 10-fold CV
#   sccs  SelfControlledCaseSeries::fitSccsModel(), two designs (season + calendar time; age + season)
#
# Configurations: develop with fitter = "default", "super" (LR only), "legacy", and optionally the
# active-set Newton branch (fitter "super", LR only). Fast first, slow last; a per-fit timeout and a
# total time budget keep the run within a few hours. Before the full runs, every fit runs once on a
# 5% subsample (smoke test); those timings drive the ETA printed after each fit.
#
# Output in <outputFolder>/<databaseLabel>/: results.csv (one row per fit), summary.csv, progress.txt,
# environment.txt, REPORT.md (aggregates only) and private/ (prepared data and coefficients: stay here).
#
#   Rscript run_benchmark.R                # everything: prepare, smoke test, full runs, report
#   Rscript run_benchmark.R report         # only rebuild summary.csv and REPORT.md

# ---- Settings (edit) ---------------------------------------------------------------------------

env <- function(name, default) { v <- Sys.getenv(name); if (nzchar(v)) v else default }
databaseLabel  <- env("BENCH_LABEL", "CCAE")
dataFolder     <- env("BENCH_DATA", file.path("cyclops_benchmark_data", databaseLabel))   # from extract_cdm.R
outputFolder   <- env("BENCH_OUT", "cyclops_benchmark_results")

libDevelop     <- env("BENCH_LIB_DEVELOP", "lib/cyclops-develop")   # Cyclops develop, installed by installCyclops() below
lib371         <- env("BENCH_LIB_371", "lib/cyclops-3.7.1")      # Cyclops 3.7.1 (CRAN), installed by installCyclops()
libNewton      <- env("BENCH_LIB_NEWTON", "lib/cyclops-newton")   # optional: the active-set Newton branch (skipped if absent)
libHades       <- env("BENCH_LIB_HADES", "")   # extra library with the HADES packages, if not in the default path

budgetHours    <- as.numeric(env("BENCH_BUDGET_HOURS", 4))       # total for the full runs; fits that would not fit are skipped (SKIPPED_BUDGET)
timeoutMinutes <- as.numeric(env("BENCH_TIMEOUT_MINUTES", 60))      # per fit (TIMEOUT)
smokeFraction  <- 0.05
maxReps        <- 3        # extra repetitions of fits under 5 minutes, while the budget allows
# PLP tasks: name:outcomeId:riskWindowEnd:removeSubjectsWithPriorOutcome, comma-separated; the data of task "plp" are
# in <dataFolder>/plpData (extract_cdm.R), of any other task in <dataFolder>/plpData_<name>
plpTasks       <- env("BENCH_PLP_TASKS", "plp:6:1095:1")
# Configurations in run order: develop default, develop super, CRAN 3.7.1, then the Newton branch (if installed; PS and
# PLP only); "legacy" (the 4.0 legacy fitter) can be added
configsToRun   <- strsplit(env("BENCH_CONFIGS", "default,super,3.7.1,newton"), ",")[[1]]
workloadsToRun <- strsplit(env("BENCH_WORKLOADS", "ps,plp,cox,sccs"), ",")[[1]]   # a subset, e.g. "ps,plp"
psThreads      <- as.integer(env("BENCH_PS_THREADS", 10))   # threads run folds in parallel: one per fold (10)
coxOutcomeId   <- NA       # NA: death (4) if it has >= 1000 events in the PS cohort, else AMI (3)

# Install develop (and the Newton branch from a source tarball, if you have one) into separate libraries:
installCyclops <- function(developRef = "develop", newtonTarball = NULL) {
  dir.create(lib371, recursive = TRUE, showWarnings = FALSE)
  # 3.7.1 is CRAN's current release (a binary on Windows and macOS); once it is archived, use
  # remotes::install_version("Cyclops", "3.7.1", lib = lib371, upgrade = "never", dependencies = FALSE)
  utils::install.packages("Cyclops", lib = lib371, repos = "https://cloud.r-project.org", dependencies = FALSE)
  if (utils::packageVersion("Cyclops", lib.loc = lib371) != "3.7.1") warning("lib371 does not hold Cyclops 3.7.1")
  dir.create(libDevelop, recursive = TRUE, showWarnings = FALSE)
  remotes::install_github("OHDSI/Cyclops", ref = developRef, lib = libDevelop, upgrade = "never", dependencies = FALSE)
  if (!is.null(newtonTarball)) {
    dir.create(libNewton, recursive = TRUE, showWarnings = FALSE)
    install.packages(newtonTarball, lib = libNewton, repos = NULL, type = "source")
  }
}

# ---- Setup -------------------------------------------------------------------------------------

# Andromeda (and R's temp files in the fit processes) in a folder next to the project, not the system temp
# drive: BENCH_TEMP, or ./tmp
tempFolder <- normalizePath(if (nzchar(Sys.getenv("BENCH_TEMP"))) Sys.getenv("BENCH_TEMP") else "tmp", mustWork = FALSE)
dir.create(tempFolder, recursive = TRUE, showWarnings = FALSE)
options(andromedaTempFolder = tempFolder)
Sys.setenv(TMPDIR = tempFolder, TMP = tempFolder, TEMP = tempFolder)   # inherited by child R processes


suppressPackageStartupMessages(library(dplyr))
if (nzchar(libHades)) .libPaths(c(normalizePath(libHades), .libPaths()))
folder <- file.path(outputFolder, databaseLabel)
private <- file.path(folder, "private")
dir.create(file.path(private, "coefs"), recursive = TRUE, showWarnings = FALSE)
resultsFile <- file.path(folder, "results.csv")
smokeFile <- file.path(private, "smoke.csv")
logFile <- file.path(folder, "progress.txt")
say <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  message(line)
  cat(line, "\n", file = logFile, append = TRUE)
}
hms <- function(s) if (is.na(s)) "?" else sprintf("%d:%02d", as.integer(s %/% 3600), as.integer((s %% 3600) %/% 60))
null2na <- function(x) if (is.null(x) || length(x) == 0) NA else x
libs <- c(develop = normalizePath(libDevelop, mustWork = FALSE), newton = normalizePath(libNewton, mustWork = FALSE),
          v371 = normalizePath(lib371, mustWork = FALSE))
haveNewton <- file.exists(file.path(libs[["newton"]], "Cyclops"))

# ---- Data preparation (once; a full and a smoke-test version of each workload) -----------------

.tidy <- function(covariates, covariateRef, analysisRef, n) {
  data <- Andromeda::andromeda(covariates = covariates, covariateRef = covariateRef, analysisRef = analysisRef)
  attr(data, "metaData") <- list(populationSize = n)
  class(data) <- "CovariateData"
  attr(class(data), "package") <- "FeatureExtraction"
  tidy <- FeatureExtraction::tidyCovariateData(data, minFraction = 0.001, normalize = TRUE, removeRedundancy = TRUE)
  Andromeda::close(data)
  tidy
}

.save <- function(outcomes, covariates, file, extra = list()) {
  out <- Andromeda::andromeda(outcomes = outcomes, covariates = covariates)
  info <- c(list(n_rows = nrow(outcomes), n_events = sum(outcomes$y > 0),
                 n_covariates = out$covariates %>% distinct(covariateId) %>% count() %>% pull(),
                 n_nonzero_entries = out$covariates %>% count() %>% pull()), extra)
  Andromeda::saveAndromeda(out, file)
  saveRDS(info, paste0(file, ".info.rds"))
  info
}

.subsample <- function(ids, fraction, seed = 1) {
  set.seed(seed)
  if (fraction >= 1) ids else sample(ids, min(length(ids), max(2000, round(length(ids) * fraction))))
}

preparePs <- function(fraction, file) {
  cm <- CohortMethod::loadCohortMethodData(file.path(dataFolder, "cmData.zip"))
  on.exit(Andromeda::close(cm))
  population <- collect(cm$cohorts)[, c("rowId", "treatment")]
  set.seed(0) # createPs(): at most 250,000 per arm
  keep <- unlist(lapply(c(1, 0), function(t) {
    ids <- population$rowId[population$treatment == t]
    if (length(ids) > 250000) sample(ids, 250000) else ids
  }))
  keep <- .subsample(keep, fraction)
  population <- population[population$rowId %in% keep, ]
  tidy <- .tidy(cm$covariates %>% filter(.data$rowId %in% local(population$rowId)), cm$covariateRef, cm$analysisRef,
                nrow(population))
  on.exit(Andromeda::close(tidy), add = TRUE)
  .save(data.frame(rowId = population$rowId, y = population$treatment), tidy$covariates, file)
}

plpTaskList <- lapply(strsplit(strsplit(plpTasks, ",")[[1]], ":"), function(x)
  list(name = x[1], design = if (x[1] == "plp") "plp" else paste0("plp_", x[1]),
       folder = if (x[1] == "plp") "plpData" else paste0("plpData_", x[1]),
       outcomeId = as.numeric(x[2]), riskWindowEnd = as.numeric(x[3]), removePrior = x[4] == "1"))
names(plpTaskList) <- sapply(plpTaskList, `[[`, "design")

preparePlp <- function(fraction, file, task) {
  plpData <- PatientLevelPrediction::loadPlpData(file.path(dataFolder, task$folder))
  population <- PatientLevelPrediction::createStudyPopulation(plpData, outcomeId = task$outcomeId,
    populationSettings = PatientLevelPrediction::createStudyPopulationSettings(
      washoutPeriod = 365, removeSubjectsWithPriorOutcome = task$removePrior, requireTimeAtRisk = FALSE,
      riskWindowStart = 1, riskWindowEnd = task$riskWindowEnd))
  y <- as.numeric(population$outcomeCount > 0)
  set.seed(42) # stratified 75% training split, as PLP's default split
  train <- unlist(lapply(c(0, 1), function(v) {
    ids <- population$rowId[y == v]
    sample(ids, round(length(ids) * 0.75))
  }))
  train <- .subsample(train, fraction)
  outcomes <- data.frame(rowId = population$rowId, y = y)
  outcomes <- outcomes[outcomes$rowId %in% train, ]
  cd <- plpData$covariateData
  tidy <- .tidy(cd$covariates %>% filter(.data$rowId %in% local(outcomes$rowId)), cd$covariateRef, cd$analysisRef,
                nrow(outcomes))
  on.exit(Andromeda::close(tidy))
  .save(outcomes, tidy$covariates, file)
}

prepareCox <- function(fraction, file) {
  cm <- CohortMethod::loadCohortMethodData(file.path(dataFolder, "cmData.zip"))
  on.exit(Andromeda::close(cm))
  outcomeId <- coxOutcomeId
  if (is.na(outcomeId)) {
    deaths <- collect(cm$outcomes) %>% filter(.data$outcomeId == 4) %>% nrow()
    outcomeId <- if (deaths >= 1000) 4 else 3
  }
  population <- CohortMethod::createStudyPopulation(cm, outcomeId = outcomeId,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      removeSubjectsWithPriorOutcome = TRUE, minDaysAtRisk = 1, riskWindowStart = 1, startAnchor = "cohort start",
      riskWindowEnd = 1825, endAnchor = "cohort start"))
  set.seed(0)
  keep <- unlist(lapply(c(1, 0), function(t) {
    ids <- population$rowId[population$treatment == t]
    if (length(ids) > 250000) sample(ids, 250000) else ids
  }))
  keep <- .subsample(keep, fraction)
  population <- population[population$rowId %in% keep, ]
  tidy <- .tidy(cm$covariates %>% filter(.data$rowId %in% local(population$rowId)), cm$covariateRef, cm$analysisRef,
                nrow(population))
  on.exit(Andromeda::close(tidy), add = TRUE)
  treatmentId <- -1 # the treatment as a covariate, not penalized
  treatment <- data.frame(rowId = population$rowId[population$treatment == 1], covariateId = treatmentId,
                          covariateValue = 1)
  tidy$treatment <- treatment
  covariates <- tidy$covariates %>% select("rowId", "covariateId", "covariateValue")
  tidy$allCovariates <- covariates
  Andromeda::appendToTable(tidy$allCovariates, tidy$treatment)
  outcomes <- data.frame(rowId = population$rowId, y = as.numeric(population$outcomeCount > 0),
                         time = population$survivalTime)
  .save(outcomes, tidy$allCovariates %>% arrange(.data$rowId), file, extra = list(treatmentId = treatmentId, outcomeId = outcomeId))
}

prepareSccs <- function(fraction, file, design) {
  sccsData <- SelfControlledCaseSeries::loadSccsData(file.path(dataFolder, "sccsData.zip"))
  on.exit(Andromeda::close(sccsData))
  population <- SelfControlledCaseSeries::createStudyPopulation(sccsData, outcomeId = 3,
    createStudyPopulationArgs = SelfControlledCaseSeries::createCreateStudyPopulationArgs(
      firstOutcomeOnly = TRUE, naivePeriod = 180))
  if (fraction < 1) {
    cases <- .subsample(unique(population$outcomes$caseId), fraction)
    population$outcomes <- population$outcomes[population$outcomes$caseId %in% cases, ]
    population$cases <- population$cases[population$cases$caseId %in% cases, ]
  }
  exposures <- c(1177480, 1115008) # ibuprofen, naproxen
  names(exposures) <- c("ibuprofen", "naproxen")
  eras <- c(lapply(names(exposures), function(n) SelfControlledCaseSeries::createEraCovariateSettings(
      label = paste("Exposure", n), includeEraIds = exposures[[n]], start = 0, end = 0, endAnchor = "era end",
      profileLikelihood = TRUE)),
    list(SelfControlledCaseSeries::createEraCovariateSettings(label = "Pre-exposure", includeEraIds = exposures,
      stratifyById = TRUE, start = -30, end = -1, endAnchor = "era start")))
  minCases <- if (fraction < 1) 100 else 10000
  args <- if (design == "season_calendar") {
    SelfControlledCaseSeries::createCreateSccsIntervalDataArgs(eraCovariateSettings = eras,
      seasonalityCovariateSettings = SelfControlledCaseSeries::createSeasonalityCovariateSettings(),
      calendarTimeCovariateSettings = SelfControlledCaseSeries::createCalendarTimeCovariateSettings(),
      minCasesForTimeCovariates = minCases)
  } else {
    SelfControlledCaseSeries::createCreateSccsIntervalDataArgs(eraCovariateSettings = eras,
      ageCovariateSettings = SelfControlledCaseSeries::createAgeCovariateSettings(),
      seasonalityCovariateSettings = SelfControlledCaseSeries::createSeasonalityCovariateSettings(),
      minCasesForTimeCovariates = minCases)
  }
  interval <- SelfControlledCaseSeries::createSccsIntervalData(population, sccsData, args)
  SelfControlledCaseSeries::saveSccsIntervalData(interval, file)
  interval <- SelfControlledCaseSeries::loadSccsIntervalData(file)
  on.exit(Andromeda::close(interval), add = TRUE)
  info <- list(n_rows = interval$outcomes %>% count() %>% pull(),
               n_events = interval$outcomes %>% summarise(n = sum(.data$y, na.rm = TRUE)) %>% pull(),
               n_covariates = interval$covariates %>% distinct(covariateId) %>% count() %>% pull(),
               n_nonzero_entries = interval$covariates %>% count() %>% pull(),
               cases = length(unique(population$outcomes$caseId)))
  saveRDS(info, paste0(file, ".info.rds"))
  info
}

datasets <- data.frame(
  workload = c("ps", rep("plp", length(plpTaskList)), "cox", "sccs", "sccs"),
  design   = c("ps", names(plpTaskList), "cox", "season_calendar", "age_season"),
  needs    = c("cmData.zip", sapply(plpTaskList, `[[`, "folder"), "cmData.zip", "sccsData.zip", "sccsData.zip"),
  stringsAsFactors = FALSE)
datasets <- datasets[datasets$workload %in% workloadsToRun, ]

prepareAll <- function() {
  for (i in seq_len(nrow(datasets))) {
    d <- datasets[i, ]
    if (!file.exists(file.path(dataFolder, d$needs))) {
      say("No ", d$needs, " in ", dataFolder, ": skipping ", d$design)
      next
    }
    for (kind in c("smoke", "full")) {
      file <- file.path(private, sprintf("%s_%s.zip", d$design, kind))
      if (file.exists(paste0(file, ".info.rds"))) next
      fraction <- if (kind == "smoke") smokeFraction else 1
      say("Preparing ", d$design, " (", kind, ")")
      t <- system.time(info <- tryCatch(switch(d$workload,
        ps = preparePs(fraction, file), plp = preparePlp(fraction, file, plpTaskList[[d$design]]), cox = prepareCox(fraction, file),
        sccs = prepareSccs(fraction, file, d$design)),
        error = function(e) { say("  FAILED: ", conditionMessage(e)); NULL }))[["elapsed"]]
      if (!is.null(info)) say(sprintf("  %s rows, %s events, %s covariates, %s non-zeros (%.1f min)",
                                      info$n_rows, info$n_events, info$n_covariates, info$n_nonzero_entries, t / 60))
    }
  }
}

# ---- The fits ----------------------------------------------------------------------------------

# A row per fit: config (lib + fitter), workload, design, variant (what is fitted), threads
fitPlan <- function() {
  lr <- c("ps", "plp")
  configs <- data.frame(config = c("default", "super", "newton", "legacy", "3.7.1"),
                        lib = c("develop", "develop", "newton", "develop", "v371"),
                        fitter = c("default", "super", "super", "legacy", NA), stringsAsFactors = FALSE)
  if (!haveNewton) configs <- configs[configs$config != "newton", ]
  if (!file.exists(file.path(libs[["v371"]], "Cyclops"))) configs <- configs[configs$config != "3.7.1", ]
  configs <- configs[configs$config %in% configsToRun, ]
  configs <- configs[order(match(configs$config, configsToRun)), ]
  rows <- list()
  for (k in seq_len(nrow(configs))) {
    cfg <- configs[k, ]
    for (i in seq_len(nrow(datasets))) {
      d <- datasets[i, ]
      if (cfg$config %in% c("super", "newton") && !(d$workload %in% lr)) next # super = default outside LR
      variants <- switch(d$workload,
        ps = list(c("cv", psThreads)), plp = list(c("cv", 3)),
        cox = list(c("fixed", 1), c("fixed", 8), c("cv", 8)), sccs = list(c("fit", 1)))
      for (v in variants) rows[[length(rows) + 1]] <- data.frame(cfg, workload = d$workload, design = d$design,
                                                                 variant = v[1], threads = as.integer(v[2]))
    }
  }
  plan <- do.call(rbind, rows)
  plan$key <- with(plan, paste(workload, design, variant, threads, config, sep = "/"))
  plan
}

.fitOne <- function(spec) {
  suppressPackageStartupMessages(library(Cyclops))
  null2na <- function(x) if (is.null(x) || length(x) == 0) NA else x
  options(andromedaTempFolder = spec$tempFolder)
  out <- list(version = as.character(utils::packageVersion("Cyclops")),
              sha = tryCatch(utils::packageDescription("Cyclops")$RemoteSha, error = function(e) NA))
  warnings <- character()
  keep <- function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
  ctl <- function(args) { if (!is.na(spec$fitter)) args$fitter <- spec$fitter; do.call(createControl, args) }
  if (spec$workload == "sccs") {
    suppressPackageStartupMessages(library(SelfControlledCaseSeries))
    interval <- loadSccsIntervalData(spec$file)
    on.exit(Andromeda::close(interval))
    control <- ctl(list(cvType = "auto", selectorType = "byPid", startingVariance = 0.1, seed = 1,
                        resetCoefficients = TRUE, noiseLevel = "quiet", threads = spec$threads))
    out$seconds <- system.time(model <- withCallingHandlers(
      fitSccsModel(interval, createFitSccsModelArgs(control = control)), warning = keep))[["elapsed"]]
    out$status <- model$status
    e <- model$estimates
    expo <- e[grepl("^Exposure", e$covariateName), ]
    out$key_estimate <- paste(sprintf("%s IRR %.3f", sub("^Exposure ", "", expo$covariateName), exp(expo$logRr)),
                              collapse = "; ")
    out$key_se_or_ci <- paste(sprintf("%.3f-%.3f", exp(expo$logLb95), exp(expo$logUb95)), collapse = "; ")
    out$beta <- setNames(e$logRr, e$covariateId)
    out$log_lik <- null2na(model$logLikelihood)
  } else {
    a <- Andromeda::loadAndromeda(spec$file)
    on.exit(Andromeda::close(a))
    if (spec$workload == "cox") {
      data <- convertToCyclopsData(a$outcomes, a$covariates, modelType = "cox", quiet = TRUE)
      treatmentId <- spec$treatmentId
      if (spec$variant == "fixed") {
        prior <- createPrior("laplace", variance = 0.01, exclude = treatmentId)
        control <- ctl(list(noiseLevel = "silent", threads = spec$threads))
      } else {
        prior <- createPrior("laplace", exclude = treatmentId, useCrossValidation = TRUE)
        control <- ctl(list(noiseLevel = "silent", cvType = "auto", fold = 10, cvRepetitions = 1, seed = 1,
                            startingVariance = 0.01, threads = spec$threads))
      }
    } else {
      data <- convertToCyclopsData(a$outcomes, a$covariates, modelType = "lr", addIntercept = TRUE, quiet = TRUE)
      prior <- createPrior("laplace", exclude = c(0), useCrossValidation = TRUE)
      control <- if (spec$workload == "ps") {
        ctl(list(noiseLevel = "silent", cvType = "auto", seed = 1, resetCoefficients = TRUE, tolerance = 2e-7,
                 cvRepetitions = 1, startingVariance = 0.01, threads = spec$threads))
      } else {
        ctl(list(cvType = "auto", fold = 3, startingVariance = 0.01, lowerLimit = 0.01, upperLimit = 20,
                 tolerance = 2e-6, cvRepetitions = 1, selectorType = "byPid", noiseLevel = "silent",
                 threads = spec$threads, maxIterations = 3000, seed = 42))
      }
    }
    out$seconds <- system.time(fit <- withCallingHandlers(
      fitCyclopsModel(data, prior = prior, control = control), warning = keep))[["elapsed"]]
    beta <- coef(fit, ignoreConvergence = TRUE)
    out$status <- fit$return_flag
    out$variance <- fit$variance[1]
    out$iterations <- null2na(fit$iterations)
    out$log_lik <- fit$log_likelihood
    out$n_nonzero_coefs <- sum(beta[names(beta) != "(Intercept)"] != 0)
    out$beta <- beta
    if (spec$workload == "cox") {
      b <- beta[as.character(treatmentId)]
      out$key_estimate <- sprintf("HR %.4f", exp(b))
      ci <- tryCatch(confint(fit, parm = treatmentId, includePenalty = TRUE), error = function(e) NULL)
      out$key_se_or_ci <- if (is.null(ci)) NA else sprintf("%.4f-%.4f", exp(ci[1, 2]), exp(ci[1, 3]))
    }
  }
  out$warnings <- paste(unique(warnings), collapse = " | ")
  out
}

runFit <- function(p, kind, rep, timeout) {
  file <- normalizePath(file.path(private, sprintf("%s_%s.zip", p$design, kind)))
  info <- readRDS(paste0(file, ".info.rds"))
  spec <- c(as.list(p[, c("workload", "design", "variant", "threads", "fitter")]),
            list(file = file, treatmentId = null2na(info$treatmentId), tempFolder = tempFolder))
  res <- tryCatch(
    callr::r(.fitOne, args = list(spec = spec), libpath = c(libs[[p$lib]], .libPaths()), timeout = timeout,
             user_profile = FALSE), # an renv project's .Rprofile would reset the library path
    error = function(e) list(status = if (inherits(e, "callr_timeout_error")) "TIMEOUT" else
      paste("ERROR:", gsub("[\r\n,]+", " ", conditionMessage(e)))))
  if (!is.null(res$version) && (p$config == "3.7.1") != (res$version == "3.7.1")) {
    stop("Config ", p$config, " loaded Cyclops ", res$version, " from ", libs[[p$lib]])
  }
  if (!is.null(res$beta) && kind == "full") {
    saveRDS(res$beta, file.path(private, "coefs", sprintf("%s_%d.rds", gsub("/", "_", p$key), rep)))
  }
  data.frame(env = databaseLabel, workload = p$workload, design = paste0(p$design, if (p$variant != "cv" &&
             p$workload == "cox") paste0("_", p$variant) else ""),
             config = p$config, cyclops_version = null2na(res$version), cyclops_sha = null2na(res$sha),
             fitter = p$fitter, threads = p$threads, rep = rep, seconds = null2na(res$seconds),
             status = null2na(res$status), iterations = null2na(res$iterations), variance = null2na(res$variance),
             n_rows = info$n_rows, n_events = info$n_events, n_covariates = info$n_covariates,
             n_nonzero_entries = info$n_nonzero_entries, n_nonzero_coefs = null2na(res$n_nonzero_coefs),
             log_lik = null2na(res$log_lik), key_estimate = null2na(res$key_estimate),
             key_se_or_ci = null2na(res$key_se_or_ci), warnings = substr(null2na(res$warnings), 1, 300),
             timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), key = p$key, kind = kind,
             stringsAsFactors = FALSE)
}

appendRow <- function(row, file) {
  write.table(row, file, sep = ",", row.names = FALSE, col.names = !file.exists(file), append = file.exists(file),
              qmethod = "double")
}
readRows <- function(file) if (file.exists(file)) read.csv(file, stringsAsFactors = FALSE) else NULL

# Expected seconds of a full fit: measured, else the smoke ratio to default times the full default,
# else the smoke time times the median full/smoke ratio so far (20x before any full fit)
estimate <- function(key, smoke, full) {
  if (!is.null(full) && any(full$key == key & !is.na(full$seconds))) return(median(full$seconds[full$key == key], na.rm = TRUE))
  s <- smoke$seconds[smoke$key == key][1]
  if (is.null(s) || is.na(s)) s <- 1
  parts <- strsplit(key, "/")[[1]]
  defaultKey <- paste(c(parts[1:4], "default"), collapse = "/")
  sd <- smoke$seconds[smoke$key == defaultKey][1]
  if (!is.null(full) && any(full$key == defaultKey & !is.na(full$seconds)) && !is.na(sd) && sd > 0) {
    return(median(full$seconds[full$key == defaultKey], na.rm = TRUE) * s / sd)
  }
  ratio <- 20
  if (!is.null(full)) {
    both <- merge(aggregate(seconds ~ key, full, median), aggregate(seconds ~ key, smoke, median), by = "key")
    if (nrow(both) > 0) ratio <- median(both$seconds.x / pmax(both$seconds.y, 0.5))
  }
  s * ratio
}

# Each configuration's library must hold the right Cyclops before anything runs
checkLibraries <- function(plan) {
  for (cfg in unique(plan$config)) {
    lib <- libs[[plan$lib[plan$config == cfg][1]]]
    v <- tryCatch(as.character(utils::packageVersion("Cyclops", lib.loc = lib)), error = function(e) NA)
    ok <- !is.na(v) && if (cfg == "3.7.1") v == "3.7.1" else utils::compareVersion(v, "4.0.0") >= 0
    if (!ok) stop(sprintf("Config %s needs Cyclops %s in %s, found %s. Run installCyclops() (see README).",
                          cfg, if (cfg == "3.7.1") "3.7.1" else ">= 4.0", lib, if (is.na(v)) "none" else v))
    say(sprintf("Config %-8s Cyclops %s from %s", cfg, v, lib))
  }
}

runAll <- function() {
  plan <- fitPlan()
  plan <- plan[file.exists(file.path(private, sprintf("%s_full.zip.info.rds", plan$design))), ]
  say(nrow(plan), " fits planned: ", paste(unique(plan$config), collapse = ", "))
  checkLibraries(plan)

  # Smoke test: every fit once on the subsample
  smoke <- readRows(smokeFile)
  for (i in seq_len(nrow(plan))) {
    p <- plan[i, ]
    if (!is.null(smoke) && p$key %in% smoke$key[!grepl("^ERROR", smoke$status)]) next # errors are retried
    row <- runFit(p, "smoke", 1, timeoutMinutes * 60)
    appendRow(row, smokeFile)
    say(sprintf("smoke %-40s %7.1f s  %s", p$key, row$seconds, row$status))
  }
  smoke <- readRows(smokeFile)
  if (!is.null(smoke)) smoke <- smoke[!duplicated(smoke$key, fromLast = TRUE), ]
  if (is.null(smoke)) { say("Nothing to run"); return(invisible(NULL)) }
  bad <- smoke[!(smoke$status %in% c("SUCCESS", "OK", "NO_CONVERGENCE")) , c("key", "status")]
  if (nrow(bad) > 0) { say("Smoke-test problems (these fits still run):"); for (j in seq_len(nrow(bad))) say("  ", bad$key[j], ": ", bad$status[j]) }

  # Full runs, fast configurations first (plan order), then extra repetitions of short fits
  start <- Sys.time()
  budget <- budgetHours * 3600
  passes <- c(list(list(rep = 1, keys = plan$key)), lapply(seq_len(maxReps - 1) + 1, function(r) list(rep = r, keys = NULL)))
  for (pass in passes) {
    full <- readRows(resultsFile)
    keys <- if (is.null(pass$keys)) {
      if (is.null(full)) character() else intersect(plan$key, full$key[full$rep == 1 & !is.na(full$seconds) &
                                                    full$seconds < 300 & full$status %in% c("SUCCESS", "OK")])
    } else pass$keys
    for (key in keys) {
      full <- readRows(resultsFile)
      # errors and budget skips are tried again; timeouts too with BENCH_RETRY_TIMEOUT=1 (and a longer timeout)
      retry <- grepl("^ERROR", full$status) | full$status == "SKIPPED_BUDGET" |
        (full$status == "TIMEOUT" & nzchar(Sys.getenv("BENCH_RETRY_TIMEOUT")))
      if (!is.null(full) && any(full$key == key & full$rep == pass$rep & !retry)) next
      p <- plan[plan$key == key, ]
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      expected <- estimate(key, smoke, full)
      remaining <- budget - elapsed
      if (expected > remaining) {
        if (pass$rep == 1) {
          say(sprintf("skip  %-40s expected %s > %s left", key, hms(expected), hms(remaining)))
          appendRow(data.frame(env = databaseLabel, workload = p$workload, design = p$design, config = p$config,
            cyclops_version = NA, cyclops_sha = NA, fitter = p$fitter, threads = p$threads, rep = 1, seconds = NA,
            status = "SKIPPED_BUDGET", iterations = NA, variance = NA, n_rows = NA, n_events = NA, n_covariates = NA,
            n_nonzero_entries = NA, n_nonzero_coefs = NA, log_lik = NA, key_estimate = NA, key_se_or_ci = NA,
            warnings = sprintf("expected %.0f s", expected), timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            key = key, kind = "full"), resultsFile)
        }
        next
      }
      row <- runFit(p, "full", pass$rep, min(timeoutMinutes * 60, max(60, remaining)))
      appendRow(row, resultsFile)
      full <- readRows(resultsFile)
      todo <- setdiff(plan$key, full$key[full$rep == 1])
      eta <- sum(vapply(todo, estimate, numeric(1), smoke = smoke, full = full))
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      say(sprintf("%-46s rep %d %8.1f s  %-10s | %d/%d first runs done, elapsed %s, ETA %s (~%s, budget ends %s)",
                  key, pass$rep, row$seconds, row$status, length(intersect(plan$key, full$key[full$rep == 1])), nrow(plan),
                  hms(elapsed), hms(min(eta, budget - elapsed)), format(Sys.time() + min(eta, budget - elapsed), "%H:%M"),
                  format(start + budget, "%H:%M")))
    }
  }
}

# ---- Summary and report ------------------------------------------------------------------------

report <- function() {
  results <- readRows(resultsFile)
  if (is.null(results)) return(invisible(NULL))
  results <- results[results$config %in% configsToRun, ]
  coefDiff <- function(key, config) {
    ref <- file.path(private, "coefs", sprintf("%s_1.rds", gsub("/", "_", sub("/[^/]+$", "/default", key))))
    this <- file.path(private, "coefs", sprintf("%s_1.rds", gsub("/", "_", key)))
    if (config == "default" || !file.exists(ref) || !file.exists(this)) return(NA)
    a <- readRDS(ref); b <- readRDS(this); ids <- union(names(a), names(b))
    a <- a[ids]; b <- b[ids]; a[is.na(a)] <- 0; b[is.na(b)] <- 0
    max(abs(a - b))
  }
  summary <- results %>%
    group_by(env, workload, design, config, threads, key) %>%
    summarise(runs = sum(!is.na(seconds)), seconds = median(seconds, na.rm = TRUE),
              status = paste(unique(status), collapse = "/"), iterations = median(iterations, na.rm = TRUE),
              variance = median(variance, na.rm = TRUE), n_rows = suppressWarnings(max(n_rows, na.rm = TRUE)),
              n_events = suppressWarnings(max(n_events, na.rm = TRUE)),
              n_covariates = suppressWarnings(max(n_covariates, na.rm = TRUE)),
              n_nonzero_coefs = median(n_nonzero_coefs, na.rm = TRUE), log_lik = median(log_lik, na.rm = TRUE),
              key_estimate = first(key_estimate), key_se_or_ci = first(key_se_or_ci), .groups = "drop") %>%
    group_by(env, workload, design, threads) %>%
    mutate(speedup_vs_371 = seconds[config == "3.7.1"][1] / seconds,
           speedup_vs_legacy = seconds[config == "legacy"][1] / seconds,
           speedup_vs_default = seconds[config == "default"][1] / seconds) %>%
    ungroup() %>%
    mutate(max_abs_coef_diff_vs_default = mapply(coefDiff, key, config)) %>%
    arrange(workload, design, threads, match(config, c("3.7.1", "legacy", "default", "super", "newton"))) %>%
    select(-key)
  summary <- as.data.frame(summary)
  summary[] <- lapply(summary, function(x) if (is.numeric(x)) ifelse(is.finite(x), x, NA) else x)
  write.csv(summary, file.path(folder, "summary.csv"), row.names = FALSE)
  write.csv(results[, setdiff(names(results), c("key", "kind"))], file.path(folder, "results_clean.csv"), row.names = FALSE)

  md <- c(sprintf("# Cyclops 4.0 benchmark: %s", databaseLabel), "",
          sprintf("Cyclops develop %s; %s", paste(unique(na.omit(results$cyclops_sha)), collapse = ", "), .machine()), "",
          "Median seconds of the Cyclops call(s) only; each fit in a fresh R process. Speed-up relative to Cyclops 3.7.1.", "")
  for (w in unique(summary$workload)) {
    s <- summary[summary$workload == w, ]
    md <- c(md, sprintf("## %s", w), "", sprintf("%s rows, %s events, %s covariates", s$n_rows[1], s$n_events[1], s$n_covariates[1]), "",
            "| design | threads | config | seconds | runs | status | speed-up vs 3.7.1 | variance | non-zero | estimate | CI | max coef diff vs default |",
            "|---|---:|---|---:|---:|---|---:|---:|---:|---|---|---:|",
            sprintf("| %s | %d | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s |", s$design, s$threads, s$config,
                    ifelse(is.na(s$seconds), "-", sprintf("%.1f", s$seconds)), s$runs, s$status,
                    ifelse(is.na(s$speedup_vs_371), "-", sprintf("%.1fx", s$speedup_vs_371)),
                    ifelse(is.na(s$variance), "-", signif(s$variance, 3)), ifelse(is.na(s$n_nonzero_coefs), "-", s$n_nonzero_coefs),
                    ifelse(is.na(s$key_estimate), "", s$key_estimate), ifelse(is.na(s$key_se_or_ci), "", s$key_se_or_ci),
                    ifelse(is.na(s$max_abs_coef_diff_vs_default), "-", signif(s$max_abs_coef_diff_vs_default, 2))), "")
  }
  failed <- results[!(results$status %in% c("SUCCESS", "OK")), c("workload", "design", "config", "status", "warnings")]
  if (nrow(failed) > 0) md <- c(md, "## Failed, timed out or skipped", "",
                                sprintf("- %s / %s / %s: %s %s", failed$workload, failed$design, failed$config, failed$status,
                                        ifelse(is.na(failed$warnings), "", failed$warnings)), "")
  writeLines(md, file.path(folder, "REPORT.md"))
  say("Wrote summary.csv and REPORT.md")
  invisible(summary)
}

.machine <- function() {
  cpu <- tryCatch({
    os <- Sys.info()[["sysname"]]
    if (os == "Linux") sub(".*: ", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1])
    else if (os == "Darwin") system("sysctl -n machdep.cpu.brand_string", intern = TRUE)
    else Sys.getenv("PROCESSOR_IDENTIFIER")
  }, error = function(e) NA_character_)
  ram <- tryCatch(if (.Platform$OS.type == "windows") NA else
    sprintf("%.0f GB RAM", as.numeric(sub("\\D+(\\d+).*", "\\1", grep("MemTotal", readLines("/proc/meminfo"), value = TRUE))) / 2^20),
    error = function(e) NA)
  sprintf("%s, %d cores, %s, %s, %s", cpu, parallel::detectCores(), ram, utils::sessionInfo()$running, R.version.string)
}

writeEnvironment <- function() {
  pk <- c("Cyclops", "callr", "Andromeda", "DatabaseConnector", "FeatureExtraction", "CohortMethod",
          "SelfControlledCaseSeries", "PatientLevelPrediction", "survival")
  versions <- sapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "-"))
  cy <- sapply(libs, function(l) tryCatch(paste(packageDescription("Cyclops", lib.loc = l)[c("Version", "RemoteSha")], collapse = " "),
                                          error = function(e) "-"))
  writeLines(c(paste("Machine:", .machine()), paste("Libraries:", paste(names(cy), cy, sep = " = ", collapse = "; ")),
               paste("Packages:", paste(pk, versions, sep = " ", collapse = ", ")),
               sprintf("Settings: budget %s h, timeout %s min, smoke fraction %s, max reps %s", budgetHours, timeoutMinutes,
                       smokeFraction, maxReps),
               "Deviations from the slide instructions: cohorts from extract_cdm.R (plain SQL);",
               "COX outcome is death where usable, else AMI; COXC (left-censored Cox) not run; super only for LR."),
             file.path(folder, "environment.txt"))
}

# ---- Main --------------------------------------------------------------------------------------

args <- commandArgs(TRUE)
if (!interactive()) {
  if (length(args) == 0 || args[1] == "all") {
    writeEnvironment()
    prepareAll()
    runAll()
    report()
  } else if (args[1] == "report") {
    report()
  } else if (args[1] == "prepare") {
    prepareAll()
  }
}
