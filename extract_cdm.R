# Cyclops 4.0 benchmark: cohort definitions and data extraction
# ==============================================================
#
# Extracts, once per database, the HADES data objects the benchmark fits on:
#
#   cmData     CohortMethodData: new users of lisinopril (target) vs hydrochlorothiazide
#              (comparator), default FeatureExtraction covariates (both drugs excluded),
#              outcomes AMI and death                 -> LR-PS and COX workloads
#   plpData    PatientLevelPrediction data: first type 2 diabetes diagnosis (target),
#              heart failure (outcome), default covariates                -> LR-PLP workload
#   sccsData   SelfControlledCaseSeries data: first AMI (outcome), ibuprofen and naproxen
#              (exposures, drug_era), 2010-2019                           -> SCCS workload
#
# The cohorts are plain OHDSI SQL (below), translated by SqlRender to your dbms, written to
# <cohortDatabaseSchema>.<cohortTable>. Nothing in the CDM schema is modified.
#
# Output: <outputFolder>/<databaseLabel>/ with cmData.zip, plpData/, sccsData.zip,
# cohort_counts.csv (aggregate) and extract.log. These stay on the server; the benchmark reads
# them without a database connection.
#
# Every step is skipped when its output exists, so a failed run can be restarted as is.
#
# Needs: DatabaseConnector, SqlRender, FeatureExtraction, CohortMethod (>= 5),
# PatientLevelPrediction (>= 6), SelfControlledCaseSeries (>= 5), Andromeda.

# ---- Settings (edit) ---------------------------------------------------------------------------

# The Databricks JDBC driver (3.x) reads results with Apache Arrow, which on Java >= 16 needs this permission.
# It must be set before Java starts: run this script with Rscript, or in a fresh R session.
options(java.parameters = "--add-opens=java.base/java.nio=org.apache.arrow.memory.core,ALL-UNNAMED")

databaseLabel        <- "CCAE"                  # name on the slides
outputFolder         <- "cyclops_benchmark_data"
cdmDatabaseSchema    <- "cdm_ccae"              # CDM schema (read only)
cohortDatabaseSchema <- "scratch"               # a schema you may write to (e.g. catalog.schema)
cohortTable          <- "cyclops_bench_cohort"  # created (or replaced) there
tempEmulationSchema  <- cohortDatabaseSchema    # Spark has no temp tables: they go here (NULL on platforms with them)
options(sqlRenderTempEmulationSchema = tempEmulationSchema)

# Databricks, as in SignalAdaptiveStudy: credentials from environment variables
#   DATABASE_SERVER    Databricks host, or a full jdbc:databricks://... connection string
#   DATABASE_PATH      the SQL warehouse's HTTP path (not needed with a full string)
#   DATABASE_PASSWORD  a personal access token (user "token")
#   DATABASE_DRIVER    folder with the JDBC driver (DatabaseConnector::downloadJdbcDrivers("spark", pathToDriver = ...))
server <- Sys.getenv("DATABASE_SERVER")
connectionString <- if (startsWith(server, "jdbc:")) server else
  sprintf("jdbc:databricks://%s:443/default;transportMode=http;ssl=1;AuthMech=3;httpPath=%s;", server,
          Sys.getenv("DATABASE_PATH"))
if (!grepl("EnableArrow=", connectionString, ignore.case = TRUE)) {
  connectionString <- paste0(sub(";?$", ";", connectionString), "EnableArrow=1;")
}
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "spark",
  connectionString = connectionString,
  user = "token",
  password = Sys.getenv("DATABASE_PASSWORD"),
  pathToDriver = Sys.getenv("DATABASE_DRIVER")
)
# DuckDB instead (set tempEmulationSchema <- NULL and maxCohortSize <- 0):
# connectionDetails <- DatabaseConnector::createDbiConnectionDetails(dbms = "duckdb", drv = duckdb::duckdb(),
#                                                                    dbdir = "path/to/cdm.duckdb")

# Size caps (0 = no cap). With a full claims database the uncapped cohorts are millions of people
# and default covariates for all of them would take very long to extract.
maxCohortSize   <- 250000   # CohortMethod: per treatment arm (createPs() fits at most 250,000 per arm)
                            # On DuckDB use 0: CohortMethod's downsampling SQL fails there
plpSampleSize   <- 250000   # PLP: random sample of the target cohort
maxSccsCases    <- 200000   # SCCS: cases of the outcome

steps <- c("cohorts", "cmData", "plpData", "sccsData")   # drop the ones you don't want

# ---- Cohort definitions ------------------------------------------------------------------------
#
# id  name                         definition
#  1  lisinopril new users         first exposure to any drug with ingredient lisinopril (1308216),
#                                  >= 365 days of prior observation; ends at end of observation
#  2  hydrochlorothiazide new users  same, ingredient 974166
#  3  acute myocardial infarction  first ever condition occurrence of 4329847 or a descendant,
#                                  within an observation period
#  4  death                        death date, within (or 1 day after) an observation period
#  5  type 2 diabetes              first ever condition occurrence of 201826 or a descendant,
#                                  >= 365 days of prior observation, age >= 18
#  6  heart failure                first ever condition occurrence of 316139 or a descendant
#
# Exposures for SCCS come from drug_era (ingredients ibuprofen 1177480, naproxen 1115008).
#
# Study designs applied on top (in the benchmark, not here):
#   PS/COX   CohortMethod: first exposure only, 365-day washout, people in both cohorts kept in the
#            first one (truncated at the second), restricted to the common period.
#            COX outcome: AMI (or death where death data are usable) in the 5 years after index.
#   PLP      outcome heart failure 1-1095 days after index, prior outcome removed.
#   SCCS     first AMI, ibuprofen and naproxen eras with a 30-day pre-exposure window, 2010-2019.

cohortSql <- "
IF OBJECT_ID('@cohort_database_schema.@cohort_table', 'U') IS NOT NULL
  DROP TABLE @cohort_database_schema.@cohort_table;

CREATE TABLE @cohort_database_schema.@cohort_table (
  cohort_definition_id BIGINT,
  subject_id BIGINT,
  cohort_start_date DATE,
  cohort_end_date DATE
);

-- 1, 2: new users of an ingredient
INSERT INTO @cohort_database_schema.@cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
SELECT CASE WHEN f.ingredient = 1308216 THEN 1 ELSE 2 END, f.person_id, f.start_date, op.observation_period_end_date
FROM (
  SELECT ca.ancestor_concept_id AS ingredient, de.person_id, MIN(de.drug_exposure_start_date) AS start_date
  FROM @cdm_database_schema.drug_exposure de
  INNER JOIN @cdm_database_schema.concept_ancestor ca ON de.drug_concept_id = ca.descendant_concept_id
  WHERE ca.ancestor_concept_id IN (1308216, 974166)
  GROUP BY ca.ancestor_concept_id, de.person_id
) f
INNER JOIN @cdm_database_schema.observation_period op
  ON f.person_id = op.person_id
 AND f.start_date >= op.observation_period_start_date
 AND f.start_date <= op.observation_period_end_date
WHERE DATEADD(DAY, 365, op.observation_period_start_date) <= f.start_date;

-- 3, 5, 6: first ever diagnosis
INSERT INTO @cohort_database_schema.@cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
SELECT f.cohort_id, f.person_id, f.start_date, op.observation_period_end_date
FROM (
  SELECT CASE WHEN ca.ancestor_concept_id = 4329847 THEN 3
              WHEN ca.ancestor_concept_id = 201826 THEN 5 ELSE 6 END AS cohort_id,
         co.person_id, MIN(co.condition_start_date) AS start_date
  FROM @cdm_database_schema.condition_occurrence co
  INNER JOIN @cdm_database_schema.concept_ancestor ca ON co.condition_concept_id = ca.descendant_concept_id
  WHERE ca.ancestor_concept_id IN (4329847, 201826, 316139)
  GROUP BY ca.ancestor_concept_id, co.person_id
) f
INNER JOIN @cdm_database_schema.observation_period op
  ON f.person_id = op.person_id
 AND f.start_date >= op.observation_period_start_date
 AND f.start_date <= op.observation_period_end_date
INNER JOIN @cdm_database_schema.person p ON f.person_id = p.person_id
WHERE f.cohort_id <> 5
   OR (DATEADD(DAY, 365, op.observation_period_start_date) <= f.start_date
       AND YEAR(f.start_date) - p.year_of_birth >= 18);

-- 4: death
INSERT INTO @cohort_database_schema.@cohort_table (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
SELECT 4, d.person_id, MIN(d.death_date), MIN(d.death_date)
FROM @cdm_database_schema.death d
INNER JOIN @cdm_database_schema.observation_period op
  ON d.person_id = op.person_id
 AND d.death_date >= op.observation_period_start_date
 AND d.death_date <= DATEADD(DAY, 1, op.observation_period_end_date)
GROUP BY d.person_id;
"

# ---- Run ---------------------------------------------------------------------------------------

folder <- file.path(outputFolder, databaseLabel)
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
logFile <- file.path(folder, "extract.log")
say <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  message(line)
  cat(line, "\n", file = logFile, append = TRUE)
}
timed <- function(label, expr) {
  say(label, " ...")
  t <- system.time(value <- expr)[["elapsed"]]
  say(label, sprintf(" done in %.1f min", t / 60))
  invisible(value)
}

say("Packages: ", paste(sapply(c("DatabaseConnector", "SqlRender", "FeatureExtraction", "CohortMethod",
                                 "PatientLevelPrediction", "SelfControlledCaseSeries", "Andromeda"),
                               function(p) paste(p, tryCatch(as.character(packageVersion(p)), error = function(e) "-"))),
                        collapse = ", "))

if ("cohorts" %in% steps && !file.exists(file.path(folder, "cohort_counts.csv"))) {
  timed("Creating cohorts", {
    connection <- DatabaseConnector::connect(connectionDetails)
    DatabaseConnector::renderTranslateExecuteSql(connection, cohortSql,
                                                 cdm_database_schema = cdmDatabaseSchema,
                                                 cohort_database_schema = cohortDatabaseSchema,
                                                 cohort_table = cohortTable)
    counts <- DatabaseConnector::renderTranslateQuerySql(connection,
      "SELECT cohort_definition_id, COUNT(*) AS entries, COUNT(DISTINCT subject_id) AS persons
       FROM @cohort_database_schema.@cohort_table GROUP BY cohort_definition_id ORDER BY cohort_definition_id;",
      cohort_database_schema = cohortDatabaseSchema, cohort_table = cohortTable, snakeCaseToCamelCase = TRUE)
    DatabaseConnector::disconnect(connection)
    counts$name <- c("1" = "lisinopril new users", "2" = "hydrochlorothiazide new users", "3" = "AMI",
                     "4" = "death", "5" = "type 2 diabetes", "6" = "heart failure")[as.character(counts$cohortDefinitionId)]
    write.csv(counts, file.path(folder, "cohort_counts.csv"), row.names = FALSE)
    print(counts)
  })
}

if ("cmData" %in% steps && !file.exists(file.path(folder, "cmData.zip"))) {
  timed("Extracting CohortMethodData", {
    covariateSettings <- FeatureExtraction::createDefaultCovariateSettings(
      excludedCovariateConceptIds = c(1308216, 974166), addDescendantsToExclude = TRUE)
    cmData <- CohortMethod::getDbCohortMethodData(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = cdmDatabaseSchema,
      targetId = 1, comparatorId = 2, outcomeIds = c(3, 4),
      exposureDatabaseSchema = cohortDatabaseSchema, exposureTable = cohortTable,
      outcomeDatabaseSchema = cohortDatabaseSchema, outcomeTable = cohortTable,
      getDbCohortMethodDataArgs = CohortMethod::createGetDbCohortMethodDataArgs(
        removeDuplicateSubjects = "keep first, truncate to second",
        firstExposureOnly = TRUE, washoutPeriod = 365, restrictToCommonPeriod = TRUE,
        maxCohortSize = maxCohortSize, covariateSettings = covariateSettings))
    CohortMethod::saveCohortMethodData(cmData, file.path(folder, "cmData.zip"))
  })
}

if ("plpData" %in% steps && !dir.exists(file.path(folder, "plpData"))) {
  timed("Extracting plpData", {
    databaseDetails <- PatientLevelPrediction::createDatabaseDetails(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = cdmDatabaseSchema, cdmDatabaseName = databaseLabel, cdmDatabaseId = databaseLabel,
      tempEmulationSchema = if (is.null(tempEmulationSchema)) cdmDatabaseSchema else tempEmulationSchema,
      cohortDatabaseSchema = cohortDatabaseSchema, cohortTable = cohortTable,
      outcomeDatabaseSchema = cohortDatabaseSchema, outcomeTable = cohortTable,
      targetId = 5, outcomeIds = 6)
    plpData <- PatientLevelPrediction::getPlpData(
      databaseDetails = databaseDetails,
      covariateSettings = FeatureExtraction::createDefaultCovariateSettings(),
      restrictPlpDataSettings = PatientLevelPrediction::createRestrictPlpDataSettings(
        firstExposureOnly = TRUE, washoutPeriod = 365,
        sampleSize = if (plpSampleSize > 0) plpSampleSize else NULL))
    PatientLevelPrediction::savePlpData(plpData, file.path(folder, "plpData"))
  })
}

if ("sccsData" %in% steps && !file.exists(file.path(folder, "sccsData.zip"))) {
  timed("Extracting SccsData", {
    sccsData <- SelfControlledCaseSeries::getDbSccsData(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = cdmDatabaseSchema,
      outcomeDatabaseSchema = cohortDatabaseSchema, outcomeTable = cohortTable, outcomeIds = 3,
      exposureDatabaseSchema = cdmDatabaseSchema, exposureTable = "drug_era",
      getDbSccsDataArgs = SelfControlledCaseSeries::createGetDbSccsDataArgs(
        exposureIds = c(1177480, 1115008),
        studyStartDates = "20100101", studyEndDates = "20191231",
        maxCasesPerOutcome = maxSccsCases))
    SelfControlledCaseSeries::saveSccsData(sccsData, file.path(folder, "sccsData.zip"))
  })
}

say("Finished. Output in ", normalizePath(folder))
