# CyclopsBenchmarks

Speed of Cyclops 4.0 (`fitter = "legacy"`, `"default"`, `"super"`, and optionally the active-set Newton branch) on
HADES workloads, timed on real OMOP CDM data. Only aggregate results leave the server.

| step | script | needs a database | output (stays on the server) |
|---|---|---|---|
| 1. cohorts + extraction | `extract_cdm.R` | yes | `cyclops_benchmark_data/<db>/`: cmData.zip, plpData/, sccsData.zip, cohort_counts.csv |
| 2. benchmark | `run_benchmark.R` | no | `cyclops_benchmark_results/<db>/`: results.csv, summary.csv, REPORT.md, environment.txt, progress.txt, private/ |

## Workloads

- **ps**: propensity model as `CohortMethod::createPs()` fits it: new users of lisinopril vs hydrochlorothiazide,
  default FeatureExtraction covariates, Laplace prior, 10-fold auto CV, 8 threads.
- **plp**: LASSO as PatientLevelPrediction fits it: type 2 diabetes -> heart failure in 1-1095 days, default
  covariates, 75% training split, 3-fold auto CV, 3 threads.
- **cox**: Cox outcome model on the PS cohort (death within 5 years where death data are usable, else AMI):
  variance 0.01 fit at 1 and 8 threads, and 10-fold CV at 8 threads.
- **sccs**: `SelfControlledCaseSeries::fitSccsModel()` for AMI with ibuprofen and naproxen, two designs
  (season + calendar time; age + season).

`super` only differs from `default` for logistic regression, so it runs for ps and plp only.

## Run

```sh
git clone https://github.com/egillax/CyclopsBenchmarks.git
cd CyclopsBenchmarks
```

1. **Extract** (edit the Settings block first: CDM schema, a writable scratch schema, connection; Databricks by
   default with the `DATABASE_*` environment variables as in SignalAdaptiveStudy):
   ```sh
   Rscript extract_cdm.R
   ```
2. **Install Cyclops develop** into a separate library (your R library is not touched):
   ```r
   source("run_benchmark.R")   # in an interactive session this only defines the functions
   installCyclops()            # develop from GitHub into lib/cyclops-develop
   ```
3. **Benchmark** (fresh R process per fit; smoke test on a 5% subsample first; then fastest configuration first,
   with an ETA after each fit; 4-hour budget and 60-minute per-fit timeout by default):
   ```sh
   BENCH_LABEL=CCAE Rscript run_benchmark.R
   ```
   Settings can also be given as environment variables: `BENCH_LABEL`, `BENCH_DATA`, `BENCH_OUT`,
   `BENCH_LIB_DEVELOP`, `BENCH_LIB_NEWTON`, `BENCH_LIB_HADES`, `BENCH_BUDGET_HOURS`, `BENCH_TIMEOUT_MINUTES`,
   `BENCH_WORKLOADS` (e.g. `ps,plp`), `BENCH_PLP_TASKS` (extra PLP tasks, see the Settings block), `BENCH_TEMP` (Andromeda and temp files; default `./tmp`, also used by `extract_cdm.R`).
   Rerunning resumes: finished fits are skipped. `Rscript run_benchmark.R report` rebuilds the summary only.

Share `summary.csv`, `results_clean.csv`, `environment.txt` and `REPORT.md` once you've checked that your data
licence allows it. Never share `private/`.
