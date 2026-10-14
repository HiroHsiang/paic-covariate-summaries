# Choosing covariate summaries for population-adjusted indirect comparisons

R code and summary results for the simulation study

> Xiang Y, Zhang C, Zhang K, Lian D, Sun H, Yu Y, Chen Y. *Choosing Covariate Summaries for Population-Adjusted Indirect Comparisons: A Simulation Study.* (manuscript submitted)

Population-adjusted indirect comparisons (MAIC, STC, G-computation and G-MAIC) obtain all information on the comparator population from the summary statistics reported in its publication. Most analyses match the means and standard deviations (SDs) of the effect modifiers. The study asks whether additional summaries (quartiles, ranges and proportions above a clinical cut-point) improve the comparison, and when they make it worse. It also evaluates two proposed methods that match SDs and quartiles within tolerances set by their sampling error: soft-constraint MAIC (sMAIC) and soft-constraint G-MAIC (sG-MAIC).

## Contents

| Folder | Content |
|---|---|
| `R/` | Functions: data-generating mechanisms (`01`), outcome models and true values (`02`), balancing weights with exact and soft constraints (`03`), population-limit calculations (`04`), data generation and published tables (`05`), MAIC family, sMAIC and variance estimators (`06`), reconstruction of marginal distributions, STC, G-computation, G-MAIC and sG-MAIC (`07`), analysis of one replicate and the parallel driver (`08`) |
| `sims/` | Numbered scripts that run the simulation and produce the tables and figures (see *Reproducing the study*) |
| `examples/` | `apply_to_published_table.R`: one simulated submission analysed with the nine main-text methods |
| `output/settings/`, `output/poplimit/` | Scenario settings and calibrated design parameters (`dgm_params.rds`), read by every run script |
| `output/summary/` | Performance of every method in every scenario (bias, empirical SE, RMSE, coverage and ESS, each with its Monte Carlo SE) |
| `output/figures/` | Figures 1–4 and Supplementary Table S1 |

The per-replicate estimates (about 450 MB) are not included. The run scripts regenerate them.

## Requirements

R 4.6.0 or later. Install the packages with

```r
source("R/00_packages.R"); ensure_packages()
```

The code uses survival, sn, estmeansd, numDeriv, foreach, doParallel and doRNG; rsimsum for the check of the performance measures; and ggplot2, patchwork and ragg for the figures. `SESSION_INFO.txt` lists the versions used. Run all scripts from the repository root.

## Quick start

```sh
Rscript examples/apply_to_published_table.R
```

The example generates one submission: IPD from an A-versus-C trial and the published B-versus-C log odds ratio with a baseline table. It then estimates the marginal A-versus-B log odds ratio with MAIC (means), MAIC (means and SDs), MAIC (all reported summaries), sMAIC, STC, G-computation, G-MAIC (means and SDs), G-MAIC and sG-MAIC. It takes about one minute. The comments in the script show how to enter a published baseline table by hand. The current code handles three continuous covariates named `x1`–`x3`, as in the simulation design.

## Reproducing the study

The scripts take two optional arguments, the number of replicates and the number of cores (for example `Rscript sims/10_main_run.R 1000 12`). Run times below are for 12 cores. All run scripts can be resumed: a scenario whose output file exists is skipped.

| Order | Script | Output | Run time |
|---|---|---|---|
| 1 | `sims/00_settings_v1.2.R` | `output/settings/` | seconds |
| 2 | `sims/01_stage1_dgp_check.R` | `output/poplimit/dgm_params.rds` (calibration of effect-modification strength) | minutes |
| 3 | `sims/10_main_run.R` | `output/main/` (96 scenarios, 1000 replicates, true values) | 13 h |
| 4 | `sims/11_coverage_subset.R` | `output/coverage/` (bootstrap intervals, 24 anchored binary scenarios, 300 replicates) | 2 h |
| 5 | `sims/15_supp_m0_softgmaic.R` | `output/supp/` (unadjusted comparison and sG-MAIC) | 2 h |
| 6 | `sims/17_patch_reconstruct.R`, `sims/19_patch_coverage.R` | `output/main_patch/`, `output/coverage_patch/` (methods that reconstruct marginal distributions) | 3 h |
| 7 | `sims/22_extra_m1q.R` | `output/extra_m1q/` (MAIC with means, quartiles and proportions) | 0.5 h |
| 8 | `sims/32_patch_stc_unanch_tte.R` | `output/main_patch_stc/` (STC for unanchored time-to-event outcomes) | 0.2 h |
| 9 | `sims/41_patch_rmst_se.R` | `output/patch_rmst_se/` (infinitesimal-jackknife SE of the weighted RMST) | 2 h on 4 cores |
| 10 | `sims/26_full_kappa1.R`, `sims/27_coverage_kappa1.R`, `sims/37_full_kappa05.R`, `sims/38_coverage_kappa05.R` | tolerance of 1 and 0.5 standard errors | 8 h |
| 11 | `sims/20_sens_shift.R` | `output/sens_shift/` (reduced overlap) | 1 h |
| 12 | `sims/24_sens_n1000.R`, `sims/28_sens_n1000_full.R` | `output/sens_n1000/` (1000 patients per trial) | several hours |
| 13 | `sims/14_main_summary.R` | `output/summary/` | 5 min |
| 14 | `sims/35_best_rmse_by_scenario.R`, `sims/36_tableS1_bias_rmse.R` | lowest RMSE by scenario; Supplementary Table S1 | 2 min |
| 15 | `sims/33_fig1_schematic.R`, `sims/31_figures_main.R` | Figures 1–4 | 1 min |

Steps 6, 8 and 9 recompute methods on the same replicates as the main run, with the same seeds, after corrections made during the study. `sims/14` replaces the corresponding rows of the main run with them. `sims/14` and `sims/35` read the per-replicate files. `sims/36`, `sims/31` and `sims/33` read only `output/summary/`, so Supplementary Table S1 and the figures can be rebuilt from the files in this repository.

## Method codes

| Code (reporting format) | Method in the paper |
|---|---|
| M0 (R1) | Unadjusted (Bucher) |
| M1m (R1) | MAIC (means) |
| M1 (R1) | MAIC (means and SDs) |
| M1 (R5) | MAIC (all reported summaries) |
| M1q (R5) | MAIC (means, quartiles and proportions) |
| M7 (R5) | sMAIC; M7k1 and M7k05 use tolerances of 1 and 0.5 standard errors |
| M2 (R1) | STC |
| M3a (R1) | G-computation |
| M5 (R1) | G-MAIC (means and SDs) |
| M5 (R5) | G-MAIC |
| M5s (R5) | sG-MAIC; M5sk1 and M5sk05 use tolerances of 1 and 0.5 standard errors |

Reporting formats: R1 mean (SD); R2 median [Q1, Q3]; R2r median (range); R3 mean (SD) and n (%) above 65; R3m mean (SD) and n (%) above 60; R4 mean (SD) and median [Q1, Q3]; R5 all statistics.

## Notes

Code comments are written in Chinese. Function and variable names are in English. Some comments refer to sections of the study protocol or to numbered decisions in the authors' working notes; these documents are not part of the repository.

## Licence

MIT. See `LICENSE`.
