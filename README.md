# Is severe care dependency in later life a metropolitan rather than a rural burden in Korea? Evidence from long-term care insurance certifications

Analysis code for a district-level study of severe care dependency in later life in South Korea: residents aged 65 and older certified at long-term care insurance grade 1 or 2 across the 229 districts (si-gun-gu), 2011-2024, analysed with Bayesian spatiotemporal models (R-INLA, BYM2).

This code is released under the MIT licence. Version 1.0.1. Archived copy: https://doi.org/10.5281/zenodo.24556722 (all versions: https://doi.org/10.5281/zenodo.24552897). Version 1.0.1 adds the master-panel script `45_master_from_refdataset.py` and fixes folder creation and multiverse deferral so that the code runs from a fresh copy.

## Contents

This repository contains code and model-settings files only. No data, adjacency graph, boundary file, result tables or CSV files are included. All source statistics are public (KOSIS, Statistics Korea); the processed district-year data are available from the author on reasonable request.

| Path | What it is |
|---|---|
| `00_선별/01_스크립트/` | KOSIS download script, screening-stage counts and the judgement check |
| `01_스크립트/` | Outcome and expected counts, analytic panel, model scripts, sensitivity suite and `_var_labels.json` (variable labels) |
| `04_분석/LTCS/` | Covariate pool builders, settings writer, model settings `sets_LTCS_*.json` (principal `std3`), run scripts and the robustness/numbers summary |
| `04_분석/LTCS/gate/code/` | Village sports facility correction and the cell-level data check |

Folder names and some column and region names are in Korean because the source tables are Korean; every script starts with an English header comment (purpose, inputs, outputs).

## Data sources (KOSIS table IDs)

| Table | Content | Used for |
|---|---|---|
| DT_35006_N009 (org 350) | Long-term care certified persons by district, sex, age group and grade | Outcome (grade 1-2, aged 65+) and judgement check |
| DT_35006_N036 (org 350) | National certifications by sex and age band | National rates for the expected counts |
| DT_35006_N005, DT_35006_N037 (org 350) | Long-term care applications (district; national by sex and age) | Screening-stage comparison outcome only |
| DT_11761_N006, DT_11761_N009, DT_11761_N002, DT_11761_N007 (org 117) | Registered persons with severe disability | Screening-stage comparison outcome only |
| DT_1YL20651E | 2023 population of the integrated Changwon city | Population override for one district-year |

Covariates are read from a reference district-year dataset (SQLite, read-only) compiled from public district statistics (KOSIS local indicators, Community Health Survey, cause-of-death statistics, health-resource and public assistance statistics, population and housing census); see the variable table below.

## Model

The principal model (`S1_narrative_full`, version `std3`) contains 10 hypothesised district characteristics (care supply and access, application and recognition pathway, deprivation, survival, walking) and 5 fixed adjustments, with a BYM2 spatial term, a first-order random walk in time and an unstructured year term; offset = log expected count (no added constant); the likelihood is chosen by WAIC among Poisson and negative binomial forced-adjustment-only models. A covariate is called robust when it is credible in the principal model, in at least 75% of the 200 multiverse subsets with a consistent sign in at least 95%, and credible with the same sign in at least 80% of the alternative analyses (`90_ltcs_summary_numbers.py`).

| Variable | Role | Description | Source |
|---|---|---|---|
| `ltcs_cases` | outcome | Residents aged 65+ certified at long-term care grade 1 or 2 | KOSIS DT_35006_N009 (NHIS long-term care insurance statistics) |
| `ltcs_pop65` | offset base | Resident population aged 65+ | Resident registration statistics |
| `E_std (E_sex, E_std_age, E_crude)` | offset | Expected count, indirect standardisation by 2 sexes x 5 age bands (alternatives: sex-only, age-only, crude) | Derived: KOSIS DT_35006_N036 national certifications by sex and age; resident population |
| `ltc_hospitals_per100k_65` | hypothesis | Long-term care hospitals per 100 000 residents aged 65+ | Health-resource statistics |
| `hospitals_per100k` | hypothesis | Hospitals and general hospitals per 100 000 residents | Health-resource statistics |
| `tertiary_educated_share` | hypothesis | University-educated share: four-year university or higher (current students and non-completers included, junior college excluded) among residents aged 6+ | Population and housing census (KOSIS DT_1PM1501, DT_1PM2001), interpolated between censuses |
| `living_alone_elderly` | hypothesis | Older single-person household share | Population and housing census |
| `mort_stroke_asr` | hypothesis | Cerebrovascular disease mortality, age-standardised | Cause-of-death statistics |
| `alz_death_sir_ma3` | hypothesis | Alzheimer disease death SIR, 3-year moving average | Cause-of-death statistics |
| `basic_livelihood_rate` | hypothesis | Public assistance (basic livelihood security) recipient rate | Public assistance statistics |
| `medical_aid1_per100k` | hypothesis | Medical aid type 1 patients per 100 000 residents | Medical aid statistics |
| `mort_allcause_asr` | hypothesis | All-cause mortality, age-standardised | Cause-of-death statistics |
| `walking_rate` | hypothesis | Walking practice rate | Community Health Survey |
| `urban_population_share, physicians_per_1000, sex_ratio, fiscal_independence, pop80_share` | fixed adjustment | Urban population share, physicians, sex ratio, fiscal independence, share aged 80+ | KOSIS local indicators; resident registration statistics |
| `hit_and_run_rate, traffic_culture_index, fires_per_10000` | negative control | Hit-and-run accident rate, traffic culture index, fires | KOSIS local indicators |

Sensitivity versions (`sets_LTCS_<version>.json`): `s577` queen contiguity graph, `s2` without the share aged 80+, `sc` crude expected counts, `s8` age-only standardised expected counts, `y1524` 2015-2024, `x2021` without 2020-2021, `y1324` 2013-2024, `sexE` sex-only expected counts plus share aged 85+, `supply` care institutions and hospital beds added.

## Environment variables

| Variable | Meaning |
|---|---|
| `ST_ROOT` | Working root holding `01_스크립트/`, `02_패널/` and `04_분석/`; the run scripts default to the repository root |
| `ST_KOSIS_RAW` | Folder of downloaded KOSIS JSON files (`<tblId>/<year>.json`) |
| `KOSIS_API_KEY` | The user's own KOSIS Open API key (download script only; no key is included) |
| `ST_FBDZIP` | Folder of resident population by sex and 5-year age group, one file per year (`<year>.csv`) |
| `ST_DATA_ROOT` | Folder holding the district-year population panel `식중독최종.csv` (229 districts) read by the district key scripts |
| `ST_REF_DATASET` | Path of the reference district-year dataset (SQLite file, opened read-only) |
| `ST_SHP` | District boundary shapefile (229 districts, 2024 boundaries); default `data/boundary/final.shp` |
| `ST_GRAPH` | Principal adjacency graph in INLA graph format (229 districts, 593 edges, boundary-file order) |
| `ST_VER`, `TAG` | Model version (`std3` or a sensitivity version) and outcome tag (`LTCS`), set by the run scripts |
| `SK_RSLOTS` | Number of concurrent R jobs in `rslot.sh` (default 3) |
| `ST_THREADS` | INLA threads per model (the run scripts use 1) |
| `ST_MV_DEFER`, `ST_MV_PARTS` | `ST_MV_DEFER=1` postpones the 200-subset multiverse in `81_fit_inla_ltcs.R`; `ST_MV_PARTS=k` with `ST_MV_RANGE=a-b` runs it in k blocks that the main run then merges |
| `ST_MV_N`, `ST_SENS_LITE` | Multiverse size and a light mode used for the sensitivity versions (set in `run_ltcs_chain.sh`) |
| `EXG` | Group of additional analyses run by `88X_ltcs_extra.R` (`A`, `B`, `C`) |

To unpack the archive on macOS, use `ditto -x -k <zip> <folder>` or Python `zipfile`; the system `unzip` may fail on the Korean folder names.

## Reproduction order

Run from `ST_ROOT` after setting the environment variables above.

| Step | Script | Output |
|---|---|---|
| 0 | `00_선별/01_스크립트/kosis_fetch.py` | Download of the KOSIS source tables (own API key in `KOSIS_API_KEY`) into `ST_KOSIS_RAW` |
| 1 | `00_선별/01_스크립트/99S_build_E.py`, `99S_judgement_gate.py` | Screening-stage counts for the candidate outcomes and the judgement check (grade composition by year, new grades in 2014 and 2018, year-to-year district rank correlation) |
| 1b | `01_스크립트/45_master_from_refdataset.py` | District-year master panel `02_패널/MASTER_시군구_연도_wide.parquet` read from table `source_master_all_years` of `ST_REF_DATASET` (the same table that `50_build_master_panel.py` builds from source-specific wide files, which are not distributed) |
| 2 | `01_스크립트/99L_build_Estd_ltcs.py` (imports `40_region_crosswalk.py` and `50_build_master_panel.py`, which read `ST_DATA_ROOT` at import, so set it even though this step does not open the file) | Outcome counts (grade 1-2 certifications, residents aged 65+), offset population and expected counts by indirect standardisation (2 sexes x 5 age bands); alternative expected counts |
| 3 | `04_분석/LTCS/gate/code/06_village_corrections.py`, `04_분석/LTCS/build_pool_from_sqlite_ltcs.py`, `build_pool_ltcs_additions.py`, `gate/code/11_ltcs_cell_gate.py` | District-year covariate pool from the reference dataset and the cell-level data check |
| 4 | `04_분석/LTCS/make_pool_sets_ltcs.py` | Model settings `sets_LTCS_*.json` and hypothesis table, hashed (SHA-256) before the covariate models were fitted |
| 5 | `ST_VER=std3 python3 01_스크립트/80_build_disease_panel_v.py LTCS` | Analytic panel (imputation counts, VIF) |
| 6 | `04_분석/LTCS/run_ltcs_chain.sh`, then `run_ltcs_chain_part2.sh` | Principal model, 200-subset multiverse, six likelihoods, period models, emerging hot spots, Knorr-Held models, Spatial+, sensitivity suite, nine sensitivity versions, extra sensitivity analyses, robustness verdicts and `numbers_ltcs.json` |

## Output map

Result tables are written to `04_분석/LTCS/results_<version>/` unless another folder is shown.

| Script | Output |
|---|---|
| `99L_build_Estd_ltcs.py` | 02_패널/LTCS_Estd3.csv |
| `build_pool_from_sqlite_ltcs.py / build_pool_ltcs_additions.py` | 02_패널/EXT_canonical_pool_LTCS_base.csv, 02_패널/EXT_canonical_pool_LTCS.csv |
| `make_pool_sets_ltcs.py` | 04_분석/LTCS/sets_LTCS_<version>.json, hypothesis_table_LTCS.csv, sets_LTCS_hash_before_fit.txt |
| `80_build_disease_panel_v.py` | 04_분석/LTCS/panel_LTCS_<version>.csv, panel_META_*, panel_INFO_*, config_* |
| `81_fit_inla_ltcs.R` | results_<version>/S1_likelihood.csv, S2 set comparison (HM_setcompare_*), S3_main_IRR.csv, S3_main_fit.csv, S3_negcontrol.csv, S6_spatial_RR.csv, S8_residual_moran.csv, S4_window_*.csv, S5_covid.csv, S6_moran_year.csv, S7_multiverse_*.csv, S9_typeI_IRR.csv |
| `81R_ltcs_refit_lik.R` | S0_repro_meta.csv, S1b_refit_twice.csv, S1c_likelihood5.csv, S1c_likelihood5_coef.csv |
| `85 / 86 / 86b / 87` | S4b_window_*.csv; S11_ehsa*.csv, S11b_moran_yearly_EB.csv, S12*_persistence*.csv; S27_ehsa_mk_*.csv; S4c_window_grid_table.csv |
| `88g_sens_prep.py / 88g_sens_suite.R / 89g_sens_combine.py` | S41_*.csv |
| `88h_ppc_varpart.R` | S42_ppc.csv, S42_pit_hist.csv, S42_variance_partition.csv |
| `88k_struct_kh.R` | S43_part_<analysis>.csv, S50_*.csv |
| `88X_ltcs_extra.R` | S60_extra_*.csv |
| `90_ltcs_summary_numbers.py` | numbers_ltcs.json, results_std3/S90_robustness_verdict.csv, result_signature_LTCS.json |

R jobs are wrapped in `rslot.sh` and run with one INLA thread (`ST_THREADS=1`) so that WAIC comparisons are reproducible.

## Software

R 4.6.0 with R-INLA 25.10.19, sf 1.1.1, spdep 1.4.2, jsonlite, ggplot2 and dplyr; Python 3.9 with pandas 2.3, numpy 2.0, scipy, statsmodels, geopandas 1.0, pymannkendall, matplotlib and pyarrow. Newer INLA versions may shift DIC/WAIC by a few points without changing the rate ratios.

## Licence

MIT (see `LICENSE`).
