#!/usr/bin/env python3
# make_pool_sets_ltcs.py - Writes the model-settings files sets_LTCS_<version>.json (principal std3; sensitivity versions s577, s2, sc, s8, y1524, x2021, y1324, sexE, supply),
#  the hypothesis table and their SHA-256 hashes (sets_LTCS_hash_before_fit.txt), fixed before the covariate models were fitted.
#  Principal set S1_narrative_full = 10 hypothesised district characteristics; 5 fixed adjustments; 3 negative controls; sensitivity-only variables.
import os, json, hashlib
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); PAN = os.path.join(ROOT, "02_패널")
POOL = "EXT_canonical_pool_LTCS.csv"
P = pd.read_csv(os.path.join(PAN, POOL), encoding="utf-8-sig"); assert len(P) == 229 * 14 and P.drop(columns=["vil_n_sqlite", "vil_n_corrected"], errors="ignore").notna().all().all()
FORCED = ["urban_population_share", "physicians_per_1000", "sex_ratio", "fiscal_independence", "pop80_share"]
NEG = ["hit_and_run_rate", "traffic_culture_index", "fires_per_10000"]
SENS = ["public_health_use_rate", "care_institutions_per100k", "beds_per_1000", "outside_care_share", "pop85_share", "obesity_rate", "smoking_rate", "high_risk_drinking", "screening_uptake", "alz_death_sir"]
HYP = ["ltc_hospitals_per100k_65", "hospitals_per100k", "tertiary_educated_share", "living_alone_elderly", "mort_stroke_asr", "alz_death_sir_ma3",
       "basic_livelihood_rate", "medical_aid1_per100k", "mort_allcause_asr", "walking_rate"]
X = lambda L: ["ext__" + v for v in L]; XK = lambda v: "ext__" + v
EXPECT = {
 "ltc_hospitals_per100k_65": "- (long-term care hospitals per 100 000 aged 65+: substitution outside long-term care insurance; alternative + = supply-induced referral; two-sided interpretation fixed)",
 "hospitals_per100k": "+ (hospitals and general hospitals per 100 000: diagnosis and referral access; alternative - = acute inpatient substitution)",
 "tertiary_educated_share": "+ (tertiary-educated share: information and application capacity; alternative - = better health; census interpolated, 2021-24 carried forward)",
 "living_alone_elderly": "+ (older single-person households: absence of family care leads to formal applications; 2011-14 carried back)",
 "mort_stroke_asr": "+ (cerebrovascular age-standardised mortality: stroke burden, severe dependency)",
 "alz_death_sir_ma3": "+ (Alzheimer disease death SIR, 3-year moving average: dementia burden)",
 "basic_livelihood_rate": "+ (public assistance recipient rate: deprivation, multimorbidity; alternative - = pro-rich gradient of formal care use; two-sided)",
 "medical_aid1_per100k": "+ (medical aid type 1 patients per 100 000: deprivation; alternative - = application barriers; two-sided)",
 "mort_allcause_asr": "- (all-cause age-standardised mortality: survival is needed to reach severe dependency)",
 "walking_rate": "- (walking practice: functional maintenance, one factor among several)",
}
EVKEY = {
 "ltc_hospitals_per100k_65": "kim2015delivery; kim2020variations; sakong2017panel; seo2021effect",
 "hospitals_per100k": "jin2022regional; sakong2017panel; seo2021effect",
 "tertiary_educated_share": "park2021determinants; park2023cognitive; kreft2016expansion",
 "living_alone_elderly": "lagergren2014comparing; sakong2017panel",
 "mort_stroke_asr": "kang2011factors",
 "alz_death_sir_ma3": "kang2011factors; park2023cognitive",
 "basic_livelihood_rate": "davies2025arealevel; wang2025incidence; lera2020socioeconomic",
 "medical_aid1_per100k": "noh2026socioeconomic; davies2025arealevel; lera2020socioeconomic",
 "mort_allcause_asr": "noh2026socioeconomic; yoshikawa2023regional",
 "walking_rate": "takahashi2022risk",
}
SETS = {"S1_narrative_full": X(HYP),
        "S2_core": X(["ltc_hospitals_per100k_65", "tertiary_educated_share", "living_alone_elderly", "mort_stroke_asr", "basic_livelihood_rate", "mort_allcause_asr"]),
        "S3_supply_access": X(["ltc_hospitals_per100k_65", "hospitals_per100k"]),
        "S4_application_path": X(["tertiary_educated_share", "living_alone_elderly", "mort_stroke_asr", "alz_death_sir_ma3"]),
        "S5_deprivation": X(["basic_livelihood_rate", "medical_aid1_per100k"]),
        "S6_survival_health": X(["mort_allcause_asr", "mort_stroke_asr", "alz_death_sir_ma3", "walking_rate"])}
pov = [{"year": 2023, "fbd_region": "경상남도창원시", "population": 1009038, "source": "KOSIS DT_1YL20651E, 2023 total for the integrated Changwon city"}]
allv = ["ltcs_cases", "ltcs_pop65"] + FORCED + NEG + SENS + HYP
LAB_M = "sex-and-age-standardised expected grade 1-2 long-term care certifications among residents aged 65+, 2 sexes x 5 age bands (65-69, 70-74, 75-79, 80-84, 85+) (national rates DT_35006_N036 2013-2024, 2013 rates for 2011-2012; district resident population by sex and age), 2011-2024"
WIN = {"W1_2011-15": [2011, 2015], "W2_2016-19": [2016, 2019], "W3_2020-24": [2020, 2024]}
def make():
    return {
     "tag": "LTCS", "outcome_column": "ext__ltcs_cases",
     "outcome_label_en": "Residents aged 65+ certified at long-term care grade 1 or 2 under Korea's long-term care insurance (NHIS long-term care statistics DT_35006_N009, annual)",
     "disease_en": "Severe care dependency (long-term care grade 1-2 certification) in later life",
     "disease_ko": "65세 이상 장기요양 1–2등급 인정(중증 돌봄 의존)",
     "pathway": "care-dependency epidemiology set: certification = need x application and assessment pathway. Care supply and access (long-term care hospitals, hospitals); application and recognition pathway (tertiary education, older people living alone, stroke and dementia burden); deprivation (public assistance, medical aid type 1); survival (all-cause mortality); walking as one factor among several",
     "pool_note": "Number of long-term care institutions is not in the reference dataset; long-term care benefit users are excluded (circular with the outcome); outside-district care share is sensitivity-only; public health centre use is sensitivity-only (2023-2024 break)",
     "selection_basis": "Outcome selected by spatial screening criteria fixed before results (adjusted Moran's I >= 0.25, permutation p < 0.001, urban-rural R2 < 0.40, BYM2 mixing >= 0.5, persistence >= 0.5)",
     "offset": "ext__ltcs_pop65", "exclude_years": [],
     "forced": X(FORCED), "negative_controls": X(NEG),
     "offset_E_column": "E_std", "offset_E_file": "LTCS_Estd3.csv", "offset_E_label": LAB_M,
     "graph": "Principal adjacency graph in INLA graph format, read from the ST_GRAPH environment variable (229 districts, 593 edges, boundary-file order)",
     "population_override": pov, "sensitivity_only_extra": X(SENS),
     "version": "std3 care-dependency epidemiology set (care supply and access, application and recognition pathway, deprivation, survival; walking as one factor among several), 2011-2024",
     "fixed_on": "hashed (SHA-256) and recorded in the analysis log before the covariate models were fitted",
     "sets": SETS, "transform": {XK(v): "none" for v in FORCED + NEG + SENS + HYP},
     "expected_direction": {XK(k): v for k, v in EXPECT.items()},
     "evidence": {XK(k): "literature keys: " + EVKEY[k] for k in HYP},
     "rules": {"vif_action": "drop", "likelihood": "Poisson or negative binomial, whichever forced-adjustment-only model has the lower WAIC", "set_choice": "Lowest WAIC; within 4 WAIC points, the set with fewer variables", "credible": "95% CrI excludes 1", "vif_drop": 10, "multiverse_specs": 200,
               "main_model_note": "Principal model = S1_narrative_full (10 hypothesised characteristics + 5 fixed adjustments); the set is not changed after fitting; candidates with VIF > 10 are removed by 80 (vif_action drop); the chain stops if a fixed adjustment keeps VIF > 10.",
               "years": [2011, 2024], "windows": WIN, "covid_post_from": 2020, "windows_note": "Three periods fixed before results were seen."},
     "note": "LTCS = residents aged 65+ certified at long-term care insurance grade 1 (95 points or more) or grade 2 (75 to under 95 points); the grade 1-2 thresholds are unchanged 2008-2025. Grade 5 (2014) and the cognitive-support grade (2018) were added below 75 points (national grade 1-2 ratio 2013 to 2014: 1.009).",
     "extra_columns": {"ext__" + v: [POOL, v] for v in allv}}
OUT = {}
VERS = {
 "std3": None,
 "s577": ("sensitivity: queen contiguity graph (577 edges)", {}),
 "s2": ("sensitivity: no fixed adjustment for the share aged 80 or older", {"forced": X(FORCED[:4])}),
 "sc": ("sensitivity: crude expected counts (residents aged 65+ x national annual rate)", {"offset_E_column": "E_crude", "offset_E_label": "crude E (residents aged 65+ x national annual rate)"}),
 "s8": ("sensitivity: age-only standardised expected counts (5 bands, 65+)", {"offset_E_column": "E_std_age", "offset_E_label": "age-only standardised E"}),
 "y1524": ("sensitivity: 2015-2024 (after grade 5 was introduced in 2014)", {"exclude_years": [2011, 2012, 2013, 2014]}),
 "x2021": ("sensitivity: 2020-2021 excluded (COVID-19 period, deferred reassessments)", {"exclude_years": [2020, 2021]}),
 "y1324": ("sensitivity: 2013-2024 (years with national age-specific rates)", {"exclude_years": [2011, 2012]}),
 "sexE": ("sensitivity: sex-only standardised expected counts plus share aged 85+ as a fixed adjustment", {"offset_E_column": "E_sex", "offset_E_label": "sex-only standardised E (residents aged 65+ by sex x national sex-specific rate)", "forced": X(FORCED + ["pop85_share"])}),
 "supply": ("sensitivity: care institutions per 100 000 and hospital beds per 1000 added as fixed adjustments", {"forced": X(FORCED + ["care_institutions_per100k", "beds_per_1000"])}),
}
for ver, spec in VERS.items():
    j = make()
    if spec:
        j["sens_note"] = spec[0]; j.update(spec[1])
        if "forced" in spec[1]: j["sensitivity_only_extra"] = [v for v in j["sensitivity_only_extra"] if v not in j["forced"]]
    OUT[ver] = j
H = []
for ver, j in OUT.items():
    fp = os.path.join(HERE, f"sets_LTCS_{ver}.json"); json.dump(j, open(fp, "w"), ensure_ascii=False, indent=1)
    H.append(f"sets_LTCS_{ver}.json sha256={hashlib.sha256(open(fp, 'rb').read()).hexdigest()}")
open(os.path.join(HERE, "sets_LTCS_hash_before_fit.txt"), "w").write("# sets_LTCS json and hypothesis table SHA-256, before fitting\n" + "\n".join(H) + "\n")
rows = [dict(variable=XK(v), role="hypothesis", expected_direction=EXPECT[v], evidence_keys=EVKEY[v], in_S1=True) for v in HYP] + [dict(variable=XK(v), role="forced_control", expected_direction="fixed adjustment (no directional hypothesis)", in_S1=True) for v in FORCED] + \
       [dict(variable=XK(v), role="negative_control", expected_direction="0 (no association expected)", in_S1=False) for v in NEG] + [dict(variable=XK(v), role="sensitivity_only", expected_direction="sensitivity only", in_S1=False) for v in SENS]
pd.DataFrame(rows).to_csv(os.path.join(HERE, "hypothesis_table_LTCS.csv"), index=False, encoding="utf-8-sig")
H.append("hypothesis_table_LTCS.csv sha256=" + hashlib.sha256(open(os.path.join(HERE, "hypothesis_table_LTCS.csv"), "rb").read()).hexdigest())
open(os.path.join(HERE, "sets_LTCS_hash_before_fit.txt"), "w").write("# sets_LTCS json and hypothesis table SHA-256, before fitting\n" + "\n".join(H) + "\n")
print("\n".join(H))
