#!/usr/bin/env python3
# 99S_prep.py - Screening stage: district-year panels for the screened candidate outcomes (cases, offset, expected counts) with the fixed adjustments from the
#  reference dataset (ST_REF_DATASET, read-only; public assistance rate 2021-2023 collapse rule: < 0.2x or > 5x the 2020/2024 mean set to missing and interpolated).
#  The LTCA and LTCC panels are the comparison (negative-control) outcomes of the sensitivity suite. Output: 04_분석/<TAG>/panel_<TAG>.csv, 02_패널/<TAG>_Estd.csv
import os, json, shutil, sqlite3, unicodedata
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); AN = os.path.join(ROOT, "04_분석"); PAN = os.path.join(ROOT, "02_패널")
DB = os.environ["ST_REF_DATASET"]
NFC = lambda s: unicodedata.normalize("NFC", str(s))
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
base = pd.read_sql("select * from disease_panel where tag='PERIO'", c).drop(columns=["tag"], errors="ignore"); base["fbd_region"] = base.fbd_region.map(NFC)
base = base.sort_values(["fbd_region", "year"]).reset_index(drop=True)
_ref = base[base.year.isin([2020, 2024])].groupby("fbd_region").basic_livelihood_rate_raw.mean(); _r = base.fbd_region.map(_ref)
_bad = base.year.isin([2021, 2022, 2023]) & ((base.basic_livelihood_rate_raw < 0.2 * _r) | (base.basic_livelihood_rate_raw > 5 * _r))
print("[99S-prep] 수급률 위반 칸", int(_bad.sum()), base.loc[_bad, ["fbd_region", "year", "basic_livelihood_rate_raw"]].values.tolist())
base.loc[_bad, "basic_livelihood_rate_raw"] = np.nan
base["basic_livelihood_rate_raw"] = base.groupby("fbd_region").basic_livelihood_rate_raw.transform(lambda s: s.interpolate(method="linear", limit_area="inside"))
_s = base.basic_livelihood_rate_raw; base["basic_livelihood_rate"] = (_s - _s.mean()) / _s.std(ddof=0)
_r2 = base.fbd_region.map(base[base.year.isin([2020, 2024])].groupby("fbd_region").basic_livelihood_rate_raw.mean()); _yy = base.year.isin([2021, 2022, 2023])
assert not (_yy & ((base.basic_livelihood_rate_raw < 0.2 * _r2) | (base.basic_livelihood_rate_raw > 5 * _r2))).any() and base.basic_livelihood_rate_raw.notna().all()
print("[99S-prep] 창원 정정 후", base.loc[base.fbd_region.eq(NFC("경상남도창원시")) & base.year.between(2020, 2024), ["year", "basic_livelihood_rate_raw"]].round(3).values.tolist())
cfg = json.loads(c.execute("select content from model_settings where tag='PERIO' and kind='config'").fetchone()[0])
cfg["forced"] = ["care_institutions_per_100k", "urban_population_share", "older_adult_share", "sex_ratio", "fiscal_independence"]
for T in ["LTCS", "LTCC", "LTCA", "DISV"]:
    X = pd.read_csv(os.path.join(PAN, f"{T}_Estd3n.csv"), encoding="utf-8-sig"); X["fbd_region"] = X.fbd_region.map(NFC); X = X[X.cases.notna()]
    D = os.path.join(AN, T); os.makedirs(os.path.join(D, "results"), exist_ok=True)
    P = base.drop(columns=["cases", "E", "E_crude", "SIR", "offset"], errors="ignore").merge(X[["year", "fbd_region", "offset", "cases", "E_crude"]], on=["year", "fbd_region"], how="inner", validate="one_to_one").rename(columns={"E_crude": "E"})
    P["cases"] = P.cases.astype(int); P["SIR"] = P.cases / P.E
    assert not P.duplicated(["year", "fbd_region"]).any() and P[["cases", "E", "offset"]].notna().all().all() and (P.cases >= 0).all() and (P.E > 0).all()
    assert P.groupby("year").fbd_region.nunique().min() >= 228
    cov = [c for c in cfg["forced"] if c in P.columns]; assert P[cov].notna().all().all(), T
    P.to_csv(os.path.join(D, f"panel_{T}.csv"), index=False, encoding="utf-8-sig")
    X[["year", "fbd_region", "E_std"]].to_csv(os.path.join(PAN, f"{T}_Estd.csv"), index=False, encoding="utf-8-sig")
    json.dump(cfg, open(os.path.join(D, f"config_{T}.json"), "w"), ensure_ascii=False, indent=1)
    shutil.copy2(os.path.join(AN, "METS_std", "results", "nb_METS_std.graph"), os.path.join(D, "results", f"nb_{T}.graph"))
    print(f"[99S-prep] {T}: 행 {len(P)}, {P.year.min()}-{P.year.max()}, 시군구 {P.fbd_region.nunique()}, cases {int(P.cases.sum()):,}")
