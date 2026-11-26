#!/usr/bin/env python3
# build_pool_ltcs_additions.py - Adds the care-dependency covariates to the base pool from the reference dataset (read-only): hospitals, public health centre use,
#  medical aid type 1 patients, all-cause and stroke age-standardised mortality, 3-year moving average of the Alzheimer death SIR;
#  hospital cells missing for cities with general districts filled from the earliest observed year (logged).
#  Output: 02_패널/EXT_canonical_pool_LTCS.csv, gate/G00b_base_hash.txt, gate/tables/G00b_hospitals_filled_cells.csv
import os, hashlib, sqlite3, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); PAN = os.path.join(ROOT, "02_패널"); GT = os.path.join(HERE, "gate")
DB = os.environ["ST_REF_DATASET"]
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest(); h0 = sha(DB)
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
P = pd.read_csv(os.path.join(PAN, "EXT_canonical_pool_LTCS_base.csv"), encoding="utf-8-sig", dtype={"sigungu_cd": str}); P["fbd_region"] = P.fbd_region.map(NFC)
M = {"general_hosp_per100k": "general_hospitals_per_100k", "public_health_use_rate": "public_health_centre_use_rate", "medical_aid1_per100k": "medical_aid1_patients_per_100k",
     "beds_per_1000": "beds_per_1000", "outside_care_share": "outside_care_share", "ltc_hospitals_per100k_all": "long_term_care_hospitals_per_100k"}
D = pd.read_sql("select tag, year, fbd_region, " + ",".join(f'"{v}_raw"' for v in M.values()) + " from disease_panel where year between 2011 and 2024", con); D["fbd_region"] = D.fbd_region.map(NFC)
log = []
for k, v in M.items():
    c = v + "_raw"; ok = D.groupby("tag")[c].apply(lambda s: s.notna().sum() == 3206); tags = list(ok[ok].index); assert tags, k
    W = D[D.tag.isin(tags)].pivot_table(index=["year", "fbd_region"], columns="tag", values=c); sp = float((W.max(1) - W.min(1)).abs().max()); assert sp < 1e-9, (k, sp)
    P = P.merge(W.iloc[:, 0].rename(k).reset_index(), on=["year", "fbd_region"], how="left", validate="one_to_one"); log.append((k, c, len(tags), sp))
RY = pd.read_sql('select year, fbd_region, "mort__전체사망_asr" as mort_allcause_asr, "mort__뇌혈관질환_asr" as mort_stroke_asr from region_year where year between 2011 and 2024', con); RY["fbd_region"] = RY.fbd_region.map(NFC)
P = P.merge(RY, on=["year", "fbd_region"], how="left", validate="one_to_one")
A = pd.read_sql("select year, fbd_region, cases, E from disease_panel where tag='ALZ' and year between 2011 and 2024", con); A["fbd_region"] = A.fbd_region.map(NFC); A["alz_death_sir"] = A.cases / A.E
assert (A.E > 0).all(); P = P.merge(A[["year", "fbd_region", "alz_death_sir"]], on=["year", "fbd_region"], how="left", validate="one_to_one")
H = pd.read_sql('select year, fbd_region, "chk__요양기관_병원" as h, "chk__요양기관_종합병원" as g from region_year where year between 2011 and 2024', sqlite3.connect(f"file:{DB}?mode=ro", uri=True)); H["fbd_region"] = H.fbd_region.map(NFC)
P = P.merge(H, on=["year", "fbd_region"], how="left", validate="one_to_one"); P["hospitals_per100k"] = (P.h + P.g) / P.population * 1e5
fm = P.hospitals_per100k.isna(); FILL = P.loc[fm, ["year", "fbd_region"]].copy()
P["hospitals_per100k"] = P.groupby("fbd_region").hospitals_per100k.transform(lambda s: s.bfill()); FILL["value_filled"] = P.loc[fm, "hospitals_per100k"].values; FILL["rule"] = "원표 시군구 행 없음(일반구만) → 2013 값 이월"
FILL.assign(var="hospitals_per100k").to_csv(os.path.join(GT, "tables", "G00b_hospitals_filled_cells.csv"), index=False, encoding="utf-8-sig"); print("병원 채움 칸", len(FILL)); P = P.drop(columns=["h", "g"])
P = P.sort_values(["fbd_region", "year"]); P["alz_death_sir_ma3"] = P.groupby("fbd_region").alz_death_sir.transform(lambda s: s.rolling(3, center=True, min_periods=2).mean()); P = P.sort_values(["year", "fbd_region"]).reset_index(drop=True)
miss = P.drop(columns=["vil_n_sqlite", "vil_n_corrected"], errors="ignore").isna().sum(); assert miss.sum() == 0, miss[miss > 0]
P.to_csv(os.path.join(PAN, "EXT_canonical_pool_LTCS.csv"), index=False, encoding="utf-8-sig")
con.close(); h1 = sha(DB)
open(os.path.join(GT, "G00b_base_hash.txt"), "w").write(f"기준 SQLite 전 {h0}\n후 {h1}\n같음={h0 == h1}\n" + "\n".join(f"{a}<-{b} 태그 {c}개 최대차 {d}" for a, b, c, d in log) + "\n")
print("pool", P.shape, "해시 같음", h0 == h1); print(P[list(M) + ["mort_allcause_asr", "mort_stroke_asr", "alz_death_sir"]].describe().T[["min", "50%", "max"]].round(2))
