#!/usr/bin/env python3
# 99S_judgement_gate.py - Judgement check for the candidate outcomes: national grade composition of long-term care certifications by year (DT_35006_N009 total rows),
#  breaks of +/-30% and years when new grades were introduced (grade 5 in 2014, cognitive-support grade in 2018); year-to-year rank correlation of district SIRs.
#  Input: ST_KOSIS_RAW, 02_패널/<TAG>_Estd3n.csv. Output: 99S_judgement_*.csv
import os, json, numpy as np, pandas as pd
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SEL = os.path.dirname(ROOT); RAW = os.path.join(SEL, "kosis_raw")
rows = []
for y in range(2010, 2026):
    D = pd.DataFrame(json.load(open(os.path.join(RAW, "DT_35006_N009", f"{y}.json"))))
    D = D[(D.C1 == "15350AC400Z1") & (D.C2_NM == "합계")]; v = pd.to_numeric(D.DT, errors="coerce").groupby(D.C4_NM).sum()
    out = {"year": y, **{k: v.get(k, np.nan) for k in ["계", "1등급", "2등급", "3등급", "4등급", "5등급", "인지지원등급", "등급외", "등급외A", "등급외B", "등급외C"]}}
    rows.append(out)
G = pd.DataFrame(rows).set_index("year")
G["인정(1-5+인지)"] = G[["1등급", "2등급", "3등급", "4등급", "5등급", "인지지원등급"]].sum(1, min_count=1)
G["1-2등급"] = G["1등급"] + G["2등급"]
for c in ["1등급", "2등급", "3등급", "4등급", "5등급", "인지지원등급"]: G[f"{c}_비율(인정중)"] = (G[c] / G["인정(1-5+인지)"]).round(4)
for c in ["1등급", "2등급", "1-2등급", "인정(1-5+인지)"]: G[f"{c}_전년비"] = (G[c] / G[c].shift()).round(3)
G["단절표시"] = ""
for c in ["1등급", "2등급", "1-2등급", "인정(1-5+인지)"]:
    G.loc[(G[f"{c}_전년비"] - 1).abs() >= 0.30, "단절표시"] += f"{c} ±30%; "
G.loc[2014, "단절표시"] += "4·5등급 신설(시행령 2014-07-01, 3등급 하한 51→60); "; G.loc[2018, "단절표시"] += "인지지원등급 신설(2018-01-01); "
G.loc[2012, "단절표시"] += "3등급 하한 55→53(2012-07-01); "; G.loc[2013, "단절표시"] += "3등급 하한 53→51(2013-07-01); "
G.to_csv(os.path.join(SEL, "99S_judgement_grade_composition.csv"), encoding="utf-8-sig")
print(G[["1등급", "2등급", "1-2등급", "인정(1-5+인지)", "1등급_전년비", "2등급_전년비", "1-2등급_전년비", "인정(1-5+인지)_전년비", "단절표시"]].to_string())
P = []
for T in ["LTCS", "LTCC", "LTCA", "DISV"]:
    X = pd.read_csv(os.path.join(ROOT, "02_패널", f"{T}_Estd3n.csv"), encoding="utf-8-sig"); X["s"] = np.log(X.cases / X.E_std)
    W = X.pivot(index="fbd_region", columns="year", values="s")
    for y in W.columns[:-1]:
        r = spearmanr(W[y], W[y + 1]).correlation
        q = (X[X.year == y + 1].set_index("fbd_region").cases / X[X.year == y].set_index("fbd_region").cases)
        P.append(dict(TAG=T, 연도쌍=f"{y}-{y+1}", rho_순위=round(r, 3), 시군구비_P10=round(q.quantile(.1), 3), 시군구비_중앙=round(q.median(), 3), 시군구비_P90=round(q.quantile(.9), 3)))
P = pd.DataFrame(P); P.to_csv(os.path.join(SEL, "99S_judgement_yearpair_rank.csv"), index=False, encoding="utf-8-sig")
print(P.pivot(index="연도쌍", columns="TAG", values="rho_순위").to_string())
print(P[(P.TAG == "DISV") & (P.연도쌍 == "2018-2019")].to_string())
