#!/usr/bin/env python3
# 11_ltcs_cell_gate.py - Cell-level data check before fitting: year-to-year jumps (neighbour ratio < 0.2 or > 5; 3-fold steps), public assistance rate 2021-2023 rule,
#  boundary-change districts, denominator breaks, interpolation traces. Output: 04_분석/LTCS/gate/tables/G02_*, G05_*, G06_*, G_summary.json
import os, json, sqlite3, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); AN = os.path.dirname(B); ROOT = os.path.dirname(os.path.dirname(AN)); T = os.path.join(B, "tables")
P = pd.read_csv(os.path.join(ROOT, "02_패널", "EXT_canonical_pool_LTCS.csv"), encoding="utf-8-sig", dtype={"sigungu_cd": str})
P["rate"] = P.ltcs_cases / P.ltcs_pop65 * 1e3
VARS = [c for c in P.columns if c not in ("year", "fbd_region", "sigungu_cd")]
yrs = list(range(2011, 2025)); rows = []
for v in VARS:
    W = P.pivot(index="fbd_region", columns="year", values=v)
    for r in W.index:
        s = W.loc[r]
        for i, y in enumerate(yrs):
            x = s[y]
            if pd.isna(x): continue
            prv = [t for t in yrs if t < y and not pd.isna(s[t])]; nxt = [t for t in yrs if t > y and not pd.isna(s[t])]
            p = s[prv[-1]] if prv else np.nan; q = s[nxt[0]] if nxt else np.nan
            fl = []
            if prv and nxt and p > 0 and q > 0:
                m = (p + q) / 2
                if x < 0.2 * m or x > 5 * m: fl.append("a_이웃평균0.2/5배")
                if (x >= 3 * p and x >= 3 * q) or (x <= p / 3 and x <= q / 3): fl.append("b_앞뒤둘다3배")
            if prv and p > 0 and (x >= 3 * p or x <= p / 3) and nxt and (q >= 3 * p or q <= p / 3): fl.append("c_3배계단유지")
            if (not prv or not nxt) and (p if prv else q) > 0:
                o = p if prv else q
                if x >= 3 * o or x <= o / 3: fl.append("d_끝단3배")
            if x == 0: fl.append("zero")
            if fl: rows.append(dict(var=v, fbd_region=r, year=y, value=x, prev=p, next=q, flags=";".join(fl)))
J = pd.DataFrame(rows); J.to_csv(os.path.join(T, "G02_칸급변_플래그.csv"), index=False, encoding="utf-8-sig")
print("[G02] 플래그 칸:", len(J)); print(J.groupby("var").flags.apply(lambda s: s.str.split(";").explode().value_counts().to_dict()).to_string()[:3000])
BL = P.pivot(index="fbd_region", columns="year", values="basic_livelihood_rate"); ref = (BL[2020] + BL[2024]) / 2
viol = [(r, y, BL.at[r, y], ref[r]) for r in BL.index for y in (2021, 2022, 2023) if BL.at[r, y] < 0.2 * ref[r] or BL.at[r, y] > 5 * ref[r]]
print("[G04] 수급률 2021–23 규칙 위반 칸:", len(viol), viol[:5])
BD = P[P.fbd_region.str.contains("군위|세종|청주|부천")].sort_values(["fbd_region", "year"])
BD.to_csv(os.path.join(T, "G05_경계변경4곳_전변수_시계열.csv"), index=False, encoding="utf-8-sig")
gw = P[P.fbd_region.str.contains("군위")].set_index("year")
G5 = pd.DataFrame({v: [gw.at[2022, v], gw.at[2023, v], gw.at[2024, v], gw.at[2023, v] / ((gw.at[2022, v] + gw.at[2024, v]) / 2) if (gw.at[2022, v] + gw.at[2024, v]) else np.nan] for v in VARS}, index=["2022", "2023", "2024", "2023/이웃평균"]).T
G5.to_csv(os.path.join(T, "G05_군위2023_전변수.csv"), encoding="utf-8-sig"); print("[G05] 군위 2023/이웃평균 0.67–1.5 밖:\n", G5[(G5["2023/이웃평균"] < 0.67) | (G5["2023/이웃평균"] > 1.5)].round(3).to_string())
print("[G05] 세종 행", int(P.fbd_region.str.contains("세종").sum()), "개(14 기대)")
U = []
for v in VARS:
    g = P.groupby("year")[v]; med = g.median(); wm = P.groupby("year").apply(lambda d: np.average(d[v], weights=d.population))
    jump = (wm / wm.shift()).dropna(); z0 = P[P[v] == 0].groupby("year").size()
    U.append(dict(var=v, min=P[v].min(), median=P[v].median(), max=P[v].max(), wmean_2011=wm.iloc[0], wmean_2024=wm.iloc[-1], max_yoy_ratio=jump.max(), min_yoy_ratio=jump.min(),
                  yoy_break_years=";".join(str(y) for y, r_ in jump.items() if r_ > 1.3 or r_ < 0.77), zero_cells=int((P[v] == 0).sum()), zero_years=";".join(f"{y}:{n}" for y, n in z0.items()), over100=int((P[v] > 100).sum())))
U = pd.DataFrame(U); U.to_csv(os.path.join(T, "G06_단위분모_연도단절.csv"), index=False, encoding="utf-8-sig"); print("[G06]\n", U.round(3).to_string())
OB = P.pivot(index="fbd_region", columns="year", values="obesity_rate"); eq = ((OB[2019] - (OB[2018] + OB[2020]) / 2).abs() < 1e-9).sum()
print("[G07] 2019 비만율 = 2018·2020 평균인 칸:", int(eq), "/ 229")
for v in ["moderate_pa_rate", "public_sport_fac_per100k", "walking_rate", "screening_uptake"]:
    W = P.pivot(index="fbd_region", columns="year", values=v); e = {y: int(((W[y] - (W[y - 1] + W[y + 1]) / 2).abs() < 1e-9).sum()) for y in range(2012, 2024)}
    print(f"[G07] {v} 앞뒤 평균과 같은 칸(보간 흔적) 연도별:", {k: n for k, n in e.items() if n > 20})
json.dump(dict(n_flag=int(len(J)), bl_violations=len(viol), obes2019_interp=int(eq)), open(os.path.join(T, "G_summary.json"), "w"))
