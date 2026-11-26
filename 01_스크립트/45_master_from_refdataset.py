# 45_master_from_refdataset.py
#  Purpose: write the district-year master panel 02_패널/MASTER_시군구_연도_wide.parquet (and .csv.gz) used by
#           80_build_disease_panel_v.py, by reading table `source_master_all_years` of the reference district-year
#           dataset (SQLite, opened read-only). This replaces the upstream wide-CSV build (50_build_master_panel.py),
#           whose source-specific inputs are not distributed; both give the same table.
#  Input:  ST_REF_DATASET (SQLite path), table source_master_all_years
#  Output: $ST_ROOT/02_패널/MASTER_시군구_연도_wide.parquet, MASTER_시군구_연도_wide.csv.gz
import os, sqlite3, unicodedata
import pandas as pd
ROOT = os.environ.get("ST_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DB = os.environ["ST_REF_DATASET"]
PAN = os.path.join(ROOT, "02_패널"); os.makedirs(PAN, exist_ok=True)
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
M = pd.read_sql_query("SELECT * FROM source_master_all_years", con); con.close()
M.columns = [unicodedata.normalize("NFC", c) for c in M.columns]
assert {"year", "key", "fbd_region", "population"} <= set(M.columns), "master skeleton columns missing"
assert M.groupby("year").key.nunique().nunique() == 1, "district count differs by year"
M.to_parquet(os.path.join(PAN, "MASTER_시군구_연도_wide.parquet"), index=False)
M.to_csv(os.path.join(PAN, "MASTER_시군구_연도_wide.csv.gz"), index=False, encoding="utf-8-sig", compression="gzip")
print(f"[45] MASTER {M.shape[0]} rows x {M.shape[1]} columns, years {M.year.min()}-{M.year.max()} -> {PAN}")
