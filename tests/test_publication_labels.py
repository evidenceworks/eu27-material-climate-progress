#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
REQ={
    "CMUR","DMC","RMC","TERR","CFOOT","D1","D2","D3","D_TOTAL",
    "FULL_2010_2023","EARLY_2010_2016","LATE_2017_2023","CMUR_even_year",
}

def main():
    d=pd.read_csv(ROOT/"docs"/"PUBLICATION_LABELS.csv",keep_default_na=False)
    assert d.internal_key.is_unique
    assert REQ.issubset(set(d.internal_key))
    bykey=d.set_index("internal_key")
    assert bykey.loc["CMUR_even_year","public_label"]=="Common even-year sensitivity"
    assert bykey.loc["CMUR_even_year","short_public_label"]=="Common even-year"
    assert bykey.loc["FULL_2010_2023","public_label"]=="Full (2010–2023)"
    assert bykey.loc["EARLY_2010_2016","public_label"]=="Early (2010–2016)"
    assert bykey.loc["LATE_2017_2023","public_label"]=="Late (2017–2023)"
    print(f"PUBLICATION_LABELS_TEST_PASS rows={len(d)}")

if __name__=="__main__":
    main()
