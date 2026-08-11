#!/usr/bin/env python3
"""Build publication-facing source CSVs and artwork from validated analysis outputs.

This module is intentionally a publication layer around the validated scientific
pipeline. It does not estimate trends, ranks, correlations, contrasts, or
robustness statistics. It reads those validated objects and reorganizes them
into publication-facing tables/figures.
"""
from __future__ import annotations

from pathlib import Path
import argparse, ast, hashlib, math, tempfile, zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.offsetbox import AnnotationBbox, DrawingArea
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
import shapefile
from shapely.geometry import shape as shapely_shape, box, Point
from shapely.ops import transform as shp_transform
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
EU27 = ["AT","BE","BG","HR","CY","CZ","DK","EE","FI","FR","DE","EL","HU","IE","IT","LV","LT","LU","MT","NL","PL","PT","RO","SK","SI","ES","SE"]
PERIOD_ORDER = ["FULL_2010_2023","EARLY_2010_2016","LATE_2017_2023"]
PERIOD_LABEL = {"FULL_2010_2023":"2010–2023","EARLY_2010_2016":"2010–2016","LATE_2017_2023":"2017–2023"}
OKABE_BLUE = "#0072B2"
OKABE_ORANGE = "#E69F00"
DARK = "#333333"
GRID = "#E7E7E7"
W178 = 178/25.4
W85 = 85/25.4
DIVERGE = LinearSegmentedColormap.from_list("blue_white_orange", [OKABE_BLUE, "#F5F5F5", OKABE_ORANGE])
RANK_CMAP = plt.get_cmap("cividis_r")

matplotlib.rcParams.update({
    "font.family":"DejaVu Sans",
    "font.size":7.0,
    "axes.titlesize":7.5,
    "axes.labelsize":7.0,
    "xtick.labelsize":6.2,
    "ytick.labelsize":6.2,
    "legend.fontsize":6.0,
    "pdf.fonttype":42,
    "ps.fonttype":42,
    "svg.fonttype":"none",
    "svg.hashsalt":"eu27-material-climate-progress-publication",
})

FIGURE_FILES = {
    "figure_1_primary_concordance_grid":"figure_1_primary_concordance_grid",
    "figure_2_rank_concordance_geometry":"figure_2_rank_concordance_geometry",
    "figure_3_country_representation_consequences":"figure_3_country_representation_consequences",
    "figure_4_geography_progress_ranks":"figure_4_geography_progress_ranks",
    "figure_5_temporal_representation_sensitivity":"figure_5_temporal_representation_sensitivity",
    "figure_6_deletion_robustness":"figure_6_deletion_robustness",
    "figure_s1_temporal_concordance_matrices":"figure_s1_temporal_concordance_matrices",
    "figure_s2_geography_representation_shifts":"figure_s2_geography_representation_shifts",
    "figure_s3_geography_progress_scores":"figure_s3_geography_progress_scores",
    "figure_s4_progress_score_geometry":"figure_s4_progress_score_geometry",
    "figure_s5_broader_robustness_summary":"figure_s5_broader_robustness_summary",
    "figure_s6_sign_flip_overlap":"figure_s6_sign_flip_overlap",
    "figure_s7_country_discordance_matrix":"figure_s7_country_discordance_matrix",
}


def _read(analysis: Path, rel: str) -> pd.DataFrame:
    p = analysis / rel
    if not p.is_file():
        raise FileNotFoundError(p)
    return pd.read_csv(p)


def _primary_country(analysis: Path) -> pd.DataFrame:
    d = _read(analysis, "supplementary/country_results_all_specs.csv")
    q = d[(d["period"] == "FULL_2010_2023") & (d["analysis_label"] == "primary")].copy()
    if len(q) != 27 or sorted(q.geo.tolist()) != sorted(EU27):
        raise RuntimeError("Primary country-result set is not exactly EU27")
    return q.sort_values("geo").reset_index(drop=True)


def _primary_grid_row(analysis: Path, correlation: str = "spearman") -> pd.Series:
    d = _read(analysis, "tables/representation_grid.csv")
    q = d[(d.period == "FULL_2010_2023") & (d.spec == "primary_common_pc") & (d.estimator == "OLS") & (d.cmur_mode == "pp") & (d.correlation == correlation)]
    if len(q) != 1:
        raise RuntimeError(f"Expected one primary {correlation} grid row")
    return q.iloc[0]


def _cell_rows_from_wide(df: pd.DataFrame) -> pd.DataFrame:
    mapping = [
        ("a_CMUR_TERR","CMUR","Territorial GHG"),
        ("b_DMC_TERR","DMC","Territorial GHG"),
        ("c_RMC_TERR","RMC","Territorial GHG"),
        ("d_CMUR_CFOOT","CMUR","Consumption-based GHG footprint"),
        ("e_DMC_CFOOT","DMC","Consumption-based GHG footprint"),
        ("f_RMC_CFOOT","RMC","Consumption-based GHG footprint"),
    ]
    rows=[]
    for _,r in df.iterrows():
        for col,mat,clim in mapping:
            rows.append({"period":r["period"],"material_representation":mat,"climate_perspective":clim,"spearman_rho":float(r[col])})
    return pd.DataFrame(rows)


def build_publication_sources(analysis: Path, output_root: Path) -> None:
    """Create publication source CSVs from the analysis outputs."""
    figsd = output_root / "figures" / "source_data"
    tables = output_root / "results" / "reference" / "tables"
    supp = output_root / "supplementary_tables"
    figsd.mkdir(parents=True, exist_ok=True); tables.mkdir(parents=True, exist_ok=True); supp.mkdir(parents=True, exist_ok=True)

    main = _read(analysis, "tables/main_results.csv")
    country = _primary_country(analysis)
    rg = _read(analysis, "tables/representation_grid.csv")
    rank_summary = _read(analysis, "tables/country_rank_change_summary.csv")
    robustness = _read(analysis, "supplementary/robustness_summary.csv")
    loco = _read(analysis, "supplementary/leave_one_country_out.csv")
    loyo = _read(analysis, "supplementary/leave_one_year_out.csv")
    cmur_even = _read(analysis, "supplementary/cmur_even_year.csv")

    # MAIN Table 2: source facts, with units/provenance aligned to the validated data dictionary.
    t2 = pd.DataFrame([
        ["CMUR","env_ac_cur","Circular material use rate","Percent","OLS slope in percentage points/year","Higher is favorable","Economy-wide circular use share","Official Eurostat series; annual 2010–2023"],
        ["DMC","env_ac_mfa","Domestic material consumption","Tonnes per capita","Negative OLS slope of log per-capita level","Lower is favorable","Domestic material consumption","Common average-population denominator used in primary analysis"],
        ["RMC","env_ac_rme","Raw material consumption (material footprint)","Tonnes per capita","Negative OLS slope of log per-capita level","Lower is favorable","Consumption-side material footprint","Official Eurostat environmental-account series; status flags summarized in Table S7"],
        ["Territorial GHG","sdg_13_10","Territorial greenhouse-gas emissions","Tonnes CO2e per capita","Negative OLS slope of log per-capita level","Lower is favorable","Territorial emissions excluding LULUCF and memo items","Common average-population denominator used in primary analysis"],
        ["Consumption-based GHG footprint","cli_gge_foot","Consumption-based greenhouse-gas footprint","Tonnes CO2e per capita","Negative OLS slope of log per-capita level","Lower is favorable","Consumption-accounting footprint","Eurostat footprint series; common average-population denominator used in primary analysis"],
    ], columns=["representation","official_dataset","construct","unit","primary_progress_metric","favorable_direction","accounting_boundary","status_provenance_note"])
    t2.to_csv(tables/"table_2_representation_contract_source.csv", index=False)

    # MAIN Table 3: six FULL cells plus ordered contrasts.
    r = main[main.period == "FULL_2010_2023"].iloc[0]
    cells = [
        ("a_CMUR_TERR","CMUR × Territorial GHG"),("b_DMC_TERR","DMC × Territorial GHG"),("c_RMC_TERR","RMC × Territorial GHG"),
        ("d_CMUR_CFOOT","CMUR × Consumption-based GHG footprint"),("e_DMC_CFOOT","DMC × Consumption-based GHG footprint"),("f_RMC_CFOOT","RMC × Consumption-based GHG footprint")]
    rows=[]
    for key,label in cells:
        rows.append(["concordance",key,label,float(r[key]),"FULL_2010_2023","Spearman",27,"Ranks use deterministic 12-decimal canonicalization before tie handling"])
    for key,label in [("D1_b_minus_a","D1: CMUR → DMC"),("D2_c_minus_b","D2: DMC → RMC"),("D3_f_minus_c","D3: Territorial → consumption-based GHG"),("D_TOTAL_f_minus_a","Overall representation contrast")]:
        rows.append(["ordered_contrast",key,label,float(r[key]),"FULL_2010_2023","Difference in Spearman rho",27,"Ordered contrast of validated concordance cells"])
    pd.DataFrame(rows,columns=["object_type","internal_key","public_label","value","period","metric","n_countries","tie_metadata"]).to_csv(tables/"table_3_primary_concordance_and_contrasts.csv",index=False)

    # MAIN Table 4: validated summary only; no new country-consequence statistic is added.
    t4 = rank_summary.rename(columns={
        "substitution":"representation_substitution","max_country":"country_at_maximum","sign_flips_n":"progress_sign_disagreements_n"})
    labelmap={"CMUR_to_DMC":"CMUR → DMC","DMC_to_RMC":"DMC → RMC","TERR_to_CFOOT":"Territorial → consumption-based GHG"}
    t4.insert(1,"public_label",t4.representation_substitution.map(labelmap))
    t4.to_csv(tables/"table_4_country_consequences.csv",index=False)

    # MAIN Table 5: long-format temporal and robustness evidence.
    t5=[]
    for _,rr in main.iterrows():
        for key,label in [("D1_b_minus_a","D1: CMUR → DMC"),("D2_c_minus_b","D2: DMC → RMC"),("D3_f_minus_c","D3: Territorial → consumption-based GHG"),("D_TOTAL_f_minus_a","Overall representation contrast")]:
            t5.append(["A_temporal",rr.period,label,float(rr[key]),float(rr[key]),0])
    # Primary and selected robustness rows for D_TOTAL.
    full=float(main[main.period=="FULL_2010_2023"].iloc[0].D_TOTAL_f_minus_a)
    t5.append(["B_full_robustness","primary_FULL","Overall representation contrast",full,full,0])
    for _,rr in robustness[robustness.estimand=="D_TOTAL_f_minus_a"].iterrows():
        t5.append(["B_full_robustness",rr.diagnostic,"Overall representation contrast",float(rr["min"]),float(rr["max"]),int(rr.sign_reversals)])
    # Common even-year is a validated standalone sensitivity.
    if len(cmur_even)==1:
        v=float(cmur_even.iloc[0].D_TOTAL_f_minus_a)
        t5.append(["B_full_robustness","CMUR_even_year","Overall representation contrast",v,v,0])
    pd.DataFrame(t5,columns=["panel","period_or_diagnostic","estimand","min_value","max_value","sign_reversals"]).drop_duplicates().to_csv(tables/"table_5_temporal_and_robustness.csv",index=False)

    # Supplement Table S1: source series and human-readable unit/role support.
    sources=pd.read_csv(ROOT/"data"/"sources.csv")
    dd=pd.read_csv(ROOT/"data"/"data_dictionary.csv")
    s1=sources[["source_id","dataset","role","raw_file","api_query","resolved_selectors","rows","countries","years"]].copy()
    def _selector_label(raw, dimension):
        try:
            items=ast.literal_eval(raw) if isinstance(raw,str) else []
            for item in items:
                if item.get("dimension")==dimension:
                    return item.get("label") or item.get("code") or ""
        except Exception:
            return ""
        return ""
    s1["unit"]=s1.resolved_selectors.map(lambda x:_selector_label(x,"unit"))
    s1.loc[s1.source_id.eq("average_population_demo_gind") & s1.unit.eq(""),"unit"]="Persons (average population)"
    s1["indicator_or_measure"]=s1["role"]
    s1["provenance_note"]="Frozen Eurostat API response included in data/raw; source query and resolved selectors recorded in data/sources.csv"
    s1=s1[["source_id","dataset","indicator_or_measure","unit","role","raw_file","api_query","rows","countries","years","provenance_note"]]
    s1.to_csv(supp/"table_s1_sources_units_provenance.csv",index=False)

    # S2 detailed temporal matrices.
    s2=_cell_rows_from_wide(main)
    s2["period_label"]=s2.period.map(PERIOD_LABEL)
    s2.to_csv(supp/"table_s2_temporal_concordance_matrices.csv",index=False)

    # S3 Kendall/Pearson companions, with related Spearman reference.
    comp=[]
    base=rg[(rg.spec=="primary_common_pc")&(rg.estimator=="OLS")&(rg.cmur_mode=="pp")]
    for period in PERIOD_ORDER:
        sr=base[(base.period==period)&(base.correlation=="spearman")].iloc[0]
        for corr in ["kendall","pearson"]:
            rr=base[(base.period==period)&(base.correlation==corr)].iloc[0]
            for key,mat,clim in [("a_CMUR_TERR","CMUR","Territorial GHG"),("b_DMC_TERR","DMC","Territorial GHG"),("c_RMC_TERR","RMC","Territorial GHG"),("d_CMUR_CFOOT","CMUR","Consumption-based GHG footprint"),("e_DMC_CFOOT","DMC","Consumption-based GHG footprint"),("f_RMC_CFOOT","RMC","Consumption-based GHG footprint")]:
                comp.append([period,mat,clim,corr,float(rr[key]),float(sr[key]),27])
    pd.DataFrame(comp,columns=["period","material_representation","climate_perspective","correlation","correlation_value","spearman_reference","n_countries"]).to_csv(supp/"table_s3_companion_correlations.csv",index=False)

    # S4 detailed ordered contrasts.
    s4=[]
    for _,rr in main.iterrows():
        for key,label in [("D1_b_minus_a","D1: CMUR → DMC"),("D2_c_minus_b","D2: DMC → RMC"),("D3_f_minus_c","D3: Territorial → consumption-based GHG"),("D_TOTAL_f_minus_a","Overall representation contrast")]:
            s4.append([rr.period,key,label,float(rr[key]),float(rr.identity_D_sum_error) if "identity_D_sum_error" in rr.index else 0.0])
    pd.DataFrame(s4,columns=["period","internal_key","public_label","value","identity_check_error"]).to_csv(supp/"table_s4_temporal_ordered_contrasts.csv",index=False)

    # S5 extended robustness: validated full-period specifications + LOCO/LOYO ranges + common even-year sensitivity.
    s5=[]
    full_rg=rg[(rg.period=="FULL_2010_2023")&(rg.correlation=="spearman")]
    for _,rr in full_rg.iterrows():
        key=f"{rr.spec}|{rr.estimator}|{rr.cmur_mode}"
        for e in ["D1_b_minus_a","D2_c_minus_b","D3_f_minus_c","D_TOTAL_f_minus_a"]:
            s5.append(["fixed_specification",key,e,float(rr[e]),float(rr[e]),0])
    for _,rr in robustness.iterrows():
        s5.append(["deletion_range",rr.diagnostic,rr.estimand,float(rr["min"]),float(rr["max"]),int(rr.sign_reversals)])
    ce=cmur_even.iloc[0]
    for e in ["D1_b_minus_a","D2_c_minus_b","D3_f_minus_c","D_TOTAL_f_minus_a"]:
        s5.append(["fixed_specification","CMUR_even_year",e,float(ce[e]),float(ce[e]),0])
    pd.DataFrame(s5,columns=["record_type","specification","estimand","min_value","max_value","sign_reversals"]).to_csv(supp/"table_s5_extended_robustness.csv",index=False)

    _read(analysis,"supplementary/accounting_bridge.csv").to_csv(supp/"table_s6_production_account_bridge.csv",index=False)
    _read(analysis,"supplementary/rmc_status.csv").to_csv(supp/"table_s7_rmc_status_provenance.csv",index=False)

    # S8 concise diagnostics (summaries, not the full Data S1 detail dump).
    s8=[]
    pcs=_read(analysis,"supplementary/published_vs_common_rank_summary.csv")
    for _,rr in pcs.iterrows():
        s8.append(["published_vs_common",rr.period,rr.representation,"spearman_progress_common_vs_published",float(rr.spearman_progress_common_vs_published),""])
        s8.append(["published_vs_common",rr.period,rr.representation,"mean_abs_rank_displacement",float(rr.mean_abs_rank_displacement),""])
        s8.append(["published_vs_common",rr.period,rr.representation,"max_abs_rank_displacement",float(rr.max_abs_rank_displacement),""])
    ties=_read(analysis,"supplementary/tie_audit.csv")
    for _,rr in ties.iterrows():
        note=f"{rr['variable']}; n_tied={int(rr['n_tied'])}; geos={rr['geos']}; raw_span={rr['raw_span']:.17g}"
        s8.append(["tie_audit",rr.period,rr.spec,"canonical_value_12dp",float(rr.canonical_value_12dp),note])
    s8.append(["CMUR_even_year","FULL_2010_2023","CMUR_even_year","D_TOTAL_f_minus_a",float(ce.D_TOTAL_f_minus_a),str(ce.years)])
    total=full_rg[(full_rg.spec=="total_scale")&(full_rg.estimator=="OLS")&(full_rg.cmur_mode=="pp")]
    if len(total)==1:
        s8.append(["total_scale","FULL_2010_2023","total_scale","D_TOTAL_f_minus_a",float(total.iloc[0].D_TOTAL_f_minus_a),"Validated total-scale sensitivity"])
    pd.DataFrame(s8,columns=["record_type","period","representation_or_spec","metric","value","notes"]).to_csv(supp/"table_s8_numeric_scale_diagnostics.csv",index=False)

    # Figure source data.
    f1=[]
    for key,mat,clim in [("a_CMUR_TERR","CMUR","Territorial GHG"),("b_DMC_TERR","DMC","Territorial GHG"),("c_RMC_TERR","RMC","Territorial GHG"),("d_CMUR_CFOOT","CMUR","Consumption-based GHG footprint"),("e_DMC_CFOOT","DMC","Consumption-based GHG footprint"),("f_RMC_CFOOT","RMC","Consumption-based GHG footprint")]:
        f1.append([mat,clim,float(r[key])])
    pd.DataFrame(f1,columns=["material_representation","climate_perspective","spearman_rho"]).to_csv(figsd/"figure_1_primary_concordance_grid.csv",index=False)

    # Figure 2: all-country rank geometry, six panels.
    rr=[]
    for mat in ["CMUR","DMC","RMC"]:
        for clim,clabel in [("TERR","Territorial GHG"),("CFOOT","Consumption-based GHG footprint")]:
            rho=float(r[["a_CMUR_TERR","b_DMC_TERR","c_RMC_TERR","d_CMUR_CFOOT","e_DMC_CFOOT","f_RMC_CFOOT"][["CMUR_TERR","DMC_TERR","RMC_TERR","CMUR_CFOOT","DMC_CFOOT","RMC_CFOOT"].index(f"{mat}_{clim}")]])
            for _,cc in country.iterrows():
                rr.append([f"{mat}_{clim}",f"{mat} × {clabel}",cc.geo,float(cc[f"rank_{mat}"]),float(cc[f"rank_{clim}"]),rho])
    pd.DataFrame(rr,columns=["panel","panel_label","geo","material_rank","climate_rank","spearman_rho"]).to_csv(figsd/"figure_2_rank_concordance_geometry.csv",index=False)

    cc_cols=["geo","rank_shift_CMUR_to_DMC","rank_shift_DMC_to_RMC","rank_shift_TERR_to_CFOOT","sign_flip_CMUR_DMC","sign_flip_DMC_RMC","sign_flip_TERR_CFOOT"]
    country[cc_cols].to_csv(figsd/"figure_3_country_representation_consequences.csv",index=False)
    country[["geo","rank_CMUR","rank_DMC","rank_RMC","rank_TERR","rank_CFOOT"]].to_csv(figsd/"figure_4_geography_progress_ranks.csv",index=False)
    main.to_csv(figsd/"figure_5_temporal_representation_sensitivity.csv",index=False)

    # Figure 6 deletion distribution: FULL + every FULL LOCO/LOYO realization, all four contrasts.
    f6=[]
    for e in ["D1_b_minus_a","D2_c_minus_b","D3_f_minus_c","D_TOTAL_f_minus_a"]:
        f6.append([e,"FULL","FULL",float(r[e])])
        for _,x in loco[loco.period=="FULL_2010_2023"].iterrows(): f6.append([e,"LOCO",str(x.omit_geo),float(x[e])])
        for _,x in loyo[loyo.period=="FULL_2010_2023"].iterrows(): f6.append([e,"LOYO",str(int(x.omit_year)),float(x[e])])
    pd.DataFrame(f6,columns=["estimand","source","omitted_unit","value"]).to_csv(figsd/"figure_6_deletion_robustness.csv",index=False)

    main.to_csv(figsd/"figure_s1_temporal_concordance_matrices.csv",index=False)
    country[cc_cols].to_csv(figsd/"figure_s2_geography_representation_shifts.csv",index=False)
    country[["geo","CMUR","DMC","RMC","TERR","CFOOT"]].to_csv(figsd/"figure_s3_geography_progress_scores.csv",index=False)

    # S4 progress-score geometry with matching rank-correlation cells.
    ps=[]
    prow=_primary_grid_row(analysis,"pearson"); srow=_primary_grid_row(analysis,"spearman")
    colkeys={"CMUR_TERR":"a_CMUR_TERR","DMC_TERR":"b_DMC_TERR","RMC_TERR":"c_RMC_TERR","CMUR_CFOOT":"d_CMUR_CFOOT","DMC_CFOOT":"e_DMC_CFOOT","RMC_CFOOT":"f_RMC_CFOOT"}
    for panel,key in colkeys.items():
        mat,clim=panel.split("_")
        for _,cc in country.iterrows():
            ps.append([panel,cc.geo,float(cc[mat]),float(cc[clim]),float(srow[key]),float(prow[key])])
    pd.DataFrame(ps,columns=["panel","geo","material_progress_score","climate_progress_score","spearman_rho","pearson_r"]).to_csv(figsd/"figure_s4_progress_score_geometry.csv",index=False)

    # S5 D_TOTAL robustness summary.
    f5=[["primary_FULL","D_TOTAL_f_minus_a",full,full,0]]
    for _,rrr in robustness[robustness.estimand=="D_TOTAL_f_minus_a"].iterrows(): f5.append([rrr.diagnostic,rrr.estimand,float(rrr["min"]),float(rrr["max"]),int(rrr.sign_reversals)])
    # Add selected fixed specs absent from robustness_summary.
    for label,spec,est,mode in [("published_pc","published_pc","OLS","pp"),("common_log_CMUR","primary_common_pc","OLS","log"),("Theil_Sen","primary_common_pc","THEIL_SEN","pp"),("total_scale","total_scale","OLS","pp")]:
        q=full_rg[(full_rg.spec==spec)&(full_rg.estimator==est)&(full_rg.cmur_mode==mode)]
        if len(q)==1:
            v=float(q.iloc[0].D_TOTAL_f_minus_a); f5.append([label,"D_TOTAL_f_minus_a",v,v,0])
    f5.append(["CMUR_even_year","D_TOTAL_f_minus_a",float(ce.D_TOTAL_f_minus_a),float(ce.D_TOTAL_f_minus_a),0])
    pd.DataFrame(f5,columns=["diagnostic","estimand","min","max","sign_reversals"]).drop_duplicates().to_csv(figsd/"figure_s5_broader_robustness_summary.csv",index=False)

    # S6 sign-flip intersections.
    sf=country.groupby(["sign_flip_CMUR_DMC","sign_flip_DMC_RMC","sign_flip_TERR_CFOOT"],dropna=False).agg(n_countries=("geo","size"),countries=("geo",lambda x:", ".join(sorted(x)))).reset_index().sort_values(["n_countries","sign_flip_CMUR_DMC","sign_flip_DMC_RMC","sign_flip_TERR_CFOOT"],ascending=[False,True,True,True])
    sf.to_csv(figsd/"figure_s6_sign_flip_overlap.csv",index=False)

    # S7 country-level discordance matrix: material rank minus climate rank.
    cd=[]
    for _,cc in country.iterrows():
        for mat in ["CMUR","DMC","RMC"]:
            for clim in ["TERR","CFOOT"]:
                cd.append([cc.geo,f"{mat}×{clim}",float(cc[f"rank_{mat}"]-cc[f"rank_{clim}"])])
    pd.DataFrame(cd,columns=["geo","cell","material_minus_climate_rank"]).to_csv(figsd/"figure_s7_country_discordance_matrix.csv",index=False)


def _natural_earth_rows(ne_zip: Path):
    with tempfile.TemporaryDirectory(prefix="natural_earth_") as td:
        with zipfile.ZipFile(ne_zip) as z:
            if z.testzip() is not None: raise RuntimeError("Natural Earth ZIP CRC failure")
            z.extractall(td)
        shp=next(Path(td).glob("*.shp"))
        rdr=shapefile.Reader(str(shp),encoding="latin1")
        fields=[f[0] for f in rdr.fields[1:]]
        for sr in rdr.iterShapeRecords():
            rec=dict(zip(fields,list(sr.record)))
            yield rec, shapely_shape(sr.shape.__geo_interface__)


def build_country_code_map(ne_zip: Path, output_root: Path) -> tuple[pd.DataFrame, dict]:
    records=list(_natural_earth_rows(ne_zip))
    mapping=[]; geoms={}
    clip=box(-12.8,34.0,35.8,72.5)
    transformer=Transformer.from_crs("EPSG:4326","EPSG:3035",always_xy=True)
    for sci in EU27:
        if sci=="FR":
            cand=[(r,g) for r,g in records if str(r.get("ADM0_A3"))=="FRA" and str(r.get("ADMIN"))=="France"]
            method="ADM0_A3=FRA + ADMIN=France (avoids French dependency records)"
            ne_a2="FR"
        elif sci=="EL":
            cand=[(r,g) for r,g in records if str(r.get("ISO_A2"))=="GR" and str(r.get("ADMIN"))=="Greece"]
            method="Eurostat EL explicitly mapped to Natural Earth ISO_A2=GR"
            ne_a2="GR"
        else:
            cand=[(r,g) for r,g in records if str(r.get("ISO_A2"))==sci]
            method="Natural Earth ISO_A2 exact match"
            ne_a2=sci
        if len(cand)!=1:
            raise RuntimeError(f"Natural Earth join for {sci} returned {len(cand)} rows")
        r,g=cand[0]
        clipped=g.intersection(clip)
        if clipped.is_empty: raise RuntimeError(f"Natural Earth geometry empty after Europe crop: {sci}")
        proj=shp_transform(transformer.transform,clipped).simplify(2500,preserve_topology=True)
        geoms[sci]=proj
        mapping.append([sci,ne_a2,r.get("ADM0_A3"),r.get("ADMIN"),method,1,0,"PASS"])
    df=pd.DataFrame(mapping,columns=["scientific_geo","natural_earth_iso_a2","natural_earth_adm0_a3","natural_earth_admin","match_method","matched_geometry_n","duplicate_matches_n","status"])
    if len(df)!=27 or df.scientific_geo.nunique()!=27 or (df.status!="PASS").any(): raise RuntimeError("Natural Earth 27/27 join check failed")
    (output_root/"cartography").mkdir(parents=True,exist_ok=True)
    df.to_csv(output_root/"cartography"/"COUNTRY_CODE_MAP.csv",index=False)
    parity=df.copy()
    parity["expected_eu27_n"]=27; parity["matched_eu27_n"]=27; parity["missing_eu27_n"]=0; parity["duplicate_country_matches_n"]=0; parity["scientific_value_recomputed"]=False
    parity.to_csv(output_root/"MAP_JOIN_PARITY.csv",index=False)
    return df, geoms


def _poly_patches(geom, **kwargs):
    geoms=list(geom.geoms) if geom.geom_type=="MultiPolygon" else [geom]
    return [MplPolygon(np.asarray(p.exterior.coords),closed=True,**kwargs) for p in geoms if not p.is_empty]


def _map_extent(geoms):
    bounds=np.array([g.bounds for g in geoms.values()])
    xmin,ymin,xmax,ymax=bounds[:,0].min(),bounds[:,1].min(),bounds[:,2].max(),bounds[:,3].max()
    dx=xmax-xmin; dy=ymax-ymin
    return xmin-.025*dx,xmax+.025*dx,ymin-.02*dy,ymax+.02*dy


def _draw_map(ax, geoms, values, cmap, norm, extent):
    for c in EU27:
        col=cmap(norm(values[c])) if c in values and pd.notna(values[c]) else "#F4F4F4"
        for patch in _poly_patches(geoms[c],facecolor=col,edgecolor="#777777",linewidth=.32): ax.add_patch(patch)
    ax.set_xlim(extent[0],extent[1]); ax.set_ylim(extent[2],extent[3]); ax.set_aspect("equal"); ax.axis("off")


def _save(fig, output_root: Path, stem: str):
    d=output_root/"figures"/"reference"; d.mkdir(parents=True,exist_ok=True)
    meta={"Creator":"eu27-material-climate-progress public renderer","CreationDate":None,"ModDate":None}
    fig.savefig(d/f"{stem}.pdf",bbox_inches="tight",pad_inches=.03,metadata=meta)
    fig.savefig(d/f"{stem}.png",dpi=300,bbox_inches="tight",pad_inches=.03,metadata={"Software":"eu27-material-climate-progress public renderer"})
    fig.savefig(d/f"{stem}.svg",bbox_inches="tight",pad_inches=.03,metadata={"Date":None,"Creator":"eu27-material-climate-progress public renderer"})
    plt.close(fig)


def _publication_labels(output_root: Path | None = None) -> dict:
    # The publication-label dictionary is a repository input, not an output-root artifact.
    # This keeps --output-root usable with a genuinely empty destination.
    d=pd.read_csv(ROOT/'docs'/'PUBLICATION_LABELS.csv',keep_default_na=False)
    return {r.internal_key:{'public':r.public_label,'short':r.short_public_label or r.public_label} for r in d.itertuples()}


def _pub(labels: dict, key: str, short: bool=True) -> str:
    if key not in labels:
        raise KeyError(f'Missing publication label: {key}')
    return labels[key]['short' if short else 'public']


def _adaptive_text(rgba) -> str:
    r,g,b=rgba[:3]
    lum=.2126*r+.7152*g+.0722*b
    return 'black' if lum>.58 else 'white'


def _rr_candidates():
    out=[]
    for r in [6,8,10,12,14,17,20,24,28,34]:
        for a in np.linspace(0,2*np.pi,16,endpoint=False): out.append((r*np.cos(a),r*np.sin(a)))
    return out


def _place_rank_labels(fig,ax,dd,font_pt=4.0):
    fig.canvas.draw(); ren=fig.canvas.get_renderer(); axbox=ax.get_window_extent(ren)
    marker=[]; placed=[]; rows=[]
    for _,r in dd.iterrows():
        px,py=ax.transData.transform((r.material_rank,r.climate_rank)); marker.append((r.geo,Bbox.from_extents(px-3.2,py-3.2,px+3.2,py+3.2)))
    lab=dd.assign(gap=(dd.material_rank-dd.climate_rank).abs()).sort_values(['gap','geo'],ascending=[False,True]); dims={}; tmp=[]
    for _,r in lab.iterrows(): tmp.append((r.geo,ax.text(r.material_rank,r.climate_rank,r.geo,fontsize=font_pt,ha='center',va='center',alpha=0)))
    fig.canvas.draw(); ren=fig.canvas.get_renderer()
    for c,t in tmp: b=t.get_window_extent(ren); dims[c]=(b.width*1.05,b.height*1.1); t.remove()
    overrides={
        ('CMUR_TERR','PL'):[(-8,-8)],('CMUR_TERR','LU'):[(-8,-7)],('DMC_TERR','LT'):[(8,7)],('DMC_TERR','HU'):[(-8,-8)],('DMC_TERR','DK'):[(8,-7)],('DMC_TERR','EL'):[(-8,7)],('DMC_TERR','NL'):[(-8,-7)],('RMC_TERR','HU'):[(-8,-8)],('RMC_TERR','BE'):[(-8,7)],('CMUR_CFOOT','LU'):[(-8,-7)],('DMC_CFOOT','HR'):[(8,7)],('DMC_CFOOT','HU'):[(-8,-8)],('DMC_CFOOT','NL'):[(-8,-7)],('DMC_CFOOT','PT'):[(-8,7)],('DMC_CFOOT','SE'):[(8,-7)],('RMC_CFOOT','ES'):[(-8,7)],('RMC_CFOOT','CY'):[(8,-7)]}
    panel=str(dd.panel.iloc[0]); texts=[]
    for _,r in lab.iterrows():
        anchor=np.array(ax.transData.transform((r.material_rank,r.climate_rank))); w,h=dims[r.geo]; chosen=None
        for dx,dy in overrides.get((panel,r.geo),[])+_rr_candidates():
            cp=anchor+np.array([dx,dy]); b=Bbox.from_bounds(cp[0]-w/2,cp[1]-h/2,w,h)
            if b.x0<axbox.x0+1 or b.x1>axbox.x1-1 or b.y0<axbox.y0+1 or b.y1>axbox.y1-1: continue
            if any(Bbox.overlaps(b,pb) for _,pb in placed): continue
            if any(Bbox.overlaps(b,mb) for _,mb in marker): continue
            chosen=(b,ax.transData.inverted().transform(cp),dx,dy); break
        if chosen is None: raise RuntimeError(f'Label placement failed: {panel}/{r.geo}')
        b,data,dx,dy=chosen
        t=ax.text(data[0],data[1],r.geo,fontsize=font_pt,ha='center',va='center',color=DARK,zorder=6,clip_on=True); t._geo=r.geo; texts.append(t); placed.append((r.geo,b))
        rows.append(dict(figure='figure_2_rank_concordance_geometry',panel=panel,country=r.geo,anchor_x=float(r.material_rank),anchor_y=float(r.climate_rank),label_center_x=float(data[0]),label_center_y=float(data[1]),distance_px=float(math.hypot(dx,dy)),label_vs_other_point=0,label_vs_other_label=0,outside_axes=0,status='PASS'))
    fig.canvas.draw(); ren=fig.canvas.get_renderer(); boxes={t._geo:t.get_window_extent(ren).expanded(1.01,1.04) for t in texts}
    for row in rows:
        b=boxes[row['country']]
        row['label_vs_other_point']=int(any(Bbox.overlaps(b,mb) for c,mb in marker if c!=row['country']))
        row['label_vs_other_label']=int(any(Bbox.overlaps(b,ob) for c,ob in boxes.items() if c!=row['country']))
        row['outside_axes']=int(not(axbox.contains(b.x0,b.y0) and axbox.contains(b.x1,b.y1)))
        if row['label_vs_other_point'] or row['label_vs_other_label'] or row['outside_axes']: row['status']='NEEDS_FIX'
    return rows


def _s2_label_offsets():
    # Directional micro-overrides in typographic points. Text alignment is chosen
    # from the sign of the offset, so a small offset keeps labels close to their
    # own sign-flip anchor rather than shifting the text center across borders.
    return {
        'CMUR → DMC':{
            'SE':(0,3),'FI':(0,-4),'EE':(6,1),'LV':(5,0),'LT':(6,-3),'DK':(-4,0),'ES':(0,4),'FR':(-5,0),
            # Manual Natural Earth label offsets preserve country ownership and legibility.
            'HU':(1,-1),'SI':(-3,0),'HR':(1.5,2),'BG':(6,0),'MT':(5,3)},
        'DMC → RMC':{
            'SE':(0,3),'DE':(0,4),
            'CZ':(-3,0),'SK':(3,0),'SI':(-4,-1),'LU':(-12,4)
        },
        'Territorial → consumption-based GHG':{
            'IE':(4,0),'LV':(5,0),'HU':(-1,-1),'RO':(3,0),'BG':(6,0)
        }
    }


def _render_s2_labels(fig,ax,geoms,d,flip,title):
    offsets=_s2_label_offsets();
    # Short leader lines improve semantic ownership in small-state and dense-country clusters.
    leaders={('CMUR → DMC','MT'),('CMUR → DMC','SI'),('CMUR → DMC','HR'),
             ('DMC → RMC','LU'),('DMC → RMC','SI')}
    q=d[d[flip].astype(bool)].copy().sort_values('geo'); marker_boxes={}; texts=[]; rows=[]
    for _,r0 in q.iterrows():
        p=geoms[r0.geo].representative_point(); px,py=ax.transData.transform((p.x,p.y)); marker_boxes[r0.geo]=Bbox.from_extents(px-3.0,py-3.0,px+3.0,py+3.0)
        ax.scatter([p.x],[p.y],marker='D',s=13,facecolors='none',edgecolors='black',linewidths=.62,zorder=6)
    for _,r0 in q.iterrows():
        c=r0.geo; p=geoms[c].representative_point(); dx,dy=offsets.get(title,{}).get(c,(4,3))
        ha='left' if dx>0 else ('right' if dx<0 else 'center'); va='bottom' if dy>0 else ('top' if dy<0 else 'center')
        leader=(title,c) in leaders
        t=ax.annotate(c,xy=(p.x,p.y),xytext=(dx,dy),textcoords='offset points',fontsize=4.0,ha=ha,va=va,color='#222222',zorder=8,clip_on=True)
        if leader:
            anchor_disp=np.array(ax.transData.transform((p.x,p.y))); text_anchor_disp=anchor_disp+np.array([dx,dy])*(fig.dpi/72.0); data_anchor=ax.transData.inverted().transform(text_anchor_disp)
            ax.plot([p.x,float(data_anchor[0])],[p.y,float(data_anchor[1])],color='#8F8F8F',lw=.22,zorder=5,solid_capstyle='round')
        t._geo=c; texts.append(t)
        rows.append(dict(figure='figure_s2_geography_representation_shifts',panel=title,country=c,marker_anchor_x=float(p.x),marker_anchor_y=float(p.y),text_bbox_x0_px=np.nan,text_bbox_y0_px=np.nan,text_bbox_x1_px=np.nan,text_bbox_y1_px=np.nan,label_center_x=np.nan,label_center_y=np.nan,distance_px=np.nan,leader='YES' if leader else 'NO',text_vs_own_marker=0,text_vs_other_marker=0,text_vs_other_text=0,manual_override='YES' if c in offsets.get(title,{}) else 'NO',human_semantic_ownership='PENDING',status='PASS'))
    fig.canvas.draw(); ren=fig.canvas.get_renderer(); axbox=ax.get_window_extent(ren); boxes={t._geo:t.get_window_extent(ren).expanded(1.01,1.04) for t in texts}
    for row in rows:
        c=row['country']; b=boxes[c]; row['text_bbox_x0_px']=float(b.x0); row['text_bbox_y0_px']=float(b.y0); row['text_bbox_x1_px']=float(b.x1); row['text_bbox_y1_px']=float(b.y1); center=np.array([(b.x0+b.x1)/2,(b.y0+b.y1)/2]); p=geoms[c].representative_point(); anchor=np.array(ax.transData.transform((p.x,p.y))); data=ax.transData.inverted().transform(center)
        row['label_center_x']=float(data[0]); row['label_center_y']=float(data[1]); row['distance_px']=float(np.linalg.norm(center-anchor))
        row['text_vs_own_marker']=int(Bbox.overlaps(b,marker_boxes[c]))
        row['text_vs_other_marker']=int(any(Bbox.overlaps(b,mb) for oc,mb in marker_boxes.items() if oc!=c))
        row['text_vs_other_text']=int(any(Bbox.overlaps(b,ob) for oc,ob in boxes.items() if oc!=c))
        outside=int(not(axbox.contains(b.x0,b.y0) and axbox.contains(b.x1,b.y1)))
        if row['text_vs_own_marker'] or row['text_vs_other_marker'] or row['text_vs_other_text'] or outside: row['status']='NEEDS_FIX'
    return rows


def render_all(output_root: Path, geoms: dict):
    sd=output_root/'figures'/'source_data'; extent=_map_extent(geoms); labels=_publication_labels(output_root); qa={'figure2_labels':[],'figure3_glyphs':[],'figure_s2_labels':[]}
    terr=_pub(labels,'TERR'); cfoot=_pub(labels,'CFOOT')

    # Figure 1 — selected compact matrix grammar: climate rows × material columns, sequential viridis.
    d=pd.read_csv(sd/'figure_1_primary_concordance_grid.csv'); mats=['CMUR','DMC','RMC']; clims=[terr,cfoot]
    arr=np.array([[float(d[(d.material_representation==m)&(d.climate_perspective==c)].spearman_rho.iloc[0]) for m in mats] for c in clims])
    fig,ax=plt.subplots(figsize=(W85,2.55)); cmap=plt.get_cmap('viridis'); norm=Normalize(0,.72)
    mesh=ax.pcolormesh(np.arange(4)-.5,np.arange(3)-.5,arr,cmap=cmap,norm=norm,shading='flat',edgecolors='white',linewidth=.65); ax.set_ylim(1.5,-.5)
    ax.set_xticks(range(3),mats); ax.set_yticks(range(2),clims)
    for i in range(2):
        for j in range(3):
            col=cmap(norm(arr[i,j])); ax.text(j,i,f'{arr[i,j]:.3f}',ha='center',va='center',fontsize=8.0,fontweight='bold',color=_adaptive_text(col))
    [sp.set_visible(False) for sp in ax.spines.values()]; ax.tick_params(length=0)
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,fraction=.055,pad=.04); cb.solids.set_rasterized(False); cb.set_label('Spearman ρ'); cb.set_ticks([0,.2,.4,.6,.7])
    fig.subplots_adjust(left=.28,right=.90,bottom=.16,top=.96); _save(fig,output_root,'figure_1_primary_concordance_grid')

    # Figure 2 — all-country rank-geometry design retained.
    d=pd.read_csv(sd/'figure_2_rank_concordance_geometry.csv'); fig,axes=plt.subplots(2,3,figsize=(W178,5.1)); order=['CMUR_TERR','DMC_TERR','RMC_TERR','CMUR_CFOOT','DMC_CFOOT','RMC_CFOOT']; letters='ABCDEF'
    for idx,(ax,panel) in enumerate(zip(axes.ravel(),order)):
        q=d[d.panel==panel].copy(); rho=float(q.spearman_rho.iloc[0]); ax.plot([1,27],[1,27],color='#BDBDBD',lw=.6,zorder=1); ax.scatter(q.material_rank,q.climate_rank,s=7,color='#4F5B66',zorder=2)
        ax.set_xlim(.1,27.9); ax.set_ylim(.1,27.9); ax.set_xticks([1,7,14,21,27]); ax.set_yticks([1,7,14,21,27]); ax.grid(color='#F1F1F1',lw=.35); ax.set_title(f'{q.panel_label.iloc[0]}\nρ = {rho:.3f}',fontsize=6.2,pad=4); ax.text(-.08,1.04,letters[idx],transform=ax.transAxes,fontsize=7,fontweight='bold',ha='left',va='bottom',clip_on=False)
        qa['figure2_labels'] += _place_rank_labels(fig,ax,q,4.0)
    fig.supxlabel('Material progress rank',fontsize=7); fig.supylabel('Climate progress rank',fontsize=7); fig.subplots_adjust(left=.065,right=.995,bottom=.08,top=.96,wspace=.14,hspace=.34); _save(fig,output_root,'figure_2_rank_concordance_geometry')

    # Figure 3 — selected country-consequence cell values with non-overlapping corner sign-disagreement glyphs.
    d=pd.read_csv(sd/'figure_3_country_representation_consequences.csv').sort_values('geo'); cols=['rank_shift_CMUR_to_DMC','rank_shift_DMC_to_RMC','rank_shift_TERR_to_CFOOT']; xlabels=['CMUR → DMC','DMC → RMC','Territorial →\nconsumption-based GHG']; flips=['sign_flip_CMUR_DMC','sign_flip_DMC_RMC','sign_flip_TERR_CFOOT']
    arr=d[cols].to_numpy(float); mx=max(abs(arr.min()),abs(arr.max())); fig,ax=plt.subplots(figsize=(W85,5.65)); norm=TwoSlopeNorm(vmin=-mx,vcenter=0,vmax=mx)
    mesh=ax.pcolormesh(np.arange(4)-.5,np.arange(len(d)+1)-.5,arr,cmap=DIVERGE,norm=norm,shading='flat',edgecolors='white',linewidth=.42); ax.set_ylim(len(d)-.5,-.5); ax.set_yticks(np.arange(len(d)),d.geo); ax.set_xticks(range(3),xlabels); ax.tick_params(axis='x',length=0,labelsize=5.25); ax.tick_params(axis='y',length=0)
    text_objs={}; glyph_boxes=[]
    for i,(_,r0) in enumerate(d.iterrows()):
        for j in range(3):
            rgba=DIVERGE(norm(arr[i,j])); t=ax.text(j,i,f'{arr[i,j]:+.0f}',ha='center',va='center',fontsize=5.15,color=_adaptive_text(rgba),zorder=5); text_objs[(i,j)]=t
            if bool(r0[flips[j]]):
                tri=MplPolygon([[j+.19,i-.46],[j+.46,i-.46],[j+.46,i-.19]],closed=True,facecolor='black',edgecolor='none',zorder=6); ax.add_patch(tri); glyph_boxes.append((i,j,tri))
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=DIVERGE),ax=ax,fraction=.032,pad=.025); cb.solids.set_rasterized(False); cb.set_label('Signed rank displacement')
    da=DrawingArea(10,10,0,0)
    da.add_artist(MplPolygon([[4,10],[10,10],[10,4]],closed=True,facecolor='black',edgecolor='none'))
    ab=AnnotationBbox(da,(.030,-.096),xycoords=ax.transAxes,frameon=False,box_alignment=(0,.5),pad=0,annotation_clip=False,zorder=7)
    ax.add_artist(ab)
    ax.text(.135,-.096,'corner mark = progress-sign disagreement',transform=ax.transAxes,fontsize=5.3,ha='left',va='center',clip_on=False)
    fig.subplots_adjust(left=.14,right=.88,bottom=.145,top=.985); fig.canvas.draw(); ren=fig.canvas.get_renderer()
    for i,j,tri in glyph_boxes:
        tb=text_objs[(i,j)].get_window_extent(ren).expanded(1.02,1.05); verts=ax.transData.transform(tri.get_xy()); gb=Bbox.from_extents(verts[:,0].min(),verts[:,1].min(),verts[:,0].max(),verts[:,1].max()); ov=int(Bbox.overlaps(tb,gb)); qa['figure3_glyphs'].append(dict(country=d.iloc[i].geo,substitution=xlabels[j],glyph_vs_number_overlap=ov,status='PASS' if ov==0 else 'NEEDS_FIX'))
    _save(fig,output_root,'figure_3_country_representation_consequences')

    # Figure 4 — selected 3+2 layout on frozen Natural Earth geometry.
    d=pd.read_csv(sd/'figure_4_geography_progress_ranks.csv'); rankcols=[('rank_CMUR','CMUR'),('rank_DMC','DMC'),('rank_RMC','RMC'),('rank_TERR',terr),('rank_CFOOT',cfoot)]; rn=Normalize(1,27); rcmap=plt.get_cmap('viridis_r')
    fig=plt.figure(figsize=(W178,5.05)); gs=fig.add_gridspec(4,6,height_ratios=[.14,1,.18,1],left=.002,right=.998,bottom=.105,top=.982,wspace=.002,hspace=.010)
    h1=fig.add_subplot(gs[0,:]); h1.axis('off'); h1.text(.5,.58,'Material representations',ha='center',va='center',fontsize=6.6,fontweight='semibold')
    axs=[fig.add_subplot(gs[1,0:2]),fig.add_subplot(gs[1,2:4]),fig.add_subplot(gs[1,4:6]),fig.add_subplot(gs[3,1:3]),fig.add_subplot(gs[3,3:5])]
    h2=fig.add_subplot(gs[2,:]); h2.axis('off'); h2.text(.5,.72,'Climate perspectives',ha='center',va='center',fontsize=6.6,fontweight='semibold')
    for ax,(col,title) in zip(axs,rankcols): _draw_map(ax,geoms,dict(zip(d.geo,d[col])),rcmap,rn,extent); ax.set_title(title if len(title)<24 else title.replace(' GHG footprint',' GHG\nfootprint'),fontsize=6.7,pad=.55)
    cax=fig.add_axes([.30,.032,.40,.025]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=rn,cmap=rcmap),cax=cax,orientation='horizontal'); cb.solids.set_rasterized(False); cb.set_ticks([1,7,14,21,27]); cb.set_label('Progress rank (1 = most favorable)')
    fig.canvas.draw()
    for a,b in [(axs[0],axs[1]),(axs[1],axs[2]),(axs[3],axs[4])]:
        x=(a.get_position().x1+b.get_position().x0)/2; y0=max(a.get_position().y0,b.get_position().y0); y1=min(a.get_position().y1,b.get_position().y1); fig.add_artist(Line2D([x,x],[y0,y1],transform=fig.transFigure,color='#C8C8C8',lw=.55,zorder=0))
    _save(fig,output_root,'figure_4_geography_progress_ranks')

    # Figure 5 — temporal line-profile design with publication labels and group hierarchy.
    d=pd.read_csv(sd/'figure_5_temporal_representation_sensitivity.csv').set_index('period').loc[PERIOD_ORDER].reset_index(); x=np.arange(6); cols=['a_CMUR_TERR','b_DMC_TERR','c_RMC_TERR','d_CMUR_CFOOT','e_DMC_CFOOT','f_RMC_CFOOT']; xlab=['CMUR','DMC','RMC','CMUR','DMC','RMC']
    fig,ax=plt.subplots(figsize=(W178,3.15)); styles=[(OKABE_BLUE,'D','-',1.65),(OKABE_ORANGE,'o','--',1.25),('#009E73','^','-.',1.25)]
    for (_,r0),(color,marker,ls,lw) in zip(d.iterrows(),styles): ax.plot(x,[r0[c] for c in cols],marker=marker,ls=ls,lw=lw,color=color,label=PERIOD_LABEL[r0.period],markersize=4)
    ax.axvline(2.5,color='#BBBBBB',lw=.65); ax.axhline(0,color='#8D8D8D',lw=.6); ax.set_xticks(x,xlab); ax.set_ylabel('Spearman ρ'); ax.set_ylim(-.28,.78); ax.grid(axis='y',color=GRID,lw=.4); ax.legend(frameon=False,ncol=3,loc='upper center',bbox_to_anchor=(.5,1.13)); [ax.spines[s].set_visible(False) for s in ['top','right']]
    ax.text(1,-.18,terr,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=6.2,fontweight='semibold'); ax.text(4,-.18,cfoot,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=6.2,fontweight='semibold')
    fig.subplots_adjust(left=.08,right=.99,bottom=.27,top=.84); _save(fig,output_root,'figure_5_temporal_representation_sensitivity')

    # Figure 6 — selected deletion-cloud design, public semantic row labels.
    d=pd.read_csv(sd/'figure_6_deletion_robustness.csv'); estim=['D1_b_minus_a','D2_c_minus_b','D3_f_minus_c','D_TOTAL_f_minus_a']; labs=['CMUR → DMC','DMC → RMC','Territorial → consumption-based GHG','Overall representation contrast']; fig,ax=plt.subplots(figsize=(W178,2.75)); rng=np.random.default_rng(1947)
    for i,e in enumerate(estim):
        dd=d[d.estimand==e]
        for source,color,marker in [('LOCO',OKABE_BLUE,'o'),('LOYO',OKABE_ORANGE,'s')]:
            vv=dd[dd.source==source].value.to_numpy(); yy=i+(rng.random(len(vv))-.5)*.20; ax.scatter(vv,yy,s=11,facecolors='none' if source=='LOYO' else color,edgecolors=color,lw=.6,marker=marker,alpha=.86,label=source if i==0 else None,zorder=2)
        full=float(dd[dd.source=='FULL'].value.iloc[0]); ax.scatter([full],[i],s=16,marker='D',facecolors='none',edgecolors='black',linewidths=.72,zorder=7,label='Primary 2010–2023' if i==0 else None)
    ax.axvline(0,color='#777777',lw=.7); ax.set_yticks(range(4),labs,fontsize=6.4); ax.invert_yaxis(); ax.set_xlabel('Difference in Spearman ρ'); ax.legend(frameon=False,ncol=3,fontsize=5.8,loc='upper center',bbox_to_anchor=(.5,1.16)); [ax.spines[s].set_visible(False) for s in ['top','right','left']]; ax.grid(axis='x',color=GRID,lw=.35,zorder=0); fig.subplots_adjust(left=.23,right=.985,bottom=.20,top=.80); _save(fig,output_root,'figure_6_deletion_robustness')

    # Figure S1 — selected temporal matrices; diverging scale is scientifically warranted.
    d=pd.read_csv(sd/'figure_s1_temporal_concordance_matrices.csv').set_index('period').loc[PERIOD_ORDER]; vals=[]
    for _,r in d.iterrows(): vals += [r.a_CMUR_TERR,r.b_DMC_TERR,r.c_RMC_TERR,r.d_CMUR_CFOOT,r.e_DMC_CFOOT,r.f_RMC_CFOOT]
    norm=TwoSlopeNorm(vmin=min(vals),vcenter=0,vmax=max(vals)); fig=plt.figure(figsize=(W178,2.72)); gs=fig.add_gridspec(1,4,width_ratios=[1,1,1,.055],wspace=.16); titles=['2010–2023\nPrimary','2010–2016\nSecondary','2017–2023\nSecondary']
    for k,(_,r) in enumerate(d.iterrows()):
        ax=fig.add_subplot(gs[0,k]); arr=np.array([[r.a_CMUR_TERR,r.b_DMC_TERR,r.c_RMC_TERR],[r.d_CMUR_CFOOT,r.e_DMC_CFOOT,r.f_RMC_CFOOT]]); ax.set_xlim(-.5,2.5); ax.set_ylim(1.5,-.5)
        for i in range(2):
            for j in range(3):
                col=DIVERGE(norm(arr[i,j])); ax.add_patch(plt.Rectangle((j-.5,i-.5),1,1,facecolor=col,edgecolor='white',lw=.7)); ax.text(j,i,f'{arr[i,j]:.3f}',ha='center',va='center',fontsize=6.2,color=_adaptive_text(col),fontweight='semibold')
        ax.set_xticks(range(3),['CMUR','DMC','RMC']); ax.set_title(titles[k],fontsize=6.5,pad=5); ax.tick_params(length=0)
        if k==0: ax.set_yticks([0,1],[terr,cfoot],fontsize=5.8)
        else: ax.set_yticks([])
        [sp.set_visible(False) for sp in ax.spines.values()]
    cax=fig.add_subplot(gs[0,3]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=DIVERGE),cax=cax); cb.solids.set_rasterized(False); cb.set_label('Spearman ρ'); fig.subplots_adjust(left=.15,right=.97,bottom=.16,top=.86); _save(fig,output_root,'figure_s1_temporal_concordance_matrices')

    # Figure S2 — Natural Earth country-label adaptation with explicit semantic label policy.
    d=pd.read_csv(sd/'figure_s2_geography_representation_shifts.csv'); panels=[('rank_shift_CMUR_to_DMC','sign_flip_CMUR_DMC','CMUR → DMC'),('rank_shift_DMC_to_RMC','sign_flip_DMC_RMC','DMC → RMC'),('rank_shift_TERR_to_CFOOT','sign_flip_TERR_CFOOT','Territorial → consumption-based GHG')]; mx=max(abs(d[[p[0] for p in panels]].to_numpy(float).min()),abs(d[[p[0] for p in panels]].to_numpy(float).max())); norm=TwoSlopeNorm(vmin=-mx,vcenter=0,vmax=mx)
    fig,axes=plt.subplots(1,3,figsize=(W178,2.75),gridspec_kw={'wspace':.025})
    fig.subplots_adjust(left=.005,right=.995,bottom=.18,top=.94)
    for ax,(col,flip,title) in zip(axes,panels):
        _draw_map(ax,geoms,dict(zip(d.geo,d[col])),DIVERGE,norm,extent); ax.set_title(title,fontsize=6.5,pad=2); qa['figure_s2_labels'] += _render_s2_labels(fig,ax,geoms,d,flip,title)
    cax=fig.add_axes([.30,.06,.40,.025]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=DIVERGE),cax=cax,orientation='horizontal'); cb.solids.set_rasterized(False); cb.set_label('Signed rank displacement'); fig.text(.72,.073,'◇ progress-sign disagreement',fontsize=5.3,ha='left',va='center')
    fig.canvas.draw()
    for a,b in [(axes[0],axes[1]),(axes[1],axes[2])]:
        x=(a.get_position().x1+b.get_position().x0)/2; fig.add_artist(Line2D([x,x],[.20,.91],transform=fig.transFigure,color='#C8C8C8',lw=.5,zorder=0))
    _save(fig,output_root,'figure_s2_geography_representation_shifts')

    # Figure S3 — selected layout discipline on Natural Earth, separate CMUR and common log-rate scales.
    d=pd.read_csv(sd/'figure_s3_geography_progress_scores.csv'); cols=[('CMUR','CMUR'),('DMC','DMC'),('RMC','RMC'),('TERR',terr),('CFOOT',cfoot)]; others=d[['DMC','RMC','TERR','CFOOT']].to_numpy(float); lim=max(abs(others.min()),abs(others.max())); cmur_lim=max(abs(d.CMUR.min()),abs(d.CMUR.max())); lnorm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim); cnorm=TwoSlopeNorm(vmin=min(d.CMUR.min(),0),vcenter=0,vmax=max(d.CMUR.max(),0))
    fig=plt.figure(figsize=(W178,5.10)); gs=fig.add_gridspec(4,6,height_ratios=[.14,1,.18,1],left=.002,right=.998,bottom=.117,top=.982,wspace=.002,hspace=.010); h1=fig.add_subplot(gs[0,:]); h1.axis('off'); h1.text(.5,.58,'Material representations',ha='center',va='center',fontsize=6.6,fontweight='semibold'); axs=[fig.add_subplot(gs[1,0:2]),fig.add_subplot(gs[1,2:4]),fig.add_subplot(gs[1,4:6]),fig.add_subplot(gs[3,1:3]),fig.add_subplot(gs[3,3:5])]; h2=fig.add_subplot(gs[2,:]); h2.axis('off'); h2.text(.5,.72,'Climate perspectives',ha='center',va='center',fontsize=6.6,fontweight='semibold')
    for i,(ax,(col,title)) in enumerate(zip(axs,cols)):
        normi=cnorm if i==0 else lnorm; _draw_map(ax,geoms,dict(zip(d.geo,d[col])),DIVERGE,normi,extent); ax.set_title(title if len(title)<24 else title.replace(' GHG footprint',' GHG\nfootprint'),fontsize=6.7,pad=.55)
    c1=fig.add_axes([.09,.037,.31,.024]); cb1=fig.colorbar(plt.cm.ScalarMappable(norm=cnorm,cmap=DIVERGE),cax=c1,orientation='horizontal'); cb1.solids.set_rasterized(False); cb1.set_label('CMUR pp/year')
    c2=fig.add_axes([.54,.037,.37,.024]); cb2=fig.colorbar(plt.cm.ScalarMappable(norm=lnorm,cmap=DIVERGE),cax=c2,orientation='horizontal'); cb2.solids.set_rasterized(False); cb2.set_label('Favorable log-rate/year (DMC/RMC/GHG)')
    fig.canvas.draw()
    for a,b in [(axs[0],axs[1]),(axs[1],axs[2]),(axs[3],axs[4])]:
        x=(a.get_position().x1+b.get_position().x0)/2; y0=max(a.get_position().y0,b.get_position().y0); y1=min(a.get_position().y1,b.get_position().y1); fig.add_artist(Line2D([x,x],[y0,y1],transform=fig.transFigure,color='#C8C8C8',lw=.55,zorder=0))
    _save(fig,output_root,'figure_s3_geography_progress_scores')

    # Figure S4 — progress-score geometry, reader-facing panel labels, no ISO labels.
    d=pd.read_csv(sd/'figure_s4_progress_score_geometry.csv'); order=['CMUR_TERR','DMC_TERR','RMC_TERR','CMUR_CFOOT','DMC_CFOOT','RMC_CFOOT']; fig,axes=plt.subplots(2,3,figsize=(W178,4.9)); letters='ABCDEF'; panel_labels={'CMUR_TERR':f'CMUR × {terr}','DMC_TERR':f'DMC × {terr}','RMC_TERR':f'RMC × {terr}','CMUR_CFOOT':f'CMUR × {cfoot}','DMC_CFOOT':f'DMC × {cfoot}','RMC_CFOOT':f'RMC × {cfoot}'}
    for idx,(ax,panel) in enumerate(zip(axes.ravel(),order)):
        q=d[d.panel==panel]; ax.scatter(q.material_progress_score,q.climate_progress_score,s=8,color=OKABE_BLUE,alpha=.78); ax.axhline(0,color='#BBBBBB',lw=.5); ax.axvline(0,color='#BBBBBB',lw=.5); ax.grid(color='#F0F0F0',lw=.3); ax.set_title(f'{panel_labels[panel]}\nρ = {q.spearman_rho.iloc[0]:.3f} | r = {q.pearson_r.iloc[0]:.3f}',fontsize=5.8,pad=3); ax.text(-.08,1.06,letters[idx],transform=ax.transAxes,fontsize=7,fontweight='bold',ha='left',va='bottom',clip_on=False); [ax.spines[s].set_visible(False) for s in ['top','right']]
    fig.supxlabel('Material progress score'); fig.supylabel('Climate progress score'); fig.subplots_adjust(left=.08,right=.99,bottom=.09,top=.96,wspace=.25,hspace=.32); _save(fig,output_root,'figure_s4_progress_score_geometry')

    # Figure S5 — robust structure retained, labels controlled by publication dictionary.
    d=pd.read_csv(sd/'figure_s5_broader_robustness_summary.csv').drop_duplicates('diagnostic',keep='last'); order=['primary_FULL','LOCO','LOYO','CMUR_even_year','published_pc','common_log_CMUR','Theil_Sen','total_scale']; d=d.set_index('diagnostic').reindex(order).dropna(how='all').reset_index(); label_map={'primary_FULL':'Primary 2010–2023','LOCO':_pub(labels,'LOCO',False),'LOYO':_pub(labels,'LOYO',False),'CMUR_even_year':_pub(labels,'CMUR_even_year',False),'published_pc':_pub(labels,'published_pc',False),'common_log_CMUR':_pub(labels,'common_log_CMUR',False),'Theil_Sen':_pub(labels,'Theil_Sen',False),'total_scale':_pub(labels,'total_scale',False)}
    fig,ax=plt.subplots(figsize=(W85,3.65)); y=np.arange(len(d)); primary=float(d[d.diagnostic=='primary_FULL']['min'].iloc[0])
    for i,r0 in d.iterrows():
        if abs(r0['max']-r0['min'])<1e-14: ax.scatter([r0['min']],[i],s=18,facecolors='white',edgecolors=OKABE_BLUE,zorder=3)
        else: ax.plot([r0['min'],r0['max']],[i,i],color=OKABE_BLUE,lw=1.35,marker='o',markersize=2.9)
    ax.axvline(primary,color=OKABE_ORANGE,lw=.9,ls='--'); ax.text(primary,1.015,f'Primary 2010–2023 = {primary:.3f}',transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=5.8,color=OKABE_ORANGE,clip_on=False); ax.set_yticks(y,[label_map[x] for x in d.diagnostic]); ax.invert_yaxis(); ax.set_xlabel('Overall representation contrast'); ax.grid(axis='x',color=GRID,lw=.4); [ax.spines[s].set_visible(False) for s in ['top','right','left']]; fig.subplots_adjust(left=.42,right=.98,bottom=.14,top=.91); _save(fig,output_root,'figure_s5_broader_robustness_summary')

    # Figure S6 — categorical UpSet structure retained; no artificial numeric x-axis.
    d=pd.read_csv(sd/'figure_s6_sign_flip_overlap.csv').reset_index(drop=True); x=np.arange(len(d)); fig=plt.figure(figsize=(W85,3.3)); gs=fig.add_gridspec(2,1,height_ratios=[1.25,1],hspace=.05); axb=fig.add_subplot(gs[0]); axm=fig.add_subplot(gs[1],sharex=axb); bars=axb.bar(x,d.n_countries,color=OKABE_BLUE,width=.60)
    for rect,val in zip(bars,d.n_countries): axb.text(rect.get_x()+rect.get_width()/2,rect.get_height()+.08,str(int(val)),ha='center',va='bottom',fontsize=5.8)
    axb.set_ylabel('Countries'); axb.set_xticks([]); axb.set_ylim(0,max(d.n_countries)+1.25); axb.grid(axis='y',color=GRID,lw=.35); [axb.spines[s].set_visible(False) for s in ['top','right','bottom']]
    rows=[('CMUR → DMC','sign_flip_CMUR_DMC'),('DMC → RMC','sign_flip_DMC_RMC'),('Territorial → consumption-based GHG','sign_flip_TERR_CFOOT')]
    for yi,(lab,col) in enumerate(rows):
        xs=[]
        for xi,val in enumerate(d[col].astype(bool)): axm.scatter([xi],[yi],s=21,facecolors='#222222' if val else 'white',edgecolors='#777777',linewidths=.55); xs.append(xi if val else None)
        active=[xi for xi,val in enumerate(d[col].astype(bool)) if val]
    for xi in range(len(d)):
        active=[yi for yi,(_,col) in enumerate(rows) if bool(d.loc[xi,col])]
        if len(active)>1: axm.plot([xi,xi],[min(active),max(active)],color='#222222',lw=.7,zorder=0)
    axm.set_yticks(range(3),[r[0] for r in rows]); axm.set_xticks(x,['' for _ in x]); axm.tick_params(axis='x',length=0); axm.set_xlabel('Observed sign-flip combinations'); axm.set_ylim(2.6,-.6); [axm.spines[s].set_visible(False) for s in ['top','right','left']]; fig.subplots_adjust(left=.32,right=.98,bottom=.16,top=.95); _save(fig,output_root,'figure_s6_sign_flip_overlap')

    # Figure S7 — CD01 grouping, public labels and subtle cell boundaries.
    d=pd.read_csv(sd/'figure_s7_country_discordance_matrix.csv'); countries=sorted(d.geo.unique()); cells=['CMUR×TERR','DMC×TERR','RMC×TERR','CMUR×CFOOT','DMC×CFOOT','RMC×CFOOT']; arr=np.array([[float(d[(d.geo==g)&(d.cell==c)].material_minus_climate_rank.iloc[0]) for c in cells] for g in countries]); mx=max(abs(arr.min()),abs(arr.max())); fig,ax=plt.subplots(figsize=(W178,4.75)); norm=TwoSlopeNorm(vmin=-mx,vcenter=0,vmax=mx); mesh=ax.pcolormesh(np.arange(7)-.5,np.arange(len(countries)+1)-.5,arr,cmap=DIVERGE,norm=norm,shading='flat',edgecolors='#E8E8E8',linewidth=.32); ax.set_ylim(len(countries)-.5,-.5); ax.set_yticks(range(len(countries)),countries); ax.set_xticks(range(6),['CMUR','DMC','RMC','CMUR','DMC','RMC']); ax.axvline(2.5,color='#9D9D9D',lw=.75); ax.text(1,1.02,terr,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=6.4,fontweight='semibold',clip_on=False); ax.text(4,1.02,cfoot,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=6.4,fontweight='semibold',clip_on=False); cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=DIVERGE),ax=ax,fraction=.025,pad=.02); cb.solids.set_rasterized(False); cb.set_label('Material rank − climate rank'); fig.subplots_adjust(left=.06,right=.93,bottom=.08,top=.91); _save(fig,output_root,'figure_s7_country_discordance_matrix')
    return qa

def build_all(analysis: Path, output_root: Path, ne_zip: Path):
    build_publication_sources(analysis, output_root)
    _,geoms=build_country_code_map(ne_zip, output_root)
    return render_all(output_root, geoms)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--analysis",type=Path,default=ROOT/"results"/"reference",help="Validated analysis-output root (tables/, supplementary/, source_data/)")
    ap.add_argument("--output-root",type=Path,default=ROOT,help="Root under which publication paths are written")
    ap.add_argument("--natural-earth",type=Path,default=ROOT/"cartography"/"source"/"ne_10m_admin_0_countries.zip")
    args=ap.parse_args()
    build_all(args.analysis,args.output_root,args.natural_earth)
    print("PUBLICATION_LAYER_BUILD_PASS figures=13 main_tables=4 supplementary_tables=8 map_join=27/27")

if __name__=="__main__":
    main()
