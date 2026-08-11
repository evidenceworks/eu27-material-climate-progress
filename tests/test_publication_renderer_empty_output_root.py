#!/usr/bin/env python3
from pathlib import Path
import hashlib, shutil, subprocess, sys, tempfile
ROOT=Path(__file__).resolve().parents[1]
STEMS=[
"figure_1_primary_concordance_grid","figure_2_rank_concordance_geometry","figure_3_country_representation_consequences","figure_4_geography_progress_ranks","figure_5_temporal_representation_sensitivity","figure_6_deletion_robustness","figure_s1_temporal_concordance_matrices","figure_s2_geography_representation_shifts","figure_s3_geography_progress_scores","figure_s4_progress_score_geometry","figure_s5_broader_robustness_summary","figure_s6_sign_flip_overlap","figure_s7_country_discordance_matrix"]
MAIN=["table_2_representation_contract_source.csv","table_3_primary_concordance_and_contrasts.csv","table_4_country_consequences.csv","table_5_temporal_and_robustness.csv"]
SUPP=["table_s1_sources_units_provenance.csv","table_s2_temporal_concordance_matrices.csv","table_s3_companion_correlations.csv","table_s4_temporal_ordered_contrasts.csv","table_s5_extended_robustness.csv","table_s6_production_account_bridge.csv","table_s7_rmc_status_provenance.csv","table_s8_numeric_scale_diagnostics.csv"]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    td=Path(tempfile.mkdtemp(prefix="publication_empty_root_")); out=td/"publication"
    try:
        assert not out.exists()
        subprocess.run([sys.executable,str(ROOT/"code"/"04_render_publication_figures.py"),"--analysis",str(ROOT/"results"/"reference"),"--output-root",str(out),"--natural-earth",str(ROOT/"cartography"/"source"/"ne_10m_admin_0_countries.zip")],check=True,cwd=ROOT)
        rels=[Path("results/reference/tables")/x for x in MAIN]+[Path("supplementary_tables")/x for x in SUPP]+[Path("figures/source_data")/f"{x}.csv" for x in STEMS]+[Path("figures/reference")/f"{x}.{e}" for x in STEMS for e in ["pdf","png","svg"]]+[Path("MAP_JOIN_PARITY.csv")]
        assert len(rels)==65
        for rel in rels:
            a=out/rel; b=ROOT/rel
            assert a.is_file(),f"missing {rel}"
            assert b.is_file(),f"reference missing {rel}"
            assert sha(a)==sha(b),f"publication parity mismatch: {rel}"
        print("PUBLICATION_EMPTY_OUTPUT_ROOT_TEST_PASS files=65 differences=0")
    finally: shutil.rmtree(td,ignore_errors=True)
if __name__=="__main__": main()
