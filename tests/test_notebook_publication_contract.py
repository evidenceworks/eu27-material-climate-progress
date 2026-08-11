#!/usr/bin/env python3
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=ROOT/"notebooks"/"02_reproduce_material_climate_results.ipynb"
    nb=json.loads(p.read_text(encoding="utf-8"))
    text="\n".join("".join(c.get("source",[])) for c in nb.get("cells",[]))
    required=[
        "code/04_render_publication_figures.py",
        "--output-root",
        "work/publication",
        "compare-analysis",
        "compare-publication",
        "package-publication",
        "eu27-material-climate-publication-results.zip",
        'receipt["git_commit"]',
        "repository_manifest_sha256",
    ]
    for token in required:
        assert token in text,f"Notebook 2 missing contract token: {token}"
    assert "code/03_build_results.py" not in text
    assert "package-results" not in text
    assert "eu27-material-climate-results.zip" not in text
    assert text.count("files.download(")==1,"Notebook 2 must download exactly one final publication-results ZIP"
    for c in nb.get("cells",[]):
        if c.get("cell_type")=="code":
            assert not c.get("outputs"),"Notebook contains saved outputs"
            assert c.get("execution_count") is None,"Notebook contains saved execution count"
    print("NOTEBOOK_PUBLICATION_CONTRACT_TEST_PASS renderer=final publication_zip=single_download saved_outputs=0")

if __name__=="__main__":
    main()
