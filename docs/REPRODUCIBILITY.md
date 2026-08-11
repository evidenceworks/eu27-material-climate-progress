# Reproducibility notes

A complete local reproduction uses the frozen source responses bundled with the repository; internet access is not required after dependencies are installed.

Recommended sequence:

```bash
python code/checks.py verify
python code/01_prepare_data.py
python tests/test_ties.py
python code/02_analyze_progress.py
python tests/test_reproduction.py
python code/04_render_publication_figures.py --analysis work/analysis --output-root work/publication
python tests/test_publication_tables.py
python tests/test_publication_figure_source_data.py
python tests/test_publication_artwork.py
python tests/test_publication_labels.py
python tests/test_map_join_parity.py
```

The publication renderer may also be run against `results/reference/` to regenerate the checked-in publication layer from the validated numerical reference files.

The Natural Earth source archive is frozen in the repository, so map reproduction does not depend on a live cartographic download.
