# Method and output map

## Computational chain

1. `data/raw/` stores the frozen Eurostat JSON-stat source responses.
2. `code/01_prepare_data.py` reconstructs the EU27 × 2010-2023 panel.
3. `code/02_analyze_progress.py` reproduces the validated progress scores, ranks, correlation grids, ordered contrasts, country consequences, temporal results, robustness analyses, accounting bridge, RMC status summary, and tie audit.
4. `code/03_build_results.py` retains the neutral numerical/result-package layer from the validated pipeline.
5. `code/04_render_publication_figures.py` is the publication layer. It reorganizes analysis outputs into the manuscript/Supplement table sources and figure source data, joins those values to Natural Earth geometry for the three map figures, and renders publication artwork.

Step 5 performs **no estimation** of progress scores, ranks, correlations, ordered contrasts, or robustness statistics.

## Numerical reference outputs

The validated computational reference files are under `results/reference/`. The detailed outputs include all-period/all-spec country results, representation grids, deletion analyses, denominator/scale diagnostics, the accounting bridge, RMC status information, and tie audit.

## Publication layer

- MAIN Tables 2-5: `results/reference/tables/table_2_...` through `table_5_...`
- Supplement Tables S1-S8: `supplementary_tables/`
- Figure source data: `figures/source_data/`
- Vector publication artwork: `figures/reference/*.pdf`
- PNG previews: `figures/reference/*.png`
- SVG vector companions: `figures/reference/*.svg`

The exact manuscript/Supplement mapping is in `PUBLICATION_OUTPUT_MAP.csv`.

## Determinism and tolerance

- Reconstructed scientific CSVs are compared with the validated numerical reference at `1e-12`.
- Panel reconstruction is compared at `1e-6` because some source values are serialized before reconstruction.
- Publication-source CSVs are regenerated and compared exactly at the parsed value/schema level.
- PNG reproduction is required to match dimensions and a defined pixel-difference tolerance; vector PDF/SVG files are checked structurally rather than by byte identity because backend metadata/font object ordering can vary by platform.
- Map join keys are exact and must remain 27/27 before map rendering.
