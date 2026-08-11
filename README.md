# EU27 material-climate progress

This repository reproduces a longitudinal comparison of national material and greenhouse-gas progress across the EU27 for 2010–2023. It includes the fixed Eurostat source responses used in the analysis, the analysis code, reference results, source tables and the renderer that generates the figures used in the manuscript and Supplementary Information.


## Run in Google Colab

The first notebook rebuilds and packages the analysis. The second uses that package to reproduce the publication tables and figures from the same Git commit.

- Prepare and package the analysis: [![Open notebook 1 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/evidenceworks/eu27-material-climate-progress/blob/main/notebooks/01_prepare_material_climate_analysis.ipynb)
- Reproduce the final publication layer: [![Open notebook 2 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/evidenceworks/eu27-material-climate-progress/blob/main/notebooks/02_reproduce_material_climate_results.ipynb)

## Scientific scope

The analysis compares progress rankings across five official representations:

- circular material use rate (CMUR);
- domestic material consumption (DMC);
- raw material consumption / material footprint (RMC);
- territorial greenhouse-gas emissions;
- consumption-based greenhouse-gas footprint.

The primary analysis period is Full (2010–2023). Early (2010–2016) and Late (2017–2023) are secondary analysis windows. The repository also contains Kendall/Pearson companions, country and year deletion checks, estimator/scale/denominator sensitivities, the production/residence accounting bridge, RMC status provenance, and the deterministic tie audit.

## Reproduce the analysis

Create an environment with Python 3.11+ and install the pinned dependencies:

```bash
python -m pip install -r requirements-lock.txt
```

Then run:

```bash
python code/checks.py verify
python code/01_prepare_data.py
python tests/test_ties.py
python code/02_analyze_progress.py
python tests/test_reproduction.py
```

The reproduced analysis is written to `work/analysis/` and is checked against the reference numerical outputs at tolerance `1e-12` (panel reconstruction uses `1e-6`, reflecting serialization of source values).

## Reproduce publication sources and figures

After the analysis step:

```bash
python code/04_render_publication_figures.py \
  --analysis work/analysis \
  --output-root work/publication \
  --natural-earth cartography/source/ne_10m_admin_0_countries.zip
```

Verify that the reproduced publication layer matches the checked-in reference publication files under the documented reproducibility rules (parsed schema/text semantics, numerical tolerance `1e-12`, pixel-identical PNG previews and structural vector checks), then package it:

```bash
python code/checks.py compare-publication --publication work/publication
python code/checks.py package-publication \
  --publication work/publication \
  --output eu27-material-climate-publication-results.zip
```

This writes the publication layer under `work/publication/`. The checked-in reference publication files are under:

- `results/reference/tables/` - manuscript Tables 2-5 quantitative/factual sources;
- `supplementary_tables/` - Tables S1-S8 sources;
- `figures/source_data/` - source data for Figures 1-6 and S1-S7;
- `figures/reference/` - publication artwork in PDF/PNG/SVG formats.

Run the publication tests with:

```bash
python tests/test_publication_tables.py
python tests/test_publication_figure_source_data.py
python tests/test_publication_artwork.py
python tests/test_publication_labels.py
python tests/test_map_join_parity.py
python tests/test_publication_renderer_empty_output_root.py
python tests/test_notebook_publication_contract.py
```

`docs/PUBLICATION_OUTPUT_MAP.csv` maps each manuscript/Supplement object to its repository source. `docs/PUBLICATION_LABELS.csv` provides the publication label dictionary.

## Cartography

Publication maps use Natural Earth Admin 0 Countries, 1:10m, version 5.1.1. The exact source archive is frozen under `cartography/source/`. Scientific country keys are joined to map geometry through `cartography/COUNTRY_CODE_MAP.csv`; `MAP_JOIN_PARITY.csv` records the 27/27 join check. Eurostat `EL` is explicitly mapped to Natural Earth `GR`, and France is matched to the sovereign-country record rather than dependency records.

Map geometry is presentation support only. Changing boundary geometry does not change progress scores, ranks, correlations, ordered contrasts, sign categories, or country categories.

## Data and licensing

See `DATA.md`, `DATA_AVAILABILITY.md`, and `THIRD_PARTY_DATA_AND_LICENSES.md`. The MIT License covers author-created software only; third-party/source data retain their own reuse terms.

## Repository structure

See `docs/METHOD_AND_OUTPUT_MAP.md` for the relationship among source snapshots, analysis outputs, publication tables, figure source data, artwork, and verification tests.
