# CHB gut microbial fermentation pathway analysis

This repository contains the code, public-data inputs, and derived outputs for a reproducible metagenomic pathway analysis of biopsy-stratified chronic hepatitis B (CHB) with two external public-data comparisons.

The analysis asks whether a selected fermentation-related MetaCyc pathway composite differs by biopsy-defined histologic injury in CHB, and whether the same finalized pathway definitions show concordant behaviour in public NAFLD fibrosis and cirrhosis datasets.

## Contents

- `analysis/module_definitions.tsv`: MetaCyc pathway definitions used in all cohorts.
- `data/prjdb36442/`: CHB pathway abundance matrix and sample manifest.
- `data/loombar2017/`: processed LoombaR_2017 pathway abundance matrix and metadata.
- `data/qinn2014/`: processed QinN_2014 pathway abundance matrix and metadata.
- `scripts/`: analysis and figure/table scripts.
- `outputs/`: rebuilt figures, tables, and analysis results.

The repository uses processed pathway-abundance matrices rather than raw sequencing reads. This keeps the analysis lightweight and reproducible on a standard computer.

## Reproduce the analysis

Install Python dependencies:

```bash
python3 -m pip install -r requirements.txt
```

Run the full workflow:

```bash
bash scripts/run_all.sh
```

The command rebuilds:

- `outputs/results/`: cohort-specific module scores, effect estimates, leave-one-out analyses, pathway-member summaries, and random-module calibrations.
- `outputs/figures/`: main figures and supplementary figure.
- `outputs/tables/`: main tables and supplementary workbook.

## Data sources

- PRJDB36442: public CHB shotgun metagenomic pathway profiles and sample grouping reported by the source study.
- LoombaR_2017 / PRJNA373901: processed pathway abundance and metadata from curatedMetagenomicData.
- QinN_2014 / PRJEB6337 / ERP005860: processed pathway abundance and metadata from curatedMetagenomicData.

## Interpretation

The CHB cohort is a discovery analysis based on a small biopsy-stratified public dataset. The NAFLD comparison is unadjusted because age category is fully collinear with fibrosis group in the exported metadata and sex, BMI, and diabetes are unavailable. The cirrhosis cohort is used as contextual evidence rather than as a direct histologic replication.

## License

Code is released under the MIT License. Public datasets remain governed by their original data-source terms.
