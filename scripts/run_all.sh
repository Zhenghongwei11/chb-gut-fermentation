#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

rm -rf outputs/results outputs/figures outputs/tables

mkdir -p outputs/results/prjdb36442/conservative \
  outputs/results/prjdb36442/expanded \
  outputs/results/prjdb36442/sensitivity \
  outputs/results/loombar2017 \
  outputs/results/qinn2014 \
  outputs/figures \
  outputs/tables

python3 scripts/prjdb36442_module_reanalysis.py \
  --merged-pathabundance data/prjdb36442/merged_pathabundance.tsv.gz \
  --manifest data/prjdb36442/sample_manifest.tsv \
  --module-definitions analysis/module_definitions.tsv \
  --membership conservative \
  --n-random 5000 \
  --out-dir outputs/results/prjdb36442/conservative

python3 scripts/prjdb36442_module_reanalysis.py \
  --merged-pathabundance data/prjdb36442/merged_pathabundance.tsv.gz \
  --manifest data/prjdb36442/sample_manifest.tsv \
  --module-definitions analysis/module_definitions.tsv \
  --membership expanded \
  --n-random 5000 \
  --out-dir outputs/results/prjdb36442/expanded

python3 scripts/prjdb36442_pathway_sensitivity.py \
  --merged-pathabundance data/prjdb36442/merged_pathabundance.tsv.gz \
  --manifest data/prjdb36442/sample_manifest.tsv \
  --module-definitions analysis/module_definitions.tsv \
  --out-dir outputs/results/prjdb36442/sensitivity

python3 scripts/score_processed_binary_modules.py \
  --cohort LoombaR_2017 \
  --pathway data/loombar2017/pathway_abundance_unstratified.tsv \
  --metadata data/loombar2017/metadata.tsv \
  --modules analysis/module_definitions.tsv \
  --out-dir outputs/results/loombar2017 \
  --reference-group F0_F2 \
  --test-group F3_F4 \
  --contrast F3_F4_minus_F0_F2

python3 scripts/score_processed_binary_modules.py \
  --cohort QinN_2014 \
  --pathway data/qinn2014/pathway_abundance_unstratified.tsv \
  --metadata data/qinn2014/metadata.tsv \
  --modules analysis/module_definitions.tsv \
  --out-dir outputs/results/qinn2014 \
  --reference-group healthy \
  --test-group cirrhosis \
  --contrast cirrhosis_minus_healthy

python3 scripts/build_figures_tables.py

printf 'Rebuilt analysis outputs in %s/outputs\n' "$ROOT"
