# Applying-data-specific-substitution-models-and-mitigating-the-effects-of-among-lineage-heterogeneity

## Chapter III - Measuring the efficacy of matched-pairs tests of symmetry and p-value adjustments in identifying tree-heterogeneous sequences in protein data sets.
This directory contains the novel scripts used (0_bin), the simulated data (1_simulatedData), and the empirical data and their resulting trees (2_empiricalData). Each file within the "2_empiricalData" includes the alignment, the data-specific model, and the resulting tree.

## chapter IV - Identifying the closest living algal relatives of land plants using data-specific amino acid substitution .models to analyse each of the three plant genomes.
This directory contains the novel scripts used (0_bin), the constructed data sets, single and combined (1_datasets), and the inferred trees and data-specific models (2_trees_and_models).

## Repository structure


├── chapter_III
|   ├── 0_bin
│   ├── 1_simulatedData
│   │   ├── 1_compoHetero_ratesHetero_jtt_stmtREV
│   │   ├── 2_compoHomo_ratesHetero_jtt_stmtREV
│   │   └── 3_compoHetero_ratesHomo_jtt_cpREV
│   └── 2_empiricalData
│       ├── 1_liu
│       ├── 2_strassert
│       └── 3_whelan
└── chapter_IV
    ├── 0_bin
    ├── 1_datasets
    │   ├── chloroplast_data
    │   ├── mitochondrial_data
    │   └── nuclear_data
    └── 2_trees_and_models
        ├── chloroplastTrees
        ├── mitochondrialTrees
        └── nuclearTrees
