# Visualize phylogenetic tree using Auspice

## Install Auspice

```
curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml
conda env create -f nextstrain.yml
conda activate nextstrain
npm install --global auspice
```

## Run script to generate interactive phylogenetic tree

```
augur refine --tree data/final_tree_incompletenames.nwk --metadata data/metadata.tsv --output-tree results/final_tree_incompletenames_tt.nwk --output-node-data results/branch_lengths.json --keep-root

augur export v2 --tree results/final_tree_incompletenames_tt.nwk --metadata data/metadata.tsv --node-data results/branch_lengths.json --colors config/colors.tsv --lat-longs config/lat_longs.tsv --auspice-config config/auspice_config.json --output auspice/seubayanus_auspice.json
```

## Visualize tree

```
auspice view --datasetDir auspice
```

Navigate to http://localhost:4000/
