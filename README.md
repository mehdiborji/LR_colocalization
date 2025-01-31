# Inferring Ligand-Receptor Interactions from Spatial Transcriptomics Data

The script takes as input a list of path to `AnnData` objects and a list putative LR pairs from other databases or methods such `CellPhoneDB` or `LIANA` and produces statistical metrics for each pair for each adata object.

The script that run the loop of permutations tests can be be launched with the following command:

```
usage: run_LR_par.py [-h] [-c CORES] [-l LRS_FILE] [-p PUCKS_FILE] [-r RESULTS_DIR] [-n N_PAIRS] [-m M_PUCKS]
```

For example
```python run_LR_par.py -c 50 -l ../liana_FL.csv -p ../FL_pucks_with_path.csv -r ../fl_results -m 20 -n 100```

if either of `N_PAIRS` and `M_PUCKS` are not provided the script will go through the entirety of pairs or pucks.

`RESULTS_DIR` is organized by two categories `csv` files for each Ligand-Receptor pair:

files with prefix `discard` contain pucks that did not reach criteria for number of beads positive for ligand gene and/or receptor gene counts. In this case the csv contains names of pucks and number of positive beads only

files with prefix `results` contains the main reports for further analysis. In this case the csv contains more columns which include results of permutation test: namely p-value, number of observed colocalization events, as well as median , mean and standard deviation of colocalization events for permutated background, moreover the results are reports across different neighborhood scales currently hardcoded to be 15, 30, 100 and 300 pixel radii of positive spots.

