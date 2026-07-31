# Riboseek-analysis

## Rfam_benchmark
This directory contains the queries and the targets used for the Rfam benchmark.

* query.fasta: contains the query sequences used for the benchmark.
* target.fasta: contains the target sequences used for the benchmark.
* target.tsv: contains if the target sequences are reverse complemented or not.
* rfam_lookup_target_clan.tsv: contains the lookup table for mapping target sequences to their respective families and clans in Rfam.
* query_difficult.list: contains the list of queries more than 20 members used for the benchmark.

To compute the ROC1-AUC, you can use the following command:
```
# First run nhmmer / blastn / or Riboseek to get the results in a tabular format (e.g., results.m8)

# If you have the results in BLAST tabular format (results.m8):
bash scripts/clean_m8.awk results.m8 > results_cleaned.m8
# If you have the results in nhmmer tabular format (results.tblout):
bash scripts/clean_hmmerout.awk results.tblout > results_cleaned.m8

# Then run the following command to compute the ROC1-AUC:
bash scripts/compute_auc.sh results_cleaned.m8 Rfam_benchmark/target.tsv Rfam_benchmark/rfam_lookup_target_clan.tsv roc1.count
awk 'NR==FNR {a[$1];next;} ($1 in a) {print}' Rfam_benchmark/query_difficult.list roc1.count > roc1_difficult.count

# Final ROC1-AUC
echo "ROC1-AUC  number_of_families"
awk '{sum+=$4} END {print sum/NR, NR}' roc1_difficult.count
```

If you want to compute precision-recall curves, you can use the following command:
```
# Use --score bits to get the precision-recall curve of vanilla Smith-Waterman
python scripts/bench_macro_pr.py --lookup Rfam_benchmark/rfam_lookup_target.tsv \
--m8 results_cleaned.m8 --mode evalue --score evalue --out pr_aln
```

## Dfam_benchmark
This directory contains the queries and the targets used for the Dfam benchmark.
* queries: contains the query sequences used for the benchmark.
* targets: contains the target sequences used for the benchmark.
* dfam_lookup_target.tsv: contains the lookup table for mapping target sequences to their respective families in Dfam.
* target_dfam.tsv: contains if the target sequences are reverse complemented or not.

With the queries and targets, run search tools (nhmmer / blastn / Riboseek) to get the results in a tabular format (e.g., results.m8), \
then run the same commands as above to compute the ROC1-AUC and precision-recall curves.

## Computing nhmmer MSA
Use `scripts/rna_to_a3m.py` to get the multiple sequence alignment (MSA) from given query sequences.\
You can run the following command:
```
python scripts/rna_to_a3m.py --input query.fasta --output_dir output_dir \
--db1 /path/to/RNAcentral.fasta --db2 /path/to/NT.fasta --max_sequences 30000
```
You need HMMER package to run the script.

## E-value
This directory contains the scripts to compute mu / lambda values for the E-value calibration, post-hoc calibration, and the script to plot the theoretical E-value with the empirical E-value.
* samplemulambda.cpp: compute mu / lambda values given the raw scores of the alignments.
* fit_x_to_y.R: fit the theoretical E-value to the empirical E-value.
* prep_plot_evalues.py: given the theoretical E-values, compute the empirical E-values and plot the theoretical E-values with the empirical E-values.