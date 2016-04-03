#!/usr/bin/env bash
#------------------------------------------------
# Shell script with examples to run enrichme.py #
#------------------------------------------------

# for further help see README.rst 
# or run ./enrichme.py --R Permute --help

#EXAMPLE DATA

# (1) Results from a genome wide association study (GWAS)
#     in Arabidopsis thaliana.
#     scores=-log10(pvalue) for SNP-positions across the genome.
#     File: example_GWAS_result.csv

# (2) Gene annotation for A. thaliana.
#     File: example_gene_annotation.csv

# (3) Gene to gene ontology (GO) category mapping.
#     File: example_gene_to_category.csv

#echo MINIMAL EXAMPLE 1
#echo "Simple single run of TopScores enrichment."
echo Test for enrichment in scores above 3 (p-value<10^-3)

../enrichme.py -R Permute -M TopScores \
               --feature_to_category example_gene_to_category.csv \
               --feature_to_category_cols gene_id go_identifier  \
               --rod example_GWAS_result.csv \
               --rod_cols chrom pos score  \
               --features example_gene_annotation.csv \
               --feature_cols chrom start end gene_id \
               --name minimal_test_TopScores  \
               --n_permut 10 \
               --top_type threshold \
               --top 3 \
               --descending \
               --max_dist 5000 \


#arguments explained:

# --<input_file>_cols ...
# For each input file, one also needs to specify the colum headers for
# the expected columns. This means the input file can have a different 
# column order or additional columns that are ignored.
# E.g., the flag --feature_cols specifies the header names of the chrom, start, end, feature_id
# columns in the feature (gene) annotation file.
#
# --name ... base name for all output files
# --n_permut ... number of permutations. # --top_type ... Kind of threshold to apply. Other options are quantile or count.
# --top ... says that genes close to scores>3 are tested
# --descending ... says that higher-score means more significant,
#                  for p-values, use descending instead
# --max_dist ... means that genes that fall within +-5000 basepairs
#                of a score are considered

# NOTE on the number of permutations:
# In practise, you need at least
# 100*(#gene categories tested) permutations to have a change
# of getting a result significant above multiple testing.
# E.g., if you test 3000 GO_categories, you need 30000 permutations.
# For this, you will want to run the method in parallel on multiple cores
# (e.g., --ncpus 16) or submit multiple runs to a compute cluster
# and reduce the results. (See file run_permute_reduce_examples.sh).



echo MINIMAL EXAMPLE 2
echo "Simple single run of Summary enrichment."
echo Test for enrichment of the mean across a category of the max scores across genes in this category.
echo Are there categories with high average gene scores?


../enrichme.py -M Summary \
               --feature_to_category example_gene_to_category_40_cats.csv \
               --feature_to_category_cols gene_id go_identifier  \
               --category_to_description example_go_to_name.csv \
               --category_to_description_cols go_identifier go_name \
               --rod example_GWAS_result.csv \
               --rod_cols chrom pos score  \
               --features example_gene_annotation.csv \
               --feature_cols chrom start end gene_id \
               --name minimal_test_Summary  \
               --feature_summary max \
               --category_summary mean \



