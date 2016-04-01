#!/usr/bin/env bash
#------------------------------------------------
# Shell script with examples to run enrichme.py #
#------------------------------------------------

# Example of mutliple runs on the same data set
# that are reduced to a single result.
# The idea is that runs can be executed in parallel
# locally or on a server.

# For more basic examples and input explanation
#  see run_minimal_examples.sh
#



echo PERMUTE REDUCE EXAMPLE 
echo " 4 parallel runs of Summary enrichment."



echo "A first dry run with 0 permutations is needed to create basic summary files."

name="permute_reduce_test"

../enrichme.py -R Permute -M Summary \
               --feature_to_category example_gene_to_category.csv \
               --feature_to_category_cols gene_id go_identifier  \
               --rod example_GWAS_result.csv \
               --rod_cols chrom pos score  \
               --features example_gene_annotation.csv \
               --feature_cols chrom start end gene_id \
               --name ${name}  \
               --n_permut 0 \
               --feature_summary max \
               --category_summary mean \


filenames=""

for i in `seq 1 5` #set to the number of cores
do
    echo "Starting permute job ${i}"
    filenames+="${name}_${i}.ranktable.tsv "
    ../enrichme.py -R Permute -M Summary \
                   --feature_to_category example_gene_to_category.csv \
                   --feature_to_category_cols gene_id go_identifier  \
                   --rod example_GWAS_result.csv \
                   --rod_cols chrom pos score  \
                   --features example_gene_annotation.csv \
                   --feature_cols chrom start end gene_id \
                   --name ${name}_${i}  \
                   --n_permut 10 \
                   --feature_summary max \
                   --category_summary mean \
                   --noinfo & 
done


# NOTE 1:
# On a compute cluster, you would modify the above by just submitting each
# run as a separate job.
# I.e., you would add the submit command before ../enrichme.py and remove
# the railing &

# NOTE 2:
# The total number of permutations you get is (number of runs) * n_permut

# NOTE 3:
# To suppress output of the parallel jobs, insert  &>/dev/null BEFORE the
# final &

echo "Wait for jobs to finish before starting reduce!"
wait
echo "Reduce..."

../enrichme.py -R Reduce \
               --feature_to_category example_gene_to_category.csv \
               --feature_to_category_cols gene_id go_identifier  \
               --name ${name}  \
               --permuts ${filenames} \
               --remove_input



