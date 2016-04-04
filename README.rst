enrichme
======================================================================

Testing for enrichment of genome-wide scores or p-values
in gene-categories.
enrichme.py implements a method that naturally corrects for
gene length, linkage disequilibrium and gene clustering by comaring the data
to a null-distribution obtained by randomly rotating the scores
against their genomic positions.


Given you have results from a genome wide association study (GWAS),
a genome wide seletion scan or any kind of analysis that attributes
scores, likelihoods or p-values to positions across the genome,
enrichme.py test whether high scores are enriched in specific gene categories.


Installation
======================================================================

Use `pip <http://pip-installer.org>`_ or easy_install::

    pip install enrichme==0.1.0


Testing
======================================================================

You can run unit tests using the command::

    python setup.py test

Enrichment test implementations
======================================================================
The program currently implements three methods:

1. Candidate  (enrichme.py -M Candidate --help)\
    Comparing a candidate gene list to a background gene list.
    This is a standard function that is done by many enrichment
    analysis tools. No correction for gene length or LD.

2. TopScores (enrichme.py -M TopScores --help)\
    Check whether top ranking scores are within or close to genes
    enriched in specific gene-categories.
    This method naturally corrects for gene-length and LD.
    This method is useful if one expects that only scores above
    a threshold contain biologically relevant information.
    No information is used from the ranking of scores that pass
    the threshold.

3. Summary (enrichme.py -M Summary --help)\
    Similar to TopScores but instead of defining a threshold
    on the scores, a summary of scores is calculated for each gene
    category. As an example, for each gene category, the program could
    calculate the mean across the category of the maximum score across
    each gene in the category.
    This method is useful if one thinks that scores contain information
    down to low value or if the relative value of scores is important
    beyond defining a simple threshold.

How to use
======================================================================

::

    ./enrichme.py --help

Example scripts with further explanations and example data can be found in ./examples/

Test for enrichment of GWAS scores above 3 (p-value<10^-3) in GO-categories::

    cd ./examples/
    ../enrichme.py  -M TopScores \
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


For each gene, calculate the max GWAS score within the gene. Then, test for enrichment of the mean of these gene-scores in  GO-categories::

    cd ./examples/ 
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



Input files
======================================================================

A. FEATURE to CATEGORY mapping (input argument --feature_to_category)\
    This file maps genetic features (usually genes) to feature categories
    such as gene-lists. This could be GO-terms or custom defined gene-lists.
    File can be tab-separated (.tsv) or comma-separates (.csv)

    .. code::
    
        $head examples/example_gene_annotation.csv
        gene_id,go_identifier
        AT1G01010,GO:0005634
        AT1G01010,GO:0006355
        AT1G01010,GO:0003677
        AT1G01010,GO:0007275
        AT1G01010,GO:0003700
        AT1G01010,GO:0043090
        AT1G01010,GO:0006888
        AT1G01020,GO:0016125
        AT1G01020,GO:0016020


B. FEATURES (input argument --features)\
    This file gives the position of features (e.g. genes)
    across the genome. Often this will be the gene
    annotation. The column pos gives the start of the feature.

    .. code::
    
        $head examples/example_gene_annotation.csv
        chrom,start,end,gene_id
        1,3631,5899,AT1G01010
        1,5928,8737,AT1G01020
        1,11649,13714,AT1G01030
        1,23146,31227,AT1G01040
        1,28500,28706,AT1G01046
        1,31170,33153,AT1G01050
        1,33379,37840,AT1G01060
        1,38752,40944,AT1G01070
        1,44677,44787,AT1G01073

C. Scores across the genome (input argument --rod)\
    This could be position of SNPs and a
    score or p-value associated with them.
    ROD stands for Reference Ordered Data.

    .. code::

        $head examples/example_GWAS_result.csv
        chrom,pos,score
        1,3102,0.09305379
        1,4648,0.30615359999999997
        1,4880,0.35306350000000003
        1,5975,0.9596856
        1,6063,0.23715001
        1,6449,0.019213928
        1,6514,0.43630862
        1,6603,0.23235813
        1,6768,0.58977395

D. [Optional] Mapping of categories to category descriptions (input argument --category_to_description)\
    This could be a csv with GO-category ids and descriptions.

    .. code::

        $head examples/example_go_to_name.csv 
        go_identifier,go_name
        GO:0000001,mitochondrion inheritance
        GO:0000002,mitochondrial genome maintenance
        GO:0000003,reproduction
        GO:0042254,ribosome biogenesis
        GO:0044183,protein binding involved in protein folding
        GO:0051082,unfolded protein binding
        GO:0000006,high-affinity zinc uptake transmembrane transporter activity
        GO:0000007,low-affinity zinc ion transmembrane transporter activity
        GO:0003756,protein disulfide isomerase activity

Output
======================================================================

The different modes provide different output files. The main output file is common for all modes, called <name>.pvals.tsv. It is a ranked table with most significantly enriched categories on top::

    go_identifier   out_of  rank    score_summary   p_value benjamini_hochberg      go_name
    GO:0000165      2000    1943    0.8731354255802085      0.02898550724637683     27.014492753623205      MAPK cascade
    GO:0000041      2000    1825    0.8348620634942308      0.08795602198900554     27.32500416458439       transition metal ion transport
    GO:0000160      2000    1800    0.9736749697560976      0.1004497751124438      23.404797601199405      phosphorelay signal transduction system
    GO:0000164      2000    1698    1.0469719100000001      0.15142428785607198     28.225487256371814      protein phosphatase type 1 complex
    GO:0000096      2000    1692    0.8680123230000001      0.15442278860569714     23.987006496751622      sulfur amino acid metabolic process
    GO:0000145      2000    1685    0.9976431777777778      0.15792103948025982     21.02605839937174       exocyst
    GO:0000159      2000    1562    0.9504652303750003      0.21939030484757627     25.558970514742636      protein phosphatase type 2A complex
    GO:0000156      2000    1558    0.9427544812820514      0.22138930534732637     22.92609250930091       phosphorelay response regulator activity

Parallel support
======================================================================

There are two ways to run this program in parallel. Per default, the program uses as many cores as available on the host machine. This can be controlled with the --ncpus option. Advanced users, who want to parallelise across multiple nodes of a compute cluster, can use the built in map/reduce framework to automatically combine results from multiple independent runs. See

::

    examples/run_permute_reduce_examples.sh
    
for an example.

Changelog
======================================================================

**enrichme** follows `semantic versioning <http://semver.org>`_.  The
first release with stable API will be 1.0.0 (soon).  Until then, you
are encouraged to specify explicitly the version in your dependency
tools, e.g.::

    pip install enrichme==0.1.0

- 0.1.0 Initial release. 
