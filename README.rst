Enrichme
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

You can run unit tests using the command:

    python setup.py test

How to use
======================================================================

 $enrichme.py -R Permute --help

 The program currently implements three methods:
 (1) Candidate  (enrichme.py -R Permute -M Candidate --help)
     Comparing a candidate gene list to a background gene list.
     This is a standard function that is done by many enrichment
     analysis tools. No correction for gene length or LD.

 (2) TopScores (enrichme.py -R Permute -M TopScores --help)
     Check whether top ranking scores are within or close to genes
     enriched in specific gene-categories.
     This method naturally corrects for gene-length and LD.
     This method is useful if one expects that only scores above
     a threshold contain biologically relevant information.
     No information is used from the ranking of scores that pass
     the threshold.

 (3) Summary (enrichme.py -R Permute -M Summary --help)
     Similar to TopScores but instead of defining a threshold
     on the scores, a summary of scores is calculated for each gene
     category. As an example, for each gene category, the program could
     calculate the mean across the category of the maximum score across
     each gene in the category.
     This method is useful if one thinks that scores contain information
     down to low value or if the relative value of scores is important
     beyond defining a simple threshold.

INPUT files:

(A) FEATURE to CATEGORY mapping (input argument --feature_to_category)
    This file maps genetic features (usually genes) to feature categories
    such as gene-lists. This could be GO-terms or custom defined gene-lists.
    File can be tab-separated (.tsv) or comma-separates (.csv)

    EXAMPLE
    >head examples/example_gene_annotation.csv
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


(B) FEATURES (input argument --features)
    This file gives the position of features (e.g. genes)
    across the genome. Often this will be the gene
    annotation. The column pos gives the start of the feature.

    EXAMPLE
    >head examples/example_gene_annotation.csv
    chrom,pos,end,gene_id
    1,3631,5899,AT1G01010
    1,5928,8737,AT1G01020
    1,11649,13714,AT1G01030
    1,23146,31227,AT1G01040
    1,28500,28706,AT1G01046
    1,31170,33153,AT1G01050
    1,33379,37840,AT1G01060
    1,38752,40944,AT1G01070
    1,44677,44787,AT1G01073

(C) Scores across the genome (input argument --rod)
    This could be position of SNPs and a
    score or p-value associated with them.
    ROD stands for Reference Ordered Data.

    EXAMPLE
    >head examples/example_GWAS_result.csv
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


Changelog
======================================================================

**enrichme** follows `semantic versioning <http://semver.org>`_.  The
first release with stable API will be 1.0.0 (soon).  Until then, you
are encouraged to specify explicitly the version in your dependency
tools, e.g.::

    pip install enrichme==0.1.0

- 0.1.0 Initial release. 
