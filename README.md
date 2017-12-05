# PhyloRank

**[This project is in active development. Documentation is fairly light. You are welcomed to use this software, but please expect it to change in non-trivial ways.]**

[![version status](https://img.shields.io/pypi/v/phylorank.svg)](https://pypi.python.org/pypi/phylorank)

PhyloRank provide functionality for calculating the relative evolutionary divergence (RED) of taxa in a tree and for finding the best placement of taxonomic labels in a tree . Other functionality is in development and is currently unsupported.

## Install

The simplest way to install this package is through pip:
> sudo pip install phylorank

This package makes use of the [numpy](http://www.numpy.org/), [matplotlib](https://matplotlib.org), [jinja2](http://jinja.pocoo.org/docs/2.10/), [mpld3](http://mpld3.github.io/), and [dendropy](https://www.dendropy.org/) Python libraries.

## Calculating RED

PhyloRank can calculate the relative evolutionary divergence (RED) of taxa in a tree in order to identify taxa at a given taxonomic rank that have conspicuous placements in the tree. This information can be used to refine the placement of taxa in the tree with the aim of taxa at the same rank having similar RED. RED values can be calculated using:
```
>phylorank outliers <input_tree> <taxonomy_file> <output_dir>
```

where <input_tree> is a tree in Newick format that is decorated with taxonomic information, the <taxonomy_file> indicates the taxonomic assignment of each genome in the tree (see File Formats below), and <output_dir> is the location to store all results. This command assumes the taxonomically-decorated tree and taxonomy file are congruent. The taxonomy file is required in order to establish taxonomic affiliations that could not be provided in the tree (e.g., a species with a single representative and consequently no internal node in the tree). The output directory contains a table indicating the RED value calculated for all taxa along with a plot indicating the RED value of taxa at each rank.

## Decorating a Tree

PhyloRank can decorate a tree with taxonomic information. This is accomplished by determining the placement for each taxa that results in the highest F-measure (see [McDonald et al., 2012](https://www.ncbi.nlm.nih.gov/pubmed/22134646)). If the provided taxonomy is incongruent with the topology of the tree, only the most cohesive lineage of a polyphyletic group will be labelled and it is possible that taxa at the same rank may be nested within each other (e.g., a lineage may contain two genus labels). The resulting decorated tree can be manually curated to obtain a taxonomy that is consistent with the tree topology. Alternatively, the resulting F-measure scores can be examined to assess the congruency between the taxonomy and tree topology. A tree can be decorated using:
```
>phylorank decorate <input_tree> <taxonomy_file> <output_tree> --skip_rd_refine
```

where <input_tree> is a tree in Newick format, <taxonom_file> indicates the taxonomic assignment of each genome in the tree (see File Formats below), and <output_tree> is the desired name of the decorated tree. The --skip_rd_refine flag indicates that the placement of taxa in the tree should not be adjusted to account for their RED. Adjusting for RED is only recommended once the taxonomy and tree topology have been established as being largely congruent. The file <output_tree>-table indicates the F-measure, precision, and recall for each taxa in the tree and the file <output_tree>-taxonomy gives the assigned taxonomy of each genome in the tree. 

## File formats

The taxonomy file is a simple tab-separated values file with two columns indicating the genome ID and Greengenes-style taxonomy string, e.g.:
```
>GCF_001687105.1    d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Yangia;s__
```

## Cite

If you find this package useful, please cite this git repository (https://github.com/dparks1134/PhyloRank)

## Copyright

Copyright © 2015 Donovan Parks. See LICENSE for further details.
