# Metagenomic-analysis
#Scripts that can be used to estimate the relative abundance of the recovered MAGs:

`1. MAGs-coverage-based-relabu.py`

This script will calculate the average coverage of different MAGs, then use the coverage to calculate the relative abundance. This script can be used if your recovered MAGs have represented the major proportion of the microbial community, e.g., recruited more than 70% of all the total reads. I do not recommend you to use this script to get the relative abundance if your MAGs only recruited small part of the sequenced reads.

Example:

	$ python MAGs-coverage-based-relabu.py -f ./MAGs -m s1-mapping-all-bins.txt -c 3 -l 8

Get help:

	$ python MAGs-coverage-based-relabu.py -h
	
This script can be used to estimate relative abundance of recovered MAGs.

-h,  : Print help

Required options:

-f   : folder containing all MAGs

-m   : Mapping matrix exported from CLC or other mapping tools, which summarize the ID of contig/scaffold, mapped reads of given   contig/scaffold, and length of this contig/scaffold. Example mapping matrix can be found in the folder of 'example-data'.

-c   : which column (number) contains the mapped reads number.

-l   : which column (number) contains the contigs/scaffolds length.

