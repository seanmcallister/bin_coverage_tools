# Bin coverage tools

A collection of tools for calculating the "amount" of a bin. Specifically, the genomic coverage and RNA read count per bin. When you call the program use ```-h``` to see the help information.

Tools:

1. ```bin_coverage_RNAseq_recruitment.pl```
2. ```bin_coverage_coassembly.pl```
3. ```bin_coverage_individualassembly.pl```

Tool 1 calculates the total number of reads recruiting to each bin. The infile depends on the idxstats output from samtools, which is generated from your favorite recruitment software bam output (bowtie, BWA, etc).

Tools 2 & 3 calculate an average coverage per bin, for use in estimating the abundance of bins within a metagenome. The input is from the script ```jgi_summarize_bam_contig_depths```, which can be found bundled with [MetaBat](https://bitbucket.org/berkeleylab/metabat). The depth here is not a simple average; rather the coverage is weighted in favor of longer contigs.
