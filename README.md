# $${\color{lightgreen}SeqDepth checker}$$

Welcome to the SeqDepth_checker app, created to evaluate the depth of coverage in Next-Generation Sequencing experiments for every gene region listed in the input BED file. 

The following options are available:

	-b, --bedfile: A BED file with the gene regions the enrichment of which should be evaluated.
	-d, --data_dir: A path to a folder with BAM files to be tested.
	-o, --output_dir: A path to a folder where the resultant Excel file is saved (default: the data_dir).
	-t, --threads: Number of CPU threads to be used (default: 1).
	-s, --bam_suffix: A character string separating a sample name from the file extension (.bam).
	-S, --sd_times: The number of standard deviations to be subtracted from each mean sequencing read coverage depth value to get a minimal read coverage depth value (default: 0).
	-T, --threshold: If the minimal sequencing read coverage depth value for a gene region is lower than the given threshold, the enrichment of this region is considered as failed (default: 5).
	-h, --help: prints this help message.
	-v, --version: prints the version of this program.
