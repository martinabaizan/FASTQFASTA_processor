# FASTQFASTA_processor

A python tool that, given a FASTA or FASTQ file, preprocesses all the sequences contained in it. The tool performs the following operations on the sequences: reverse-complement, trim, and adaptor-removal. 
Depending on the command line options used, the tool will perform a different preprocessing operation.

Usage:
> python FASTQFASTA_processor.py --input [FILE]--output [FILE] --operation [OPERATION]
> [OPERATION] in {rc,trim,adaptor-removal}
