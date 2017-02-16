ECgnfinder.py is aimed at finding potential e.coli genes based on e.coli codon usage.
Three files should be in a same folder to make sure the running of the program. They are ECgnfinder.py, gene_1000.fasta and Given_DNA.fasta.
Input: python ECgnfinder.py -i [INPUT_FILE]
Output: print out a codon usage table, print out every predicted ORF's likelihood. Two fasta files, one with predicted genes one with translated proteins.

Sample usage: python ECgnfinder.py -i Given_DNA.fasta
