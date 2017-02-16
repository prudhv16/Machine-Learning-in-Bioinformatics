This program pruning.py is used to implement Felsenstein's pruning algorithm on aligned sequences.
Input: Two files, one is a phylogenetic tree in newick standard format; one is multiple alignment of nucleotide sequences without gaps.
Output: The alignment as a dictionary; 
        A phylo tree in newick standard format; 
        Three substitution matrix; 
        The best ancestor sequence with the best probability of each column;
        The total probability of getting the best ancestor sequence.

Sample Usage: python pruning.py Tree.nwk Alignment.txt
The best ancestor sequence as well as the probability will be printed out on the screen.