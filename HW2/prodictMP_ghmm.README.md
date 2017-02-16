prodictMP_ghmm.py is to apply Generalized Hidden Markov Model (GHMM) on the prediction of transmembrane domains of a protein.

There are two files in the folder which can be used for running the program. They are prodictMP_ghmm.py, TMseq.ffa

Input: python prodictMP_ghmm.py -i [INPUT_FILE] -o [OUTPUT_FILE]

Output: A output file for the probability of building GHMM matrix, maximum score and the predicted result
Sample usage: python prodictMP_ghmm.py -i TMseq.ffa -o output_ghmm.txt  
