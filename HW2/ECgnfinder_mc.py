#!/usr/bin/python
# Importing libraries
import re,sys,getopt
import numpy as np
from string import maketrans

bio_file = "test.fasta"

def FASTA(filename):
  try:
    f = file(filename)
  except IOError:                     
    print "The file, %s, does not exist" % filename
    return
  order = []
  sequences = {} 
  for line in f:
    if line.startswith('>'):
      name = line[1:].rstrip('\n')
      order.append(name)
      sequences[name] = ''
    else:
      sequences[name] += line.rstrip('\n').rstrip('*')
  return sequences
       
fasta = FASTA(bio_file)
mc = []
for v in fasta.values():
    mc.append(v)
#print mc

#1st Markov Chain
AA,AT,AC,AG,TA,TT,TC,TG,CA,CT,CC,CG,GA,GT,GC,GG = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

for element in mc:
    seq = list(element)
        
    for i in range(len(seq) - 1):
        if seq[i] == 'A' and seq[i+1] == 'A': 
            AA += 1
        elif seq[i] == 'A' and seq[i+1] == 'T':
            AT += 1
        elif seq[i] == 'A' and seq[i+1] == 'C':
            AC += 1
        elif seq[i] == 'A' and seq[i+1] == 'G':
            AG += 1
        elif seq[i] == 'T' and seq[i+1] == 'A':
            TA += 1
        elif seq[i] == 'T' and seq[i+1] == 'T':
            TT += 1
        elif seq[i] == 'T' and seq[i+1] == 'C':
            TC += 1
        elif seq[i] == 'T' and seq[i+1] == 'G':
            TG += 1
        elif seq[i] == 'C' and seq[i+1] == 'A':
            CA += 1
        elif seq[i] == 'C' and seq[i+1] == 'T':
            CT += 1
        elif seq[i] == 'C' and seq[i+1] == 'C':
            CC += 1
        elif seq[i] == 'C' and seq[i+1] == 'G':
            CG += 1
        elif seq[i] == 'G' and seq[i+1] == 'A':
            GA += 1
        elif seq[i] == 'G' and seq[i+1] == 'T':
            GT += 1
        elif seq[i] == 'G' and seq[i+1] == 'C':
            GC += 1
        else:
            GG += 1    
    
markov_chain = {}
markov_chain['AA'] = float(AA) / (AA+AT+AC+AG)
markov_chain['AT'] = float(AT) / (AA+AT+AC+AG)
markov_chain['AC'] = float(AC) / (AA+AT+AC+AG)
markov_chain['AG'] = float(AG) / (AA+AT+AC+AG)
markov_chain['TA'] = float(TA) / (TA+TT+TC+TG)
markov_chain['TT'] = float(TT) / (TA+TT+TC+TG)
markov_chain['TC'] = float(TC) / (TA+TT+TC+TG)
markov_chain['TG'] = float(TG) / (TA+TT+TC+TG)
markov_chain['CA'] = float(CA) / (CA+CT+CC+CG)
markov_chain['CT'] = float(CT) / (CA+CT+CC+CG)
markov_chain['CC'] = float(CC) / (CA+CT+CC+CG)
markov_chain['CG'] = float(CG) / (CA+CT+CC+CG)
markov_chain['GA'] = float(GA) / (GA+GT+GC+GG)
markov_chain['GT'] = float(GT) / (GA+GT+GC+GG)
markov_chain['GC'] = float(GC) / (GA+GT+GC+GG)
markov_chain['GG'] = float(GG) / (GA+GT+GC+GG)


print "The 1st Markov chain library for ecoli is as below: "
print markov_chain



def calLikelihood(seq_list):
    p0 = 1000
    p1 = 2.5*2.5*2.5
    for i in range(3, len(seq_list)):
        s = seq_list[i-1]+seq_list[i]
        p0 = p1 * markov_chain[s] * 5
        p1 *= 1.25	 
    likelihood = np.log10(p0/p1)
    #print p0
    #print p1
    #print likelihood
    return likelihood
    
def translation(sequence):
    codon_table = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 
    'TAT':'Y', 'TAC':'Y', 'TAA':'_', 'TAG':'_', 'TGT':'C', 'TGC':'C', 'TGA':'_', 'TGG':'W',
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCG':'A', 'GCA':'A',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    }
        
    protein_sequence = ''    
    for n in range(0, len(sequence), 3): #translate the sequence from position0
        if codon_table.has_key(sequence[n:n+3]) == True: #make sure it doesn't return errors if there are less than 3 bases left
            protein_sequence += codon_table[sequence[n:n+3]] #translate based on the dictionary above
            
    return protein_sequence
	
if __name__ == "__main__":
  if len(sys.argv) != 3:
		print "USAGE: python ECgnfinder_mc.py -i [INPUT_FILE]"
		exit(-1)

  if sys.argv[1] == "-i":
		if sys.argv[2].endswith('.fasta'):
			input_file = sys.argv[2]
		else:
			print "The input file should be in fasta format"
			exit(-1)
  else:
		print "USAGE: python ECgnfinder_mc.py -i [INPUT_FILE]"
		exit(-1)
		 
  dataset = FASTA(input_file)
  lh = {}
  for key in dataset.keys():
    seq = list(dataset[key])
    likelihood = calLikelihood(seq)   
    lh[key] = likelihood

  print "\n"
  print "Please find two output files with predicted genes as genes.fasta and protein.fasta in the folder."
  print "\n"
  filename1 = 'output_gene.fasta'
  f1 = open(filename1, 'w')
  filename2 = 'output_protein.fasta'
  f2 = open(filename2, 'w')

  count = 0
  for key in lh.keys():
      if  lh[key] > 0 :
          count += 1
          print "The %s sequence is predicted as an ecoli gene withh likelihood score %s" % (key, lh[key])
          print "Total of %d sequences are predicted as ecoli genes." % count
          f1.write(">")
          f1.write(key)
          f1.write("\n")
          f1.write(dataset[key])
          f1.write("\n")
          protein = translation(dataset[key])
          f2.write(">")
          f2.write(dataset[key])
          f2.write("\n")
          f2.write(protein)
          f2.write("\n")
      
