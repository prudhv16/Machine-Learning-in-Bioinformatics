#!/usr/bin/python
# Importing libraries
import re,sys,getopt
import numpy as np
from string import maketrans

bio_file = open("gene_1000.fasta","r")
lines = bio_file.readlines()

gene_seq = []
count = 0
for line in lines:
    if count == 1001:
        break
    elif line[0] == '>':
        count += 1
        pass
    else:
        #print line
        gene_seq.append(line)

gene_seq= [x.strip() for x in gene_seq]
seq = "".join(gene_seq)
gen = {}

c = []
for x in range(0,len(seq),3):
    if seq[x:x+3] in gen:
        gen[seq[x:x+3]] += 1
    else :
        gen[seq[x:x+3]] = 1

sum = np.sum(gen.values())



codan_usage = {}
codan_usage['TTT'] = float(gen['TTT'])/sum * 1000
codan_usage['TTC'] = float(gen['TTC'])/sum * 1000

codan_usage['TTA'] = float(gen['TTA'])/sum * 1000
codan_usage['TTG'] = float(gen['TTG'])/sum * 1000
codan_usage['CTT'] = float(gen['CTT'])/sum * 1000
codan_usage['CTC'] = float(gen['CTC'])/sum * 1000
codan_usage['CTA'] = float(gen['CTA'])/sum * 1000
codan_usage['CTG'] = float(gen['CTG'])/sum * 1000

codan_usage['ATT'] = float(gen['ATT'])/sum * 1000
codan_usage['ATC'] = float(gen['ATC'])/sum * 1000
codan_usage['ATA'] = float(gen['ATA'])/sum * 1000

codan_usage['ATG'] = float(gen['ATG'])/sum * 1000

codan_usage['GTT'] = float(gen['GTT'])/sum * 1000
codan_usage['GTC'] = float(gen['GTC'])/sum * 1000
codan_usage['GTA'] = float(gen['GTA'])/sum * 1000
codan_usage['GTG'] = float(gen['GTG'])/sum * 1000

codan_usage['TCT'] = float(gen['TCT'])/sum * 1000
codan_usage['TCC'] = float(gen['TCC'])/sum * 1000
codan_usage['TCA'] = float(gen['TCA'])/sum * 1000
codan_usage['TCG'] = float(gen['TCG'])/sum * 1000
codan_usage['AGT'] = float(gen['AGT'])/sum * 1000
codan_usage['AGC'] = float(gen['AGC'])/sum * 1000


codan_usage['CCT'] = float(gen['CCT'])/sum * 1000
codan_usage['CCC'] = float(gen['CCC'])/sum * 1000
codan_usage['CCA'] = float(gen['CCA'])/sum * 1000
codan_usage['CCG'] = float(gen['CCG'])/sum * 1000

codan_usage['ACT'] = float(gen['ACT'])/sum * 1000
codan_usage['ACC'] = float(gen['ACC'])/sum * 1000
codan_usage['ACA'] = float(gen['ACA'])/sum * 1000
codan_usage['ACG'] = float(gen['ACG'])/sum * 1000

codan_usage['GCT'] = float(gen['GCT'])/sum * 1000
codan_usage['GCC'] = float(gen['GCC'])/sum * 1000
codan_usage['GCA'] = float(gen['GCA'])/sum * 1000
codan_usage['GCG'] = float(gen['GCG'])/sum * 1000

codan_usage['TAT'] = float(gen['TAT'])/sum * 1000 
codan_usage['TAC'] = float(gen['TAC'])/sum * 1000

codan_usage['TAA'] = float(gen['TAA'])/sum * 1000 
codan_usage['TAG'] = float(gen['TAG'])/sum * 1000
codan_usage['TGA'] = float(gen['TGA'])/sum * 1000

codan_usage['CAT'] = float(gen['CAT'])/sum * 1000 
codan_usage['CAC'] = float(gen['CAC'])/sum * 1000

codan_usage['CAA'] = float(gen['CAA'])/sum * 1000 
codan_usage['CAG'] = float(gen['CAG'])/sum * 1000

codan_usage['AAT'] = float(gen['AAT'])/sum * 1000
codan_usage['AAC'] = float(gen['AAC'])/sum * 1000

codan_usage['AAA'] = float(gen['AAA'])/sum * 1000 
codan_usage['AAG'] = float(gen['AAG'])/sum * 1000

codan_usage['GAT'] = float(gen['GAT'])/sum * 1000 
codan_usage['GAC'] = float(gen['GAC'])/sum * 1000

codan_usage['GAA'] = float(gen['GAA'])/sum * 1000
codan_usage['GAG'] = float(gen['GAG'])/sum * 1000

codan_usage['TGT'] = float(gen['TGT'])/sum * 1000 
codan_usage['TGC'] = float(gen['TGC'])/sum * 1000

codan_usage['TGG'] = float(gen['TGG'])/sum * 1000

codan_usage['CGT'] = float(gen['CGT'])/sum * 1000
codan_usage['CGC'] = float(gen['CGC'])/sum * 1000
codan_usage['CGA'] = float(gen['CGA'])/sum * 1000
codan_usage['CGG'] = float(gen['CGG'])/sum * 1000
codan_usage['AGA'] = float(gen['AGA'])/sum * 1000
codan_usage['AGG'] = float(gen['AGG'])/sum * 1000

codan_usage['GGT'] = float(gen['GGT'])/sum * 1000
codan_usage['GGC'] = float(gen['GGC'])/sum * 1000
codan_usage['GGA'] = float(gen['GGA'])/sum * 1000
codan_usage['GGG'] = float(gen['GGG'])/sum * 1000

print "The codon usage library for ecoli is as below with each codon's frequency per 1000 base: "
print codan_usage

pattern = re.compile(r'(?=(ATG(?:...)*?(?=TAG|TGA|TAA)(?:...)))')

def readFasta(filename):
#Define a function to read every line in the file.
       fh = open(filename, 'r')
       line = fh.readline()
       name = ''
       sequence = ''
       while line:
	     line = line.rstrip('\n')
             if '>' in line:
		   name = line
             else:
		   sequence = sequence + line
             line = fh.readline()
       return sequence
       
def findORF(dna):
   return pattern.findall(dna), pattern.findall(revcomp(dna))

def revcomp(dna):
    return dna[::-1].translate(maketrans("ATGC", "TACG"))

def calLikelihood(orf_list):
    codon_usage_table = codan_usage
    p0_list = []
    p1_list = []
    likelihood = []
    for orf in orf_list:
        p1 = 2.7736331578153628 
        p0 = 1.5625
        for n in range(3, len(orf), 3):
	        if codon_usage_table.has_key(orf[n:n+3]):
		        p1 *=  codon_usage_table[orf[n:n+3]] / 10
		        p0 *=  1.5625		 
        p0_list.append(p0)

        p1_list.append(p1)

        tmp = p0/p1
        likelihood.append(np.log10(tmp))

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
		print "USAGE: python ECgnfinder.py -i [INPUT_FILE]"
		exit(-1)

  if sys.argv[1] == "-i":
		if sys.argv[2].endswith('.fasta'):
			input_file = sys.argv[2]
		else:
			print "The input file should be in fasta format"
			exit(-1)
  else:
		print "USAGE: python ECgnfinder.py -i [INPUT_FILE]"
		exit(-1)
		 
  dataset = readFasta(input_file)
  p_orf, p_orf2 = findORF(dataset)
  possible_orf = []
  for i in range(len(p_orf)):
      if len(p_orf[i]) > 110:
          possible_orf.append(p_orf[i])
          
  for i in range(len(p_orf2)):
      if len(p_orf2[i]) > 110:
          possible_orf.append(p_orf2[i])
       
  likelihood = calLikelihood(possible_orf)
  print "\n"
  print "Please find two output files with predicted genes as genes.fasta and protein.fasta in the folder."
  print "\n"
  filename1 = 'output_gene.fasta'
  f1 = open(filename1, 'w')
  filename2 = 'output_protein.fasta'
  f2 = open(filename2, 'w')

  
  for i in range(len(likelihood)):
      print "ORF %d: %s" % (i+1, likelihood[i])
      if  likelihood[i] > 0:
          f1.write(">ORF %d: " % (i+1))
          f1.write("\n")
          f1.write(possible_orf[i])
          f1.write("\n")
          protein = translation(possible_orf[i])
          f2.write(">Protein %d: " % (i+1))
          f2.write("\n")
          f2.write(protein)
          f2.write("\n")
      
