import sys
import numpy as np

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
print gen

codan_usage = {}
codan_usage['TTT'] = float(gen['TTT'])/(gen['TTT']+gen['TTC'])
codan_usage['TTC'] = float(gen['TTC'])/(gen['TTC']+gen['TTT'])

codan_usage['TTA'] = float(gen['TTA'])/(gen['TTA']+gen['TTG']+gen['CTT']+gen['CTC']+gen['CTA']+gen['CTG'])
codan_usage['TTG'] = float(gen['TTG'])/(gen['TTA']+gen['TTG']+gen['CTT']+gen['CTC']+gen['CTA']+gen['CTG'])
codan_usage['CTT'] = float(gen['CTT'])/(gen['TTA']+gen['TTG']+gen['CTT']+gen['CTC']+gen['CTA']+gen['CTG'])
codan_usage['CTC'] = float(gen['CTC'])/(gen['TTA']+gen['TTG']+gen['CTT']+gen['CTC']+gen['CTA']+gen['CTG'])
codan_usage['CTA'] = float(gen['CTA'])/(gen['TTA']+gen['TTG']+gen['CTT']+gen['CTC']+gen['CTA']+gen['CTG'])
codan_usage['CTG'] = float(gen['CTG'])/(gen['TTA']+gen['TTG']+gen['CTT']+gen['CTC']+gen['CTA']+gen['CTG'])

codan_usage['ATT'] = float(gen['ATT'])/(gen['ATT']+gen['ATC']+gen['ATA'])
codan_usage['ATC'] = float(gen['ATC'])/(gen['ATT']+gen['ATC']+gen['ATA'])
codan_usage['ATA'] = float(gen['ATA'])/(gen['ATT']+gen['ATC']+gen['ATA'])

codan_usage['ATG'] = float(gen['ATG'])/(gen['ATG'])

codan_usage['GTT'] = float(gen['GTT'])/(gen['GTT']+gen['GTC']+gen['GTA']+gen['GTG'])
codan_usage['GTC'] = float(gen['GTC'])/(gen['GTT']+gen['GTC']+gen['GTA']+gen['GTG'])
codan_usage['GTA'] = float(gen['GTA'])/(gen['GTT']+gen['GTC']+gen['GTA']+gen['GTG'])
codan_usage['GTG'] = float(gen['GTG'])/(gen['GTT']+gen['GTC']+gen['GTA']+gen['GTG'])

codan_usage['TCT'] = float(gen['TCT'])/(gen['TCT']+gen['TCC']+gen['TCA']+gen['TCA']+gen['AGT']+gen['AGC'])
codan_usage['TCC'] = float(gen['TCC'])/(gen['TCT']+gen['TCC']+gen['TCA']+gen['TCA']+gen['AGT']+gen['AGC'])
codan_usage['TCA'] = float(gen['TCA'])/(gen['TCT']+gen['TCC']+gen['TCA']+gen['TCA']+gen['AGT']+gen['AGC'])
codan_usage['TCG'] = float(gen['TCG'])/(gen['TCT']+gen['TCC']+gen['TCA']+gen['TCA']+gen['AGT']+gen['AGC'])
codan_usage['AGT'] = float(gen['AGT'])/(gen['TCT']+gen['TCC']+gen['TCA']+gen['TCA']+gen['AGT']+gen['AGC'])
codan_usage['AGC'] = float(gen['ACG'])/(gen['TCT']+gen['TCC']+gen['TCA']+gen['TCA']+gen['AGT']+gen['AGC'])


codan_usage['CCT'] = float(gen['CCT'])/(gen['CCT']+gen['CCC']+gen['CCA']+gen['CCG'])
codan_usage['CCC'] = float(gen['CCC'])/(gen['CCT']+gen['CCC']+gen['CCA']+gen['CCG'])
codan_usage['CCA'] = float(gen['CCA'])/(gen['CCT']+gen['CCC']+gen['CCA']+gen['CCG'])
codan_usage['CCG'] = float(gen['CCG'])/(gen['CCT']+gen['CCC']+gen['CCA']+gen['CCG'])

codan_usage['ACT'] = float(gen['ACT'])/(gen['ACT']+gen['ACC']+gen['ACA']+gen['ACG'])
codan_usage['ACC'] = float(gen['ACC'])/(gen['ACT']+gen['ACC']+gen['ACA']+gen['ACG'])
codan_usage['ACA'] = float(gen['ACA'])/(gen['ACT']+gen['ACC']+gen['ACA']+gen['ACG'])
codan_usage['ACG'] = float(gen['ACG'])/(gen['ACT']+gen['ACC']+gen['ACA']+gen['ACG'])

codan_usage['GCT'] = float(gen['GCT'])/(gen['GCT']+gen['GCC']+gen['GCA']+gen['GCG'])
codan_usage['GCC'] = float(gen['GCC'])/(gen['GCT']+gen['GCC']+gen['GCA']+gen['GCG'])
codan_usage['GCA'] = float(gen['GCA'])/(gen['GCT']+gen['GCC']+gen['GCA']+gen['GCG'])
codan_usage['GCG'] = float(gen['GCG'])/(gen['GCT']+gen['GCC']+gen['GCA']+gen['GCG'])

codan_usage['TAC'] = float(gen['TAT'])/(gen['TAT']+gen['TAC']) 
codan_usage['TAC'] = float(gen['TAC'])/(gen['TAT']+gen['TAC'])

codan_usage['TAA'] = float(gen['TAA'])/(gen['TAA']+gen['TAG']+gen['TGA']) 
codan_usage['TAG'] = float(gen['TAG'])/(gen['TAA']+gen['TAG']+gen['TGA'])
codan_usage['TGA'] = float(gen['TGA'])/(gen['TAA']+gen['TAG']+gen['TGA'])

codan_usage['CAT'] = float(gen['CAT'])/(gen['CAT']+gen['CAC']) 
codan_usage['CAC'] = float(gen['CAC'])/(gen['CAT']+gen['CAC'])

codan_usage['CAA'] = float(gen['CAA'])/(gen['CAA']+gen['CAG']) 
codan_usage['CAG'] = float(gen['CAG'])/(gen['CAA']+gen['CAG'])

codan_usage['AAT'] = float(gen['AAT'])/(gen['AAT']+gen['AAC']) 
codan_usage['AAC'] = float(gen['AAC'])/(gen['AAT']+gen['AAC'])

codan_usage['AAA'] = float(gen['AAA'])/(gen['AAA']+gen['AAG']) 
codan_usage['AAG'] = float(gen['AAG'])/(gen['AAA']+gen['AAG'])

codan_usage['GAT'] = float(gen['GAT'])/(gen['GAT']+gen['GAC']) 
codan_usage['GAC'] = float(gen['GAC'])/(gen['GAT']+gen['GAC'])

codan_usage['GAA'] = float(gen['GAA'])/(gen['GAA']+gen['GAG']) 
codan_usage['GAG'] = float(gen['GAG'])/(gen['GAA']+gen['GAG'])

codan_usage['TGT'] = float(gen['TGT'])/(gen['TGT']+gen['TGC']) 
codan_usage['TGC'] = float(gen['TGC'])/(gen['TGT']+gen['TGC'])

codan_usage['TGG'] = float(gen['TGG'])/(gen['TGG'])

codan_usage['CGT'] = float(gen['CGT'])/(gen['CGT']+gen['CGC']+gen['CGA']+gen['CGG']+gen['AGA']+gen['AGG'])
codan_usage['CGC'] = float(gen['CGC'])/(gen['CGT']+gen['CGC']+gen['CGA']+gen['CGG']+gen['AGA']+gen['AGG'])
codan_usage['CGA'] = float(gen['CGA'])/(gen['CGT']+gen['CGC']+gen['CGA']+gen['CGG']+gen['AGA']+gen['AGG'])
codan_usage['CGG'] = float(gen['CGG'])/(gen['CGT']+gen['CGC']+gen['CGA']+gen['CGG']+gen['AGA']+gen['AGG'])
codan_usage['AGA'] = float(gen['AGA'])/(gen['CGT']+gen['CGC']+gen['CGA']+gen['CGG']+gen['AGA']+gen['AGG'])
codan_usage['AGG'] = float(gen['AGG'])/(gen['CGT']+gen['CGC']+gen['CGA']+gen['CGG']+gen['AGA']+gen['AGG'])

codan_usage['GGT'] = float(gen['GGT'])/(gen['GGT']+gen['GGC']+gen['GGA']+gen['GGG'])
codan_usage['GGC'] = float(gen['GGC'])/(gen['GGT']+gen['GGC']+gen['GGA']+gen['GGG'])
codan_usage['GGA'] = float(gen['GGA'])/(gen['GGT']+gen['GGC']+gen['GGA']+gen['GGG'])
codan_usage['GGG'] = float(gen['GGG'])/(gen['GGT']+gen['GGC']+gen['GGA']+gen['GGG'])

print codan_usage

