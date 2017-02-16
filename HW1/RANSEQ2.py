import sys, random
import numpy as np

#define the script command line
if __name__ == "__main__":
  if len(sys.argv) != 7:
    print "USAGE: python RANSEQ2.py -i [INPUT_FILE] -n [Number_Of_Sequence] -o [OUTPUT_FILE]"
    exit(-1)

  if sys.argv[1] == "-i":
    if sys.argv[2].endswith('.fasta'):
      input_file = sys.argv[2]
    else:
      print "The input file should be in fasta format"
      exit(-1)
  else:
    print "USAGE: python RANSEQ2.py -i [INPUT_FILE] -n [Number_Of_Sequence] -o [OUTPUT_FILE]"
    exit(-1)

  if sys.argv[3] == "-n":
    num = int(sys.argv[4])
  else:
    print "USAGE: python RANSEQ2.py -i [INPUT_FILE] -n [Number_Of_Sequence] -o [OUTPUT_FILE]"
    exit(-1)

  if sys.argv[5] == "-o":
    output_file = sys.argv[6]
  else:
    print "USAGE: python RANSEQ2.py -i [INPUT_FILE] -n [Number_Of_Sequence] -o [OUTPUT_FILE]"
    exit(-1)



#read a fasta file and return the sequence as a string
def FASTA(filename):
  f = file(filename)
  order = []
  sequence = ""
    
  for line in f:
    if line.startswith('>'):
      name = line[1:].rstrip('\n')
      name = name.replace('_', ' ')
      order.append(name)
    else:
      sequence += line.rstrip('\n').rstrip('*')
            
  return order, sequence


strain, seq = FASTA(input_file)
sequence = list(seq)
length = len(sequence)
f = open(output_file, 'w')


#count nucleotide frequency
def Nucleotide_Freq(sequen):
    adenine = 0
    thymine = 0
    cytosine = 0
    guanine = 0

    for i in sequen:
        if i == 'A':
            adenine += 1
        elif i == 'T':
            thymine += 1
        elif i == 'C':
            cytosine += 1
        else:
            guanine += 1  

    freq_A = float(adenine) / length
    freq_T = float(thymine) / length
    freq_C = float(cytosine) / length
    freq_G = float(guanine) / length

    freq = {}
    freq['A'] = freq_A
    freq['T'] = freq_T
    freq['C'] = freq_C
    freq['G'] = freq_G
    freq = str(freq)

    mean = (freq_A + freq_T + freq_C + freq_G)/4
    SD = np.std([freq_A, freq_T, freq_C, freq_G])

    f.write( "The frequency of each nucleotide is: ")
    f.write( freq )
    f.write( "\n" )
    f.write( "\n" )

    return freq_A, freq_T, freq_C, freq_G


f.write("The input strain is: %s \n" % seq) 
Ade, Thy, Cyt, Gua = Nucleotide_Freq(sequence)

#generate required number of output sequences with fixed frequency

SD_Ade, SD_Thy, SD_Cyt, SD_Gua = [], [], [], []

for i in range(num):
    #permutate the template string
    ran_seq = random.sample(sequence, len(sequence))
    ran_seq_string = ''.join(ran_seq)
    f.write( "The %dth permuted sequence is: " % (i+1) )
    f.write( ran_seq_string )
    f.write( "\n" )
    A, T, C, G = Nucleotide_Freq(ran_seq)
    SD_Ade.append(A)
    SD_Thy.append(T)
    SD_Cyt.append(C)
    SD_Gua.append(G)

mean_A = np.mean(SD_Ade)
mean_T = np.mean(SD_Thy)
mean_C = np.mean(SD_Cyt)
mean_G = np.mean(SD_Gua)

SD_A = np.std(SD_Ade)
SD_T = np.std(SD_Thy)
SD_C = np.std(SD_Cyt)
SD_G = np.std(SD_Gua)

print SD_G

f.write("The mean frequency of Adenine is %s, Thymine is %s, Cytosine is %s, Guannine is %s." % (mean_A, mean_T, mean_C, mean_G))
f.write("\n")
f.write("The SD of Adenine is %s, Thymine is %s, Cytosine is %s, Guannine is %s." % (SD_A, SD_T, SD_C, SD_G))

print "Please open the output file %s for results." % output_file
    