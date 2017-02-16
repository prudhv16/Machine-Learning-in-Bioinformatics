import sys
import numpy as np
#reading from file

def preprocessing(file_name):
	#file = open("TMseq.ffa",'r')
	file = open(file_name,'r')
	lines = file.readlines()
	for x in lines:
    		x = x.strip()

	protein = []
	pos = []
	for x in lines:
    		protein1, pos1 = x.split('\t')
    		protein.append(protein1)
    		pos1 = pos1.strip()
    		pos.append(pos1)
	return protein,pos

#Calculating initial probability

def initialprob(pos):
    initial_prob = {}
    initial = {}
    sum = 0
    for x in pos:
        if x[0] == '.':
            i = 0
            while x[i] == '.':
                i += 1
            if x[i] in initial_prob:
                initial_prob[x[i]] += 1
                sum += 1
            else:
                initial_prob[x[i]] = 1
                sum += 1
        elif x[0] in initial_prob:
        	initial_prob[x[0]] += 1
        	sum += 1
    	else:
    		initial_prob[x[0]] = 1
    		sum += 1
    for key in initial_prob:
        initial[key] = initial_prob[key]/float(sum)
    initial.update({'M':0})
    return initial

#Calculating Transition probability

def trans_prob(pos):
	transition_prob = {}
	mem = 0
	inner = 0
	outer = 0
	for y in pos:
    		for k in range(len(y)-1):
        		if y[k] in 'iMo' and y[k+1] in 'iMo' and y[k] != y[k+1]:
            			state = y[k] + y[k+1]
            			if state in transition_prob:
                			transition_prob[state] += 1
            			else :
                			transition_prob[state] = 1
           			if y[k] == 'i': 
                			inner += 1 
            			elif y[k] == 'M':
                			mem +=1 
            			else:
                			outer += 1
	for k in transition_prob:
    		if k[0] == 'i':
        		transition_prob[k] = float(transition_prob[k])/inner
    		elif k[0] == 'M':
        		transition_prob[k] = float(transition_prob[k])/mem
   		else:
        		transition_prob[k] = float(transition_prob[k])/outer
	transition_prob.update({'ii':0})
	transition_prob.update({'MM':0})
	transition_prob.update({'oo':0})
	transition_prob.update({'oi':0})
	transition_prob.update({'io':0})
	return transition_prob

#calculating emission probability

def emis_prob(protein,pos):
	emission_prob = {}
	inner = 0
	mem = 0
	outer = 0
	for x,y in zip(protein,pos):
	    for i in range(len(x)):
	        if y[i] in 'iMo':
        	    state = y[i] + x[i]
        	    if state in emission_prob:
        	        emission_prob[state] += 1
        	    else:
        	        emission_prob[state] = 1
        	    if y[i] == 'i':
        	        inner += 1
        	    elif y[i] == 'M':
        	        mem += 1
        	    else:
        	        outer += 1
	for k in emission_prob:
	    if k[0] == 'i':
        	emission_prob[k] = float(emission_prob[k])/inner
	    elif k[0] == 'M':
        	emission_prob[k] = float(emission_prob[k])/mem
	    else:
	        emission_prob[k] = float(emission_prob[k])/outer
	return emission_prob

def length_d(pos):
    outer = []
    inner = []
    men = []
    for x in pos:
        count_o = 0
        count_M = 0
        count_i = 0
        x = list(x)
        for i in range(len(x) - 1):
            if x[i] == 'o' and x[i+1] == 'o':
                count_o += 1  
            elif x[i] == 'o' and x[i+1] != 'o':
                count_o += 1
                outer.append(count_o)
                count_o = 0
            elif x[i] == 'M' and x[i+1] == 'M':
                count_M += 1
            elif x[i] == 'M' and x[i+1] != 'M':
                count_M += 1
                men.append(count_M)
                count_M = 0
            elif x[i] == 'i' and x[i+1] == 'i':
                count_i += 1
            elif x[i] == 'i' and x[i+1] != 'i':
                count_i += 1
                inner.append(count_i)
                count_i = 0    
    return outer, inner, men



def backtrack(test, emission_prob, pro_matrix, pos_matrix):
    length = len(test) - 1
    annotation = ''
    score = max(pro_matrix[0][length],pro_matrix[1][length],pro_matrix[2][length])
    if score == pro_matrix[0][length]:
        tmp = pos_matrix[0][length]
        state = 'M'
    elif score == pro_matrix[1][length]:
        tmp = pos_matrix[1][length]
        state = 'i'
    else:M
        tmp = pos_matrix[2][length]
        state = 'o'
    repeat = length - tmp[1]
    annotation += state*repeat
    pre_state = tmp[0]
    pre_pos = tmp[1]
    
    while True:
        tmp = pos_matrix[pre_state][pre_pos]
        if pre_state == 0:
            state = 'M'
        elif pre_state == 1:
            state = 'i'
        else:
            state = 'o'
        repeat = pre_pos - tmp[1]
        annotation += state*repeat
        pre_state = tmp[0]
        pre_pos = tmp[1]
        if pre_pos == 0:
            annotation += state
            break
            
    annotation = annotation[::-1]
    return annotation
    
def print_outcome(ghmm):
    print "    ",
    for i in range(len(ghmm)): print "%7d" % i,
    print

    for y in ghmm[0].keys():
        print "%.5s: " % y,
        for t in range(len(ghmm)):
            print "%.7s" % ("%f" % ghmm[t][y]),
        print


	
	
if len(sys.argv) != 5:
		print "USAGE: python prodictMP_ghmm.py -i [INPUT_FILE] -o [OUTPUT_FILE]"
		exit(-1)

if sys.argv[1] == "-i":
		if sys.argv[2].endswith('.ffa'):
			input_file = sys.argv[2]
		else:
			print "The input file should be in ffa format"
			exit(-1)
			
if sys.argv[3] == "-o":
    output_file = sys.argv[4]  
else:
		print "USAGE: python prodictMP_ghmm.py -i [INPUT_FILE] -o [OUTPUT_FILE]"
		exit(-1)
		
output = open(output_file, 'w')
protein,pos = preprocessing(input_file)
initial_prob = initialprob(pos)
transition_prob = trans_prob(pos)
emission_prob = emis_prob(protein,pos)


#use line1 protein sequence as test
test = list("MENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA")
lengths = len(test)/10 + 1
lengthss = len(test)/2 + len(test)%2 
#get length duration
#Call the length_d function and hist each state using R studio 
#get the frequency data of outer in bin 10 range
#get the frequency data of inner in bin 10 range
#get the frequency data of Membrane in bin 2 range
#length_dura_o, length_dura_i, length_dura_M = length_d(pos) 
outer_l = [0.35616438, 0.27397260, 0.15068493, 0.07762557, 0.04109589, 0.04109589, 
0.00913242, 0.00456621, 0.00456621, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 
0.00456621, 0.00000001, 0.00456621, 0.00000001, 0.00456621, 0.00000001, 0.00456621, 
0.00456621, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00913242, 0.00000001, 
0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00913242]

inner_l = [0.265873016, 0.349206349, 0.222222222, 0.083333333, 0.027777778, 0.019841270, 
0.003968254, 0.007936508, 0.000000001, 0.003968254, 0.000000001, 0.000000001, 0.000000001, 
0.000000001, 0.000000001, 0.003968254, 0.003968254, 0.000000001, 0.000000001, 0.000000001, 
0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 
0.000000001, 0.000000001, 0.003968254, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 
0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 
0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.003968254]

Mem_l = [0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 
0.008810573, 0.022026432, 0.090308370, 0.235682819, 0.359030837, 0.147577093, 0.103524229, 
0.017621145, 0.006607930, 0.000000001, 0.006607930, 0.000000001, 0.000000001, 0.002202643]

if len(outer_l) <= lengths:
    for i in range(len(outer_l),lengths):
        outer_l.append(0.000000001)
if len(inner_l) <= lengths:
    for i in range(len(inner_i),lengths):
        inner_l.append(0.000000001)
if len(Mem_l) <= lengthss:
    for i in range(len(Mem_l),lengthss):
        Mem_l.append(0.000000001)        
length = {}
length['o'] = outer_l
length['i'] = inner_l
length['M'] = Mem_l
#end of length duration

output.write()



#print protein
#print pos
#print initial_prob
#print transition_prob
#print emission_prob
#print "All outer annotation length:"
#print length_dura_o
#print "All inner annotation length:"
#print length_dura_i
#print "All Membrane annotation length:"
#print length_dura_M




