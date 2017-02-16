#!/usr/bin/python
# Importing libraries
import sys
import numpy as np


#define the script command line
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "USAGE: python pruning.py tree_file alignment_file"
        exit(-1)
    else:
        tree_file = sys.argv[1]
        alignment_file = sys.argv[2]


#Read and parse the information from the alignment file
compare_species = {} 
f1 = file(alignment_file)
for line in f1:
    line = line.strip("\n")
    species = line.split(":",1)[0]
    seq = line.split(":",1)[1]
    compare_species[species] = list(seq)
print compare_species	


# Newick format for phylogenetic tree
# Read and parse the information from the tree file
f2 = file(tree_file)
for line in f2:
    phylo_tree = line
tree = list(phylo_tree)
print tree

# Following is the P(0.1) substitution matrix
p01 = np.array( ((0.9, 0.05, 0.025, 0.025), 
                 (0.05, 0.9, 0.025, 0.025), 
                 (0.025, 0.025, 0.9, 0.05), 
                 (0.025, 0.025, 0.05, 0.9)) )
print "The P(0.1) substitution matrix is: "
print p01

# Building substitution matrix based on matrix multiplication
def matrix_mutiplication(matrix1, matrix2):
  new_matrix = np.dot(matrix1, matrix2)
  return new_matrix

# Computing p(0.2) and p(0.3) substitution matrix
p02 = matrix_mutiplication(p01,p01)
print "The P(0.2) substitution matrix is: "
print p02
p03 = matrix_mutiplication(p02,p01)
print "The P(0.3) substitution matrix is: "
print p03

def compute_prob(left_child, right_child, left_branch, right_branch):
  parent = [0, 0, 0, 0]
  
  if left_branch == '0.1':
    left_matrix = p01
  elif left_branch == '0.2':
    left_matrix = p02
  else:
    left_matrix = p03

  if right_branch == '0.1':
    right_matrix = p01
  elif right_branch == '0.2':
    right_matrix = p02
  else:
    right_matrix = p03

  for j in range(4): 
    #compute the parent[i] probability
    left = 0
    right = 0
    for i in range(4):
      left = left + left_child[i] * left_matrix[i][j]
    for i in range(4):
      right = right + right_child[i] * right_matrix[i][j]
    parent[j] = left * right

  return parent  

def new_species(nucleotide):
  if nucleotide == 'A':
    leaf_list = [1,0,0,0]
  elif nucleotide == 'C':
    leaf_list = [0,1,0,0]
  elif nucleotide == 'G':
    leaf_list = [0,0,1,0]
  else:
    leaf_list = [0,0,0,1]
  return leaf_list

def pruning(tree, species_list, result_alignment, result_prob):   # the format of species_list [ ['species_name', [A,C,G,T]], ... ]
  index = -1
  compare_list = []
  branch_list = []
  while True:
    index = index + 1
    #print index
    #print compare_list
    #print branch_list
    if tree[index] == ';':
      break

    if tree[index] == '(':
      if tree[index+1] != '(':
        index = index + 1
        name = ''
        while tree[index] != ':':
          name = name + tree[index]
          index = index + 1
        for element in species_list:
          if element[0] == name:
            compare_list.append(element[1])
        index = index - 1

    if tree[index] == ':':
      index = index + 1
      branch = ''
      while tree[index] != ',' and tree[index] != ')':
        branch = branch + tree[index]
        index = index + 1
      branch_list.append(branch)
      index = index - 1
      #print branch_list

    if tree[index] == ',':
      if tree[index+1] != '(':
        index = index + 1
        name = ''
        while tree[index] != ':':
          name = name + tree[index]
          index = index + 1
        for element in species_list:
          if element[0] == name:
            compare_list.append(element[1])
        index = index - 1

    if tree[index] == ')':
        right_child = compare_list.pop()
        right_branch = branch_list.pop()
        left_child = compare_list.pop()
        left_branch = branch_list.pop()
        compare_list.append(compute_prob(left_child, right_child, left_branch, right_branch))
        #print compare_list

  #print compare_list
  find_max = 0
  get_nucle = 0
  for i in range(4):
    if compare_list[0][i] > find_max:
      find_max = compare_list[0][i]
      get_nucle = i
  result_prob.append(find_max)
  if get_nucle == 0:
    result_alignment.append('A')
  elif get_nucle == 1:
    result_alignment.append('C')
  elif get_nucle == 2:
    result_alignment.append('G')
  else:
    result_alignment.append('T')

def iteration(compare_species, tree):
  result_alignment = []
  result_prob = []
  tmp = compare_species.values()
  for i in range(len(tmp[0])):
    species_list = []
    for key in compare_species.keys():
      info = []
      info.append(key)
      info.append(new_species(compare_species[key][i]))
      species_list.append(info)
      #print species_list
    pruning(tree, species_list, result_alignment, result_prob)
  return result_alignment, result_prob


result_alignment, result_prob = iteration(compare_species, tree)
print "The best alignment is: "
print result_alignment
print "The best probability for each column is: "
print result_prob
print "The total probability is: "
total = 1
for i in range(len(result_prob)):
    total = total * result_prob[i]
print total