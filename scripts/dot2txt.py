#
# @file    dot2txt.py
# @author  Chirag Jain <cjain7@gatech.edu>
# @brief   convert dot formatted compressed DBG from splitMEM to our .txt format
#            note: splitMEM.cc was modified to print vertex DNA labels as well
#                  (set OPT_DisplaySeq=1 and OPT_SeqToDisplay= <int>::max() )           
#            note: After reading code in splitMEM, we realized that it replaces each 'N' 
#                  character with 'A', and it builds the concatenated string by merging  
#                  adjacent fasta records while putting a 'N' character in between, and 
#                  at the very end
# @usage   python dot2txt.py "MEM or kmer size" "dot file"  >  ".txt file"

import sys
import re

k=int(sys.argv[1])
dotFile=sys.argv[2]

vertexLabels = []
edges = []

#######
# Read all vertex labels and directed edges
######

with open(dotFile) as fp:  
  currentVertexId = -1
  for line in fp:
    tokens = line.split()
    if len(tokens) == 2: # must be a vertex label
      currentVertexId += 1
      assert int(tokens[0]) == currentVertexId
      label = re.sub(r'[^ACGTN\$]', '', tokens[1])  #get rid of ambiguous character
      label = re.sub(r'[\$]', 'N', label) # get rid of ambiguous characters 
      vertexLabels.append(label)
    elif len(tokens) == 3:
      if tokens[1] == "->": # must be an edge
        assert int(tokens[0]) == currentVertexId
        edges.append((int(tokens[0]), int(tokens[2])))

#Remove duplicates 
edges = list(dict.fromkeys(edges))

#####
# Remove overlapping prefix from selected vertices
# Algorithm:  if a vertex has in-degree > 0, then 
#             remove first (k-1) characters from its label
#####

trimVertexLabel = [0] * len(vertexLabels)
for (u,v) in edges:
  assert u >= 0
  assert u < len(vertexLabels)
  assert v >= 0
  assert v < len(vertexLabels)
  trimVertexLabel[v] = 1

for i in range(0, len(vertexLabels)):
  if trimVertexLabel[i] == 1:
    vertexLabels[i] = vertexLabels[i][k-1:]

#####
# Compute out-neighbors of each vertex
####

out_index = {}
for (u,v) in edges:
  out_index[u] = []

for (u,v) in edges:
  out_index[u].append(v)

####
# Finally print the contents
####

print len(vertexLabels)
for u in range(0, len(vertexLabels)):
  if u in out_index:   #print out-neighbors
    for v in out_index[u]:
      print v,
  assert len(vertexLabels[u]) > 0
  print vertexLabels[u] #print label
