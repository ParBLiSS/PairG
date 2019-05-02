#
# @file    dot2txt.py
# @author  Chirag Jain <cjain7@gatech.edu>
# @brief   convert dot formatted compressed DBG from splitMEM to our .txt format
#            note: splitMEM.cc was modified to print vertex DNA labels as well
#                  (set OPT_DisplaySeq=1 and OPT_SeqToDisplay= <int>::max() )           
#            note: After reading code in splitMEM, we realized that it replaces 'N' 
#                  with 'A' character, and it builds suffix tree while merging adjacent 
#                  fasta strings with a 'N' character in between
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
      label = re.sub(r'[^ACGTN]', '', tokens[1])
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
  trimVertexLabel[v] = 1

for i in range(0, len(vertexLabels)):
  if trimVertexLabel[i] == 1:
    vertexLabels[i] = vertexLabels[i][k-1:]

#####
# Compute out-neighbors of each vertex
####

edge_index = {}
for (u,v) in edges:
  edge_index[u] = []

for (u,v) in edges:
  edge_index[u].append(v)

####
# Finally print the contents
####

print len(vertexLabels)
for u in range(0, len(vertexLabels)):
  if u in edge_index:
    for v in edge_index[u]:
      print v,
  print vertexLabels[u]


