import bctDevV2
import numpy as np
from tqdm import tqdm
from struct import pack
import sys

#####################################################################################################################

infile='cut.0'
outfile='cut2.0'
trees=bctDevV2.read_binary(infile)

print("Num trees =", len(trees))

#####################################################################################################################

outputTrees=[]
# Iterate through 'straight' trees
for treenum in tqdm(trees):
    temp = trees[treenum].copy()
    
    temp['M_Crit200'] *= 2.7e9/1e10    
    outputTrees.append(temp)  

#####################################################################################################################

struct_keys = ['DescendentOffset', 'FirstProgenitorOffset', 'NextProgenitorOffset','FirstHaloInFOFGroupOffset', 'NextHaloInFOFGroupOffset','Len', 'M_Mean200', 'M_Crit200', 'M_TopHat','Pos_x', 'Pos_y', 'Pos_z', 'Vel_x', 'Vel_y', 'Vel_z','VelDisp', 'Vmax', 'Spin_x', 'Spin_y', 'Spin_z','MostBoundID->Coretag', 'SnapNum', 'FileNr','SubhaloIndex', 'SubHalfMass']
struct_format = "<iiiiiiffffffffffffffqiiif"
header_format = "<{}i"

with open(outfile, "wb") as fh: # Write list of output trees to output file using Eve's tree structure
    NTrees = len(outputTrees)
    totNHalos = 0
    for i in range(len(outputTrees)):
        totNHalos += len(outputTrees[i])
    values = [NTrees, totNHalos]
    for i in range(len(outputTrees)):
        values.append(len(outputTrees[i]))
    header_this = header_format.format(NTrees + 2)
    fh.write(pack(header_this, *values))
    print('Wrote header for {} fof outputTrees with {} halos to {}'.format(NTrees, totNHalos, outfile))

    for i in range(len(outputTrees)):
        for j in range(len(outputTrees[i])):
            values = list(outputTrees[i][j])
            fh.write(pack(struct_format, *values))

#####################################################################################################################
