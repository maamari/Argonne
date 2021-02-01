import bctDevV2
import numpy as np
from tqdm import tqdm
from struct import pack
import sys

#####################################################################################################################

infile='../../lgalLarge/MergerTrees/LJDS/treedata/_trees_099.0'
outfile='../../lgalLarge/MergerTrees/LJDS/treedata/trees_099.0'
#infile='cutMill5002.0'
#outfile='../../lgalMill/MergerTrees/MR/treedata/trees_063.5'
trees=bctDevV2.read_binary(infile)

print("Num trees =", len(trees))

#####################################################################################################################

outputTrees=[]
# Iterate through 'straight' trees
for treenum in tqdm(trees):
    temp = trees[treenum].copy()
    
    for num in range(len(temp)):
        if temp[num]['M_Crit200']!=0.0:
            #temp[num]['Len']=1e-10
            #temp[num]['M_Crit200']=1e-10
            #temp[num]['Vmax']=1e-10
            #temp[num]['VelDisp']=1e-10
            temp[num]['Spin_x']=1e-10
            temp[num]['Spin_y']=1e-10
            temp[num]['Spin_z']=1e-10
            #temp[num]['Vel_x']=1e-10

    #temp['M_Crit200'] *= 2.7e9/1e10   
    #temp['M_Crit200'] = np.where(temp['M_Crit200']<500, 500, temp['M_Crit200'])
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
