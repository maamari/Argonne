import bctDevV2 
import numpy as np 
from tqdm import tqdm 
from struct import pack 
import sys 

#outfile='../../lgalLarge/MergerTrees/LJDS/treedata/trees_099.0'
#outfile='massiveTrees_099.0'
infile='../../lgalMill/MergerTrees/MR/treedata/_trees_063.5'
outfile='cutMill2.0'
trees=bctDevV2.read_binary(infile) 
 
massiveTrees = []
count=0
for i in tqdm(range(len(trees))):
    #if trees[i][0]['Len'] > 185:
    #    count+=1
    #    if count%7==0:
    #        massiveTrees.append(trees[i]) 
    if len(trees[i]):
        if trees[i][0]['Len'] > 500:
            massiveTrees.append(trees[i])

print("Num trees =", len(massiveTrees)) 

struct_keys = ['DescendentOffset', 'FirstProgenitorOffset', 'NextProgenitorOffset','FirstHaloInFOFGroupOffset', 'NextHaloInFOFGroupOffset','Len', 'M_Mean200', 'M_Crit200', 'M_TopHat','Pos_x', 'Pos_y', 'Pos_z', 'Vel_x', 'Vel_y', 'Vel_z','VelDisp', 'Vmax', 'Spin_x', 'Spin_y', 'Spin_z','MostBoundID->Coretag', 'SnapNum', 'FileNr','SubhaloIndex', 'SubHalfMass']
struct_format = "<iiiiiiffffffffffffffqiiif"
header_format = "<{}i"

with open(outfile, "wb") as fh:
    NTrees = len(massiveTrees)
    totNHalos = 0
    for i in range(len(massiveTrees)):
        totNHalos += len(massiveTrees[i])
    values = [NTrees, totNHalos]
    for i in range(len(massiveTrees)):
        values.append(len(massiveTrees[i]))
    header_this = header_format.format(NTrees + 2)
    fh.write(pack(header_this, *values))
    print('Wrote header for {} fof massiveTrees with {} halos to {}'.format(NTrees, totNHalos, outfile))

    for i in range(len(massiveTrees)):
        for j in range(len(massiveTrees[i])):
            values = list(massiveTrees[i][j])
            fh.write(pack(struct_format, *values))

