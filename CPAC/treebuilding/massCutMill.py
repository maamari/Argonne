import bct as bctDevV2
import numpy as np
from tqdm import tqdm
from struct import pack
import sys

#####################################################################################################################

infile='binaries/mergedMill500.0'
outfile='cutMill5003.0'

trees=bctDevV2.read_binary(infile)

count = 0
branchingTrees = []
straightTrees = []
outputTrees = []
cutLen=500
cutMass=1
for i in range(len(trees)):
    if (len(trees[i])):
        if (trees[i]['NextHaloInFOFGroupOffset']>-1).any():  
            branchingTrees.append(i)
        else:   
            straightTrees.append(i) 
print("Num trees =", len(straightTrees))

#####################################################################################################################

# Iterate through 'straight' trees
for treenum in tqdm(straightTrees):
    temp = trees[treenum].copy()
    
    delIndices = []
    cutFlag=False
    for index in range(len(temp)):
        if cutFlag:
            delIndices.append(index)
        elif temp[index]['Len']<cutLen:
            cutFlag=True
            delIndices.append(index)
        elif temp[index]['M_Crit200']<=cutMass:
            cutFlag=True
            delIndices.append(index)
    temp = np.delete(temp, delIndices) 
    
    if len(temp):
        temp['NextProgenitorOffset'] = -1
        temp[len(temp)-1]['FirstProgenitorOffset'] = -1
        count+=1
        outputTrees.append(temp)  

#####################################################################################################################

# Iterate through 'branching' tree
for i, treenum in enumerate(tqdm(branchingTrees)):  
    if i%2==0 or i%7==0: continue
    temp = trees[treenum].copy()  
    if len(temp)<=1: continue
    temp['NextProgenitorOffset']=-1
    temp['NextHaloInFOFGroupOffset']=-1

    #####################################################################################################################    
    
    newBranchStart=[]
    newBranchEnd=[]
    newBranchStart.append(0)
    for index in range(len(temp)-1):
        if temp[index]['FirstProgenitorOffset']==-1: 
            newBranchStart.append(index+1)  
            newBranchEnd.append(index) 
    newBranchEnd.append(len(temp)-1)   

    #####################################################################################################################

    delIndices = []
    for index in range(len(temp)):
        if index>=newBranchEnd[0]:
            delIndices.append(index)
    temp = np.delete(temp, delIndices)
   
    delIndices = []
    cutFlag=False
    for index in range(len(temp)):
        if cutFlag:
            delIndices.append(index)
        elif temp[index]['Len']<cutLen:
            cutFlag=True
            delIndices.append(index)
        elif temp[index]['M_Crit200']<=cutMass:
            cutFlag=True
            delIndices.append(index)            
    temp = np.delete(temp, delIndices)

    for index in range(len(temp)):
        temp[index]['FirstHaloInFOFGroupOffset']=index

    if len(temp):
        temp['NextProgenitorOffset'] = -1
        temp[len(temp)-1]['FirstProgenitorOffset'] = -1
        count+=1
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
