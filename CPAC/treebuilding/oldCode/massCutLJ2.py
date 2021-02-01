import bct as bctDevV2
import numpy as np
from tqdm import tqdm
from struct import pack
import sys

#####################################################################################################################

infile='mergedLJ_.0'
outfile='cutLJ_.0'
trees=bctDevV2.read_binary(infile)

count = 0
branchingTrees = []
straightTrees = []
outputTrees = []

cutLen=500
cutMass=65
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
        elif temp[index]['M_Crit200']>0 and temp[index]['M_Crit200']<cutMass:
            cutFlag=True
            delIndices.append(index)
    temp = np.delete(temp, delIndices) 
    
    if len(temp):
        temp['NextProgenitorOffset'] = -1
        temp[len(temp)-1]['FirstProgenitorOffset'] = -1
        count+=1
        outputTrees.append(temp)  

#####################################################################################################################

print(count)
count = 0
# Iterate through 'branching' tree
for i, treenum in enumerate(tqdm(branchingTrees)):  
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

    # Iterate through end of branches
    delIndices = []
    temp2 = temp.copy()
    for index in range(len(newBranchEnd)):    
        
        end=newBranchEnd[index]
        start=newBranchStart[index]
        while start!=end+1:
            if len(delIndices):
                fhMatch = np.where(temp2['FirstHaloInFOFGroupOffset']==start)[0]
                dMatch = np.where(temp2['DescendentOffset']==start)[0]
                fPMatch = np.where(temp2['FirstProgenitorOffset']==start)[0]

                if fhMatch.size: temp['FirstHaloInFOFGroupOffset'][fhMatch] = start-len(delIndices)
                if dMatch.size: temp['DescendentOffset'][dMatch] = start-len(delIndices)
                if fPMatch.size: temp['FirstProgenitorOffset'][fPMatch] = start-len(delIndices)
            start+=1
        
        cutFlag=False
        end=newBranchEnd[index] 
        start=newBranchStart[index] 
        while start!=end+1:
            if cutFlag==True:   
                temp['NextProgenitorOffset'] = np.where(temp['NextProgenitorOffset']==start, temp[start]['NextProgenitorOffset'], temp['NextProgenitorOffset']) 
                delIndices.append(start)
            elif temp[start]['Len']<cutLen:   
                cutFlag=True
                temp['NextProgenitorOffset'] = np.where(temp['NextProgenitorOffset']==start, temp[start]['NextProgenitorOffset'], temp['NextProgenitorOffset']) 

                if start!=newBranchStart[index]:
                    temp[start-1]['FirstProgenitorOffset'] = -1
                delIndices.append(start)  
            
            '''
            elif temp[start]['M_Crit200']>0 and temp[start]['M_Crit200']<cutMass:
                cutFlag=True
                temp['NextProgenitorOffset'] = np.where(temp['NextProgenitorOffset']==start, temp[start]['NextProgenitorOffset'], temp['NextProgenitorOffset'])
                
                if temp[start]['FirstHaloInFOFGroupOffset']==start: fhMatch = np.where(temp2['FirstHaloInFOFGroupOffset']==start)[0]
                if fhMatch.size: 
                    for element in fhMatch:
                        delIndices.append(element)

                if start!=newBranchStart[index]:
                    temp[start-1]['FirstProgenitorOffset'] = -1
                delIndices.append(start)
            '''
            start+=1
        
    temp = np.delete(temp, delIndices)
    
    #####################################################################################################################
    
    #print(count)
    try:
        for i in range(len(temp)):  # Iterate through tree to set NextProgenitors (by default set to -1 for all cores, this loop updates this property after merging)
            progs = sorted(temp[temp['DescendentOffset']==i], key=lambda tree:tree['Len'], reverse=True)    # Identify all progenitors of current core, sorted by mass
            if len(progs)>1:    # If more than one progenitor exists
                firstProg = temp[i]['FirstProgenitorOffset']
                nextProgs = []
                for prog in progs:  # Iterate through the progenitors in this list
                    if np.where(temp==prog)[0][0]==firstProg: continue  # If the current progenitor is the firstProgenitor, skip
                    else: nextProgs.append(prog)    # Otherwise, append to a list of nextProgenitors
                
                if np.where(temp==nextProgs[0])[0].size: temp[firstProg]['NextProgenitorOffset'] = np.where(temp==nextProgs[0])[0][0]    # Set firstProgenitor's nextProgenitor to the most massive nextProgen
                progIndices = []
                for i in range(len(nextProgs)-1):   # Then iterate through all other nextProgenitors
                    if np.where(temp==nextProgs[i])[0].size: currentIndex = np.where(temp==nextProgs[i])[0][0]
                    if np.where(temp==nextProgs[i+1])[0].size: progenIndex = np.where(temp==nextProgs[i+1])[0][0]
                    temp[currentIndex]['NextProgenitorOffset'] = progenIndex # Set current nextProgenitor's nextProgenitor to the next nextPrognenitor in the list of nextProgenitors :)
    
    #####################################################################################################################

        for i in range(len(temp)):
            nextHalos = sorted(temp[temp['FirstHaloInFOFGroupOffset']==i], key=lambda tree:tree['Len'], reverse=True)    # Identify all progenitors of current core, sorted by mass
    
            if len(nextHalos)>1:    # If more than one halo in group
                firstHalo = i
                nHalos = []
                for halo in nextHalos:  # Iterate through the progenitors in this list
                    if np.where(temp==halo)[0][0]==i: continue  # If the current progenitor is the firstProgenitor, skip
                    else: nHalos.append(halo)    # Otherwise, append to a list of nextProgenitors

                if np.where(temp==nHalos[0])[0].size: temp[firstHalo]['NextHaloInFOFGroupOffset'] = np.where(temp==nHalos[0])[0][0]    # Set firstProgenitor's nextProgenitor to the most massive nextProgen
                nhIndices = []
                for i in range(len(nHalos)-1):   # Then iterate through all other nextProgenitors
                    if np.where(temp==nHalos[i])[0].size: currentIndex = np.where(temp==nHalos[i])[0][0]
                    if np.where(temp==nHalos[i+1])[0].size: nextHaloIndex = np.where(temp==nHalos[i+1])[0][0]
                    temp[currentIndex]['NextHaloInFOFGroupOffset'] = nextHaloIndex # Set current nextProgenitor's nextProgenitor to the next nextPrognenitor in the list of nextProgenitors :)
    except:
        continue
    count+=1
    outputTrees.append(temp)    # Append merged tree to the list of output trees
print(count)

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
