import bct as bctDevV2 
import numpy as np 
from tqdm import tqdm 
from struct import pack 
import sys 

#####################################################################################################################

infile='fof_group_SV_lgal_SOD/trees_099.0.vector' 
outfile='mergedLJ.0'

infile='../../lgalMill/MergerTrees/MR/treedata/trees_063.4'
outfile='mergedMill500_2.0'
trees=bctDevV2.read_binary(infile) 
 
branchingTrees = [] 
straightTrees = []
for i in range(len(trees)): # Iterate through all trees
    if len(trees[i]):
        if trees[i][0]['Len'] > 0:
            if (trees[i]['NextHaloInFOFGroupOffset']>-1).any():    # If at any point there are multiple subhalos in group
                branchingTrees.append(i)    # Add the current tree to the list of trees to be fed through the merging code
            else:   # Otherwise we have a 'straight' tree and do not need to merge
                straightTrees.append(i)    # Add the current tree to list of trees to be sent straight to output

#####################################################################################################################

print("Num trees =", len(branchingTrees)) 
outputTrees = []    # Array of trees to be rewritten to new vector file

for treenum in tqdm(straightTrees): # Iterate through 'straight' trees
    outputTrees.append(trees[treenum].copy())   # Add to output array without merging
    
for treenum in tqdm(branchingTrees):    # Iterate through 'branching' trees
    temp = trees[treenum].copy()    # Create temporary copy of current tree
    debug=False 

    if debug: print(temp) 
    newBranchStart=[]  
    newBranchEnd=[] 
    for index in range(len(temp)-1):    # Iterate through current tree
        if temp[index]['FirstProgenitorOffset']==-1:    # When first progenitor offset==-1, we are at the end of the current branch
            newBranchStart.append(index+1)  # Append the index of the next core to a 'newBranchStart' array
            newBranchEnd.append(index) # Append the current index to a 'newBranchEnd' array
    newBranchEnd.pop(0) # Pop the first newBranchEnd element, as this is the end of the main progenitor line
    newBranchEnd.append(len(temp)-1)    # Append the length of the tree-1 to the newBranchEnd array, as this is the end of the final branch
     
    breakFlag=False
    delIndices = [] # Satellites to be deleted  
    delIndicesDict = {} 
    minLen = 500 # Min length of satellites   
    for index in range(len(newBranchEnd)):    # Iterate through end of branches
         
        fhTBD = [] 
        current=newBranchEnd[index] # Current halo in branch (begin at end)  
        start=newBranchStart[index] # First halo in branch  
        while current!=start-1:   # Move through branch towards start (while current!=start-1), this loop is used to prepare the tree for the removal of satellites
            if temp[current]['Len']<minLen and temp[current]['M_Crit200']==0.:   # If current halo<mass limit and not a central   
                if current!=newBranchEnd[index]:    # If current is not the end of branch   
                    temp[current+1]['DescendentOffset']=temp[current]['FirstHaloInFOFGroupOffset']  # Set progen's desc to core's central
                temp['NextHaloInFOFGroupOffset'] = np.where(temp['NextHaloInFOFGroupOffset']==current, temp[current]['NextHaloInFOFGroupOffset'], temp['NextHaloInFOFGroupOffset']) # Find any instances of current core as set as nextHalo to current core's nextHalo  
                
                firstHalos = np.where(temp['FirstHaloInFOFGroupOffset']==current)[0] # First halo not central edge case
                if firstHalos.size:  
                    breakFlag=True 
                
                delIndices.append(current)  # Add current core to list of cores to be deleted, can't delete immediately because we have to keep indexing until wehave iterate through all branches
            current-=1  # Move to next core in the current branch  
        
        # Skip FH edge case
        if breakFlag: break
        
        temp2=temp.copy() 
        current=newBranchEnd[index] 
        start=newBranchStart[index] 
        while current!=start-1:     # After finished iterating through branch (current==start-1) to flag satellites below mass cut for deletion, loop through branch once more to update the indexes of cores that are not going to be deleted (this loop must be done after the previous pass through the branch has completed because we need to know the number of cores needing deletion in order to update the indexes)
            
            if current not in delIndices:   # Ieterate through branch once again, if a core was not flagged for deletion
                fhMatch = np.where(temp2['FirstHaloInFOFGroupOffset']==start)[0]  # Find all instances of start core in tree (firstHaloInFOFGroup, nextHaloInFOFGroup, Descendent)
                dMatch = np.where(temp2['DescendentOffset']==start)[0]
#                fPMatch = np.where(temp2['FirstProgenitorOffset']==start)[0]

                if fhMatch.size: temp['FirstHaloInFOFGroupOffset'][fhMatch] = start-len(delIndices)   # Replace all instances of start core's index by index - # cores to be deleted                
                if dMatch.size: temp['DescendentOffset'][dMatch] = start-len(delIndices)
 #               if fPMatch.size: temp['FirstProgenitorOffset'][fPMatch] = start-len(delIndices)

                if temp[start]['FirstProgenitorOffset'] != -1: # Do the same for firstProgenitor  
                    temp[start]['FirstProgenitorOffset'] = start-len(delIndices)+1  # But instead replace by index - # cores to be deleted + 1 to abide by indexing conventions            
            current-=1 
    if breakFlag: continue

    temp = np.delete(temp, delIndices) # At this point we can safely delete the satellites below the mass cut from the tree
 
    if debug: print("\n",temp) 
    for i in range(len(temp)):  # Iterate through tree to set NextProgenitors (by default set to -1 for all cores, this loop updates this property after merging)
        progs = sorted(temp[temp['DescendentOffset']==i], key=lambda tree:tree['Len'], reverse=True)    # Identify all progenitors of current core, sorted by mass
 
        if len(progs)>1:    # If more than one progenitor exists
            firstProg = temp[i]['FirstProgenitorOffset']
            nextProgs = [] 
            for prog in progs:  # Iterate through the progenitors in this list
                if np.where(temp==prog)[0][0]==firstProg: continue  # If the current progenitor is the firstProgenitor, skip
                else: nextProgs.append(prog)    # Otherwise, append to a list of nextProgenitors
 
            temp[firstProg]['NextProgenitorOffset'] = np.where(temp==nextProgs[0])[0][0]    # Set firstProgenitor's nextProgenitor to the most massive nextProgen
            progIndices = [] 
            for i in range(len(nextProgs)-1):   # Then iterate through all other nextProgenitors
                currentIndex = np.where(temp==nextProgs[i])[0][0]   
                progIndex = np.where(temp==nextProgs[i+1])[0][0] 
                temp[currentIndex]['NextProgenitorOffset'] = progIndex # Set current nextProgenitor's nextProgenitor to the next nextPrognenitor in the list of nextProgenitors :) 

    if debug: print("\n",temp) 
    outputTrees.append(temp)    # Append merged tree to the list of output trees

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
