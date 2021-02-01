import bctDevV2 
import numpy as np 
from tqdm import tqdm 
from struct import pack 
import sys 

#####################################################################################################################

infile='trees_099.0.vector' 
outfile='mergedLJ2.0'

#infile='../../lgalMill/MergerTrees/MR/treedata/trees_063.5'
#outfile='mergedMill.0'
trees=bctDevV2.read_binary(infile) 
 
count=0
for i in range(len(trees)):
    if len(trees[i]):
        if trees[i][0]['Len'] > 185:
            count+=1
cutFactor=count//1000

branchingTrees = [] 
straightTrees = []
count=0
for i in range(len(trees)): # Iterate through all trees
    if len(trees[i]):
        if trees[i][0]['Len'] > 500:    # If tree satisfies some minimum z=0 length (optional, set to > 0 if you do not want to use minimum mass cut)
            if (trees[i]['NextHaloInFOFGroupOffset']>-1).any():    # If at any point there are multiple subhalos within a FOF group
                branchingTrees.append(i)    # Add the current tree to the list of trees to be fed through the merging code
            else:   # Otherwise we have a 'straight' tree and do not need to merge
                straightTrees.append(i)    # So add the current tree to a separate list of trees to be sent straight to output without merging

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
    minLen = 500  # Min length of satellites   
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
                    firstHalos = firstHalos[firstHalos!=current]                     
                    firstHalos = np.append(firstHalos,current) 
                    for halo in firstHalos: fhTBD.append(halo) 
                
                delIndices.append(current)  # Add current core to list of cores to be deleted, can't delete immediately because we have to keep indexing until wehave iterate through all branches
            current-=1  # Move to next core in the current branch  
        
        # Skip FH edge case
        if breakFlag: break
        
        temp2=temp.copy() 
        current=newBranchEnd[index] 
        start=newBranchStart[index] 
        while current!=start-1:     # After finished iterating through branch (current==start-1) to flag satellites below mass cut for deletion, loop through branch once more to update the indexes of cores that are not going to be deleted (this loop must be done after the previous pass through the branch has completed because we need to know the number of cores needing deletion in order to update the indexes)
            
            # FH edge case
            delIndicesDict[current]=len(delIndices)
            if fhTBD:  
                badHalo = fhTBD[-1]  
                fhTBD.remove(badHalo)  
                sortedFH = sorted(temp[fhTBD], key=lambda tree:tree['M_Crit200'], reverse=True)   
                fhMatch = np.where(temp['FirstHaloInFOFGroupOffset']==badHalo)[0]  
                fhMatch = fhMatch[fhMatch!=badHalo]  
                temp['FirstHaloInFOFGroupOffset'][fhMatch]=np.where(temp==sortedFH[0])[0][0]-delIndicesDict[np.where(temp==sortedFH[0])[0][0]]  
                for val in fhMatch: fhTBD.remove(val) 
            
            if current not in delIndices:   # Ieterate through branch once again, if a core was not flagged for deletion
                fhMatch = np.where(temp2['FirstHaloInFOFGroupOffset']==current)[0]  # Find all instances of current core in tree (firstHaloInFOFGroup, nextHaloInFOFGroup, Descendent)
                if fhMatch.size: temp['FirstHaloInFOFGroupOffset'][fhMatch] = current-len(delIndices)   # Replace all instances of current core's index by index - # cores to be deleted
                 
                nhMatch = np.where(temp2['NextHaloInFOFGroupOffset']==current)[0] 
                dMatch = np.where(temp2['DescendentOffset']==current)[0] 
                if nhMatch.size: temp['NextHaloInFOFGroupOffset'][nhMatch] = current-len(delIndices) 
                if dMatch.size: temp['DescendentOffset'][dMatch] = current-len(delIndices) 
                 
                if temp[current]['FirstProgenitorOffset'] != -1: # Do the same for firstProgenitor  
                    temp[current]['FirstProgenitorOffset'] = current-len(delIndices)+1  # But instead replace by index - # cores to be deleted + 1 to abide by indexing conventions 
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
    print(len(outputTrees)-1, treenum)
    # outputTrees.append(trees[treenum].copy())

    # Drawing
    #trees[treenum]=temp
    #treegraph, clusters, nodes, node_names, first_progenitors, next_progenitors = bctDevV2.drawforest(trees, treenum, xname='red')                                   

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

