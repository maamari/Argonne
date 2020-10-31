import bctDevV2 
import numpy as np 
from tqdm import tqdm 
from struct import pack 
import sys 

def validateBranch(tree, node): 
    if node[0]==-1: 
        if node[3]==0: 
            return 
        else: 
            sys.exit("Bad branch") 
            return 
    node = tree[node[0]] 
    validateBranch(tree, node) 
 
outfile='trees_099.0.vector' 
trees=bctDevV2.read_binary(outfile) 
 
debug=False 
branchingTrees = [] 
for i in range(len(trees)): 
    #if ((trees[i]['NextHaloInFOFGroupOffset']>-1).any()):
    if (trees[i][0]['Len']>0) and ((trees[i]['NextHaloInFOFGroupOffset']>-1).any()): 
        branchingTrees.append(i) 
 
print("Num trees =", len(branchingTrees)) 
massiveTrees = [] 
for treenum in tqdm(branchingTrees): 
    temp = trees[treenum] 
   
    if debug: print(temp)
    for index in range(len(temp)-1): 
        if abs(temp[index+1]['FirstHaloInFOFGroupOffset']-temp[index]['FirstHaloInFOFGroupOffset']) > 1: 
            temp[index+1]['DescendentOffset'] = temp[index]['FirstHaloInFOFGroupOffset'] 
    temp = temp[temp['FirstHaloInFOFGroupOffset']==np.arange(len(temp))] 
    temp['NextHaloInFOFGroupOffset'] = -1 
 
    if debug: print("\n",temp) 
    for index in range(len(temp)): 
        if temp[index]['FirstHaloInFOFGroupOffset']!=index: 
            temp['DescendentOffset'] = np.where(temp['DescendentOffset']==temp[index]['FirstHaloInFOFGroupOffset'],index,temp['DescendentOffset']) 
            temp['FirstProgenitorOffset'] = np.where(temp['FirstProgenitorOffset']==temp[index]['FirstHaloInFOFGroupOffset'],index,temp['FirstProgenitorOffset']) 
    temp['FirstHaloInFOFGroupOffset'] = np.arange(len(temp)) 
 
    if debug: print("\n",temp) 
    for index in range(len(temp)-1): 
        if temp[index+1]['DescendentOffset']!=-1 and not len(np.where(temp[index+1]['DescendentOffset']==temp['FirstHaloInFOFGroupOffset'])[0]): 
            temp[index+1]['DescendentOffset']=temp[index]['FirstHaloInFOFGroupOffset'] 
 
    if debug: print("\n",temp) 
    for index in range(len(temp)): 
        if temp[index]['DescendentOffset'] != -1: 
            descendent = temp[index]['DescendentOffset'] 
            if temp[descendent]['FirstProgenitorOffset']!=index: 
                temp[descendent]['FirstProgenitorOffset'] = index 
 
    if debug: print("\n",temp) 
    for i in range(len(temp)): 
        progs = sorted(temp[temp['DescendentOffset']==i], key=lambda tree:tree['Len'], reverse=True) 
        
        if len(progs)>1: 
            firstProg = temp[i]['FirstProgenitorOffset']     
            progs_ = [] 

            for prog in progs: 
                if prog[3]==firstProg: continue 
                else: progs_.append(prog) 
            
            temp[firstProg]['NextProgenitorOffset']=progs_[0][3] 
            for i in range(len(progs_)-1): 
                progs[i]['NextProgenitorOffset'] = progs_[i+1]['FirstHaloInFOFGroupOffset'] 
                temp[[prog[3] for prog in progs_]] = progs_ 
    
    for i in range(len(temp)):
        temp[i]['M_Crit200']*=(2.7e9/1.e10)

    if debug: print("\n",temp) 
    massiveTrees.append(temp)

for index in tqdm(range(len(massiveTrees))):
    tree = massiveTrees[index]
    ends = np.where(tree['FirstProgenitorOffset']==-1)[0]

    for end in ends:
        validateBranch(tree, tree[end])

outfile='massiveTrees_099.0'
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

