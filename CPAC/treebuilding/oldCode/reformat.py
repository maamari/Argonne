import bctDevV2 
import numpy as np 
from tqdm import tqdm 
from struct import pack

#outfile='/scratch/cpac/kmaamari/CoreTreesDev/trees_099.0.vector' 
outfile='trees_099.0.vector'
trees=bctDevV2.read_binary(outfile) 
  
branchingTrees = [] 
for i in range(len(trees)): 
    if((trees[i]['NextHaloInFOFGroupOffset']>-1).any()): 
        branchingTrees.append(i) 

for treenum in tqdm(branchingTrees):
    temp = trees[treenum]
    
    for index in range(len(temp)-1):
        if temp[index+1]['FirstHaloInFOFGroupOffset']-temp[index]['FirstHaloInFOFGroupOffset'] > 1: 
        #if temp[index+1]['DescendentOffset']-temp[index]['DescendentOffset'] > 1:
            temp[index+1]['DescendentOffset']=temp[index]['FirstHaloInFOFGroupOffset']
    temp = temp[temp['FirstHaloInFOFGroupOffset']==np.arange(len(temp))]
    temp['NextHaloInFOFGroupOffset']=-1

    for index in range(len(temp)):
        if temp[index]['FirstHaloInFOFGroupOffset']!=index:
            temp['DescendentOffset'] = np.where(temp['DescendentOffset']==temp[index]['FirstHaloInFOFGroupOffset'],index,temp['DescendentOffset'])
            temp['FirstProgenitorOffset'] = np.where(temp['FirstProgenitorOffset']==temp[index]['FirstHaloInFOFGroupOffset'],index,temp['FirstProgenitorOffset'])
    temp['FirstHaloInFOFGroupOffset']=np.arange(len(temp))

    for index in range(len(temp)):
        if temp[index]['DescendentOffset'] != -1:
            descendent = temp[index]['DescendentOffset']
            try: 
                if temp[descendent]['FirstProgenitorOffset']!=index: 
                    temp[descendent]['FirstProgenitorOffset'] = index
            except: 
                #print(treenum,"failed")
                continue

    for i in range(len(temp)):
        progs = sorted(temp[temp['DescendentOffset']==i], key=lambda tree:tree['Len'], reverse=True)
        if len(progs)>1:
            for i in range(len(progs)-1):
                progs[i]['NextProgenitorOffset']=progs[i+1]['FirstHaloInFOFGroupOffset']
                temp[[prog[3] for prog in progs]]=progs
  
    trees[treenum] = temp
    #treegraph, clusters, nodes, node_names, first_progenitors, next_progenitors = bctDevV2.drawforest(trees, treenum, compressed=False, xname='red')

outfile='trees_099.0'
struct_keys = ['DescendentOffset', 'FirstProgenitorOffset', 'NextProgenitorOffset','FirstHaloInFOFGroupOffset', 'NextHaloInFOFGroupOffset','Len', 'M_Mean200', 'M_Crit200', 'M_TopHat','Pos_x', 'Pos_y', 'Pos_z', 'Vel_x', 'Vel_y', 'Vel_z','VelDisp', 'Vmax', 'Spin_x', 'Spin_y', 'Spin_z','MostBoundID->Coretag', 'SnapNum', 'FileNr','SubhaloIndex', 'SubHalfMass']
struct_format = "<iiiiiiffffffffffffffqiiif"
header_format = "<{}i"

with open(outfile, "wb") as fh:
    # write global information for file
    Ntrees = len(trees)
    totNHalos = 0
    for i in range(len(trees)):
        totNHalos += len(trees[i])        
    values = [Ntrees, totNHalos]
    for i in range(len(trees)):
        values.append(len(trees[i]))
    header_this = header_format.format(Ntrees + 2)
    fh.write(pack(header_this, *values))
    print('Wrote header for {} fof trees with {} halos to {}'.format(Ntrees, totNHalos, outfile))

    for i in range(len(trees)):
        for j in range(len(trees[i])):
            values = list(trees[i][j])
            fh.write(pack(struct_format, *values))
