import bctDevV2
from tqdm import tqdm

#outfileLJ='/scratch/cpac/kmaamari/CoreTreesLJDS/trees_099.0.vector'
outfileLJ='/home/kmaamari/lgalaxies/lgal/MergerTrees/LJDS/treedata/trees_099.0'
#outfileLJ='massiveTrees_099.0'
LJtrees=bctDevV2.read_binary(outfileLJ)

IDindex=20
fileNRindex=22
snapIndex=21
numTrees=len(LJtrees)

for i in tqdm(range(100,-1,-1)):
    snapnum=i
    filename='../../lgal/MCMC/Samples/cut_optimalsample_allz_nh_100{}.dat'.format(snapnum)
    sampleDAT=open(filename, 'w')
    count=0
    for treenum in range(numTrees):
        try:
            ID = LJtrees[treenum][99-snapnum][IDindex]
            fileNR = LJtrees[treenum][99-snapnum][fileNRindex]
            count+=1
        except:
            continue
    for treenum in range(numTrees):
        try:
            ID = LJtrees[treenum][99-snapnum][IDindex]
            fileNR = LJtrees[treenum][99-snapnum][fileNRindex]
            weight=100
            count+=1
        except: 
            continue
    sampleDAT.close()
