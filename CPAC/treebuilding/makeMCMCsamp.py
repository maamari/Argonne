import build_core_trees
from tqdm import tqdm

#outfileLJ='/scratch/cpac/kmaamari/CoreTreesLJDS/trees_099.0.vector'
outfileLJ='/home/kmaamari/lgalaxies/lgal/MergerTrees/LJDS/treedata/trees_099.0'
LJtrees=build_core_trees.read_binary(outfileLJ)

IDindex=20
fileNRindex=22
snapIndex=21
numTrees=1000

for i in reversed(range(100)):
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
    print("\t",count,file=sampleDAT)
#    print('\t',92,file=sampleDAT)
    for treenum in range(numTrees):
        try:
            ID = LJtrees[treenum][99-snapnum][IDindex]
            fileNR = LJtrees[treenum][99-snapnum][fileNRindex]
            weight=1000
            #print("\t",ID,"\t",treenum,"\t",fileNR,"\t",weight)
            print("\t",ID,"\t",treenum,"\t",fileNR,"\t",weight,file=sampleDAT)
            count+=1
        except: 
            continue
    sampleDAT.close()
    print(i,"\n")
