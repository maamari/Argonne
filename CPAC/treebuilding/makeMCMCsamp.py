import bctDevV2
from tqdm import tqdm

lgalVersion = 'lgal'
#outfileLJ='/scratch/cpac/kmaamari/CoreTreesLJDS/trees_099.0.vector'
outfileLJ='/home/kmaamari/lgalaxies/{}/MergerTrees/LJDS/treedata/trees_099.0'.format(lgalVersion)
#outfileLJ='massiveTrees_099.0'
LJtrees=bctDevV2.read_binary(outfileLJ)

IDindex=20
fileNRindex=22
snapIndex=21
numTrees=len(LJtrees)

for i in tqdm(range(100,-1,-1)):
    snapnum=i
    filename='/home/kmaamari/lgalaxies/{}/MCMC/Samples/cut_optimalsample_allz_nh_100{}.dat'.format(lgalVersion,snapnum)
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
            #weight=LJtrees[treenum][0]['Len']
            #weight=100
            length=LJtrees[treenum][0]['Len']
            if length<=20:
                weight=238.28321251704585
            elif length>20 and length<=50:
                weight=493.1376626463074
            elif length>50 and length<=100:
                weight=1339.3577537518154
            elif length>100 and length<=200:
                weight=2464.3705463182896
            elif length>200 and length<=500:
                weight=3958.9792511328405
            elif length>500 and length<=1000:
                weight=11416.781292984868
            elif length>1000 and length<=5000:
                weight=15132.178669097539
            elif length>5000:
                weight=146902.65486725664
            print("\t",ID,"\t",treenum,"\t",fileNR,"\t",weight,file=sampleDAT)
            count+=1
        except:
            #print("Failed!")
            continue
    sampleDAT.close()

print(outfileLJ, filename)
