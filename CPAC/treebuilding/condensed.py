# Written by Eve K.

import sys
import os
import h5py
import re
import psutil
import glob
import math
import numpy as np
#import numba
from time import time
import struct
from struct import pack
from struct import unpack, calcsize
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import pydot
#import pydotplus as pydot

#cc_template = '09_03_2019.AQ.{}.corepropertiesextend.hdf5'
cc_template = 'm000p-{}.corepropertiesextend.hdf5'
#cc_template = '{}.coreproperties_complete.hdf5'
binfile_template = 'trees_099.{}'
#coredir = '../CoreCatalogs'
coredir = '/scratch/cpac/kmaamari/outputLJDS'
#coredir = '/home/kmaamari/lgalaxies/h5generation/genericio/python/output'
#coredir = '../CoreCatalogsLJR/LJ_CoreCatalog_004_reduced'
#treedir = '../CoreTrees'
#treedir = '../CoreTrees/test'
#treedir = '../CoreTrees/new_snaps'
#treedir = '../CoreTrees/relative'
#treedir = '../CoreTrees/fof_group'
treedir = '/scratch/cpac/kmaamari/CoreTreesLJDS'
mkey = 'm_evolved_0.9_0.001'
#mkey = 'm_evolved_0.8_0.02'
#mkey = 'm_evolved_0.9_0.005'
outfile_template = re.sub('propertiesextend', 'trees', cc_template)

coretag = 'core_tag'
foftag = 'fof_halo_tag'
coremass = 'coremass'
infall_fof_mass = 'infall_fof_halo_mass'
infall_tree_node_mass = 'infall_tree_node_mass'
first_snap = 499
last_snap = 43
#last_snap = 475
first_row = 0
last_row = 99 
#last_row = 2 
Descendent = 'Descendent'
FirstProgenitor = 'FirstProgenitor'
NextProgenitor = 'NextProgenitor'
FirstHaloInFOFGroup = 'FirstHaloInFOFGroup'
NextHaloInFOFGroup = 'NextHaloInFOFGroup'
MostBoundID_Coretag = 'MostBoundID->Coretag'
SubhaloIndex = 'SubhaloIndex'
ParentHaloTag = 'ParentHaloTag'
Offset = 'Offset'
DescendentOffset = Descendent + Offset
FirstProgenitorOffset = FirstProgenitor + Offset
NextProgenitorOffset = NextProgenitor + Offset
FirstHaloInFOFGroupOffset = FirstHaloInFOFGroup + Offset
NextHaloInFOFGroupOffset = NextHaloInFOFGroup + Offset
Len = 'Len'
Zero = 'Zero'
SnapNum = 'SnapNum'

# properties for storing in dict or matrices
core_pointers = [Descendent, DescendentOffset, FirstProgenitor, NextProgenitor,
                 FirstProgenitorOffset, NextProgenitorOffset]
#core_pointers = [Descendent, DescendentOffset, FirstProgenitor, FirstProgenitorOffset]
sibling_pointers = [FirstHaloInFOFGroup, NextHaloInFOFGroup, FirstHaloInFOFGroupOffset, NextHaloInFOFGroupOffset]
#sibling_pointers = []
core_properties_float = {'Pos_x':'x', 'Pos_y':'y', 'Pos_z':'z',
                         'Vel_x':'vx', 'Vel_y':'vy', 'Vel_z':'vz',
                         'VelDisp':'vel_disp',
                         'Vmax': 'infall_fof_halo_max_cir_vel',
                         'M_Crit200': infall_fof_mass,
                         'M_Mean200': Zero,
                         'M_TopHat': Zero,
                         'Spin_x': 'infall_sod_halo_angmom_x',
                         'Spin_y': 'infall_sod_halo_angmom_y',
                         'Spin_z': 'infall_sod_halo_angmom_z',
                         'SubHalfMass': Zero,
                        }

core_properties_int = {MostBoundID_Coretag:'core_tag',
                       'FileNr': Zero,
                       SubhaloIndex: Zero,
                      }

derived_properties_int = [SnapNum, Len, ParentHaloTag]
#derived_properties_int = [Len]

# vector properties for storing in matrix
integer_properties = core_pointers + sibling_pointers + derived_properties_int + list(core_properties_int.keys())
float_properties = list(core_properties_float.keys())
#float_properties = ['M_Crit200']

#serial properties for storing in dict
core_properties = dict(zip(float_properties, [core_properties_float[f] for f in float_properties]))
core_properties.update(core_properties_int)
derived_properties = derived_properties_int

no_int = -999
no_float = -999.
particle_mass = 1.15e9
particle_mass_AQ = 1.15e9 #M_sun/h
particle_mass_LJ = 1.15e9 #M_sun/h
particle_mass_MT = 1.15e9 #M_sun/h

# header file contains Ntrees, totNHalos, TreeNHalos
header_format = "<{}i"

struct_keys = [DescendentOffset, FirstProgenitorOffset, NextProgenitorOffset,
               FirstHaloInFOFGroupOffset, NextHaloInFOFGroupOffset,
               Len, 'M_Mean200', 'M_Crit200', 'M_TopHat',
               'Pos_x', 'Pos_y', 'Pos_z', 'Vel_x', 'Vel_y', 'Vel_z',
               'VelDisp', 'Vmax', 'Spin_x', 'Spin_y', 'Spin_z',
               MostBoundID_Coretag, SnapNum, 'FileNr',
               SubhaloIndex, 'SubHalfMass',
              ]
struct_format = "<iiiiiiffffffffffffffqiiif"
"""
#test
struct_keys = [DescendentOffset, FirstProgenitorOffset, NextProgenitorOffset,
               FirstHaloInFOFGroupOffset, NextHaloInFOFGroupOffset,
               Len, 'M_Crit200', MostBoundID_Coretag,
#               SnapNum, SubhaloIndex,
              ] 

#struct_format = "<iiiiiifqii"
struct_format = "<iiiiiifq"
"""

# cc = build_core_trees.get_core_snapshot( '../CoreCatalogs', snapshot)
'''
def get_core_snapshot(coredir, snapshot, template=cc_template):
    #print("\n1: get_core_snapshot")
    fn = os.path.join(coredir, cc_template.format(snapshot))
    data= {}
    if os.path.exists(fn):
        h5 = h5py.File(fn, 'r')
        coredata = h5['coredata']
        keys = [k for k in list(coredata.keys()) if 'm_evolved' not in k]    

        for k in keys + [mkey]:
            print(k)
            if k=='infall_fof_halo_mass': # change to mevolved
                data[k] = coredata[k][()]
                print("Len before:",len(data[k]))
                delIndices = np.where(data[k]<100*particle_mass)[0]
                #temp = np.delete(data[k].copy(),delIndices)
                #data[k]=temp
                print("Len after:",len(data[k]))
        for k in keys + [mkey]:
            if k=='infall_fof_halo_mass':
                continue
            else:
                if k=='infall_tree_node_mass':
                    data[k] = coredata[k][()]
                else:
                    data[k] = coredata[k][()]
                temp = np.delete(data[k].copy(),delIndices)
            data[k]=temp
    else:
        print('{} not found'.format(fn))
    
    #print("\t1: get_core_snapshot")
    return delIndices,data

def get_core_snapshot(coredir, snapshot, template=cc_template):
    fn = os.path.join(coredir, cc_template.format(snapshot))
    data= {}
    print(fn)
    if os.path.exists(fn):
        h5 = h5py.File(fn, 'r')
        coredata = h5['coredata']
        keys = [k for k in list(coredata.keys()) if 'm_evolved' not in k]

        for k in keys + [mkey]:
            if k=='infall_fof_halo_mass': # change to mevolved
                data[k] = coredata[k][()]
                print("Len before:",len(data[k]))
                delIndices = np.where(data[k]<100)[0]
                for i in delIndices:
                    data[k][i]=-999
                print("Len after:",len(data[k]))
        for k in keys + [mkey]:
            if k=='infall_fof_halo_mass':
                continue
            else:
                if k=='infall_tree_node_mass':
                    data[k] = coredata[k][()]
                else:
                    data[k] = coredata[k][()]
                for i in delIndices:
                    data[k][i]=-999

    else:
        print('{} not found'.format(fn))

    return delIndices, data
'''
def get_core_snapshot(coredir, snapshot, template=cc_template):
    fn = os.path.join(coredir, cc_template.format(snapshot))
    data= {}
    if os.path.exists(fn):
        h5 = h5py.File(fn, 'r')
        coredata = h5['coredata']
        keys = [k for k in list(coredata.keys()) if 'm_evolved' not in k]
        for k in keys + [mkey]:
            data[k] = coredata[k][()]
    else:
        print('{} not found'.format(fn))

    return data

def add_coremass_column(corecat):
    #print("\n2: add_coremass_column")
    mask = (corecat['central']==1) # get centrals
    central_mass = corecat[infall_tree_node_mass][mask] # get fof mass
    corecat[coremass] = corecat[mkey] # get evolved masses
    corecat[coremass][mask] = central_mass
    #print("\t2: add_coremass_column")
    return corecat

def add_offsets(coretrees, col, offset, count):
    #add offset to coretree column if element != -1
    for p in [Descendent, FirstProgenitor]: #, NextProgenitor]:
        mask = coretrees[p][:, col] > -1
        coretrees[p + Offset][:, col][mask] += offset

    return coretrees

def clean(corecat, sorted_coretags):
    #print("\n3: clean")
    mask = np.in1d(corecat[coretag], sorted_coretags, assume_unique=True)
    print('Truncating core catlog to {}/{} entries'.format(np.count_nonzero(mask),
                                                           len(corecat[coretag])))
    for p in list(corecat.keys()):
        # print("(Before mask) Property =",p,":",corecat[p])
        corecat[p] = corecat[p][mask]
        # print("(After mask) Property =",p,":",corecat[p])

    #print("\t3: clean")
    return corecat

def add_snapshot_to_trees(coretrees, corecat, current_snap, snap_index=0,
                          print_int=100, coretags_fofm=None, argsorted_fofm=None,
                          sorted_indices=None,
                          ncore_min=0, ncore_max=None, vector=False): 
    #print("\n4: add_snapshot_to_trees")
    assert coretags_fofm is not None, "No sibling-ordered tags supplied"
    assert sorted_indices is not None, "No sibling-ordered sorted indices supplied"
    #mass_order = np.flip(corecat[coremass].argsort())   #args in descending order

    # sort indices in this snapshot in decreasing foftag-mass order
    if current_snap < first_snap:
        sorted_indices_this = np.flip(np.lexsort(((corecat[coremass], corecat[foftag]))))
        coretags_fofm_this = corecat[coretag][sorted_indices_this]
    else:
        sorted_indices_this = sorted_indices
        coretags_fofm_this = coretags_fofm
    foftags_fofm_this = corecat[foftag][sorted_indices_this]  
        
    # get unique parents; index points to first occurence of tag; inverse reconstructs array
    parent_halos, index, inverse, counts = np.unique(foftags_fofm_this,
                                                     return_index=True, return_inverse=True,
                                                     return_counts=True)
    # initialize times
    stimes =  []
    atimes = []
    if vector:
        assert argsorted_fofm is not None, "No sibling-ordered argsorted indices supplied"
        # get required reordering of coretags in current snapshot for insertion into matrix
        # cannot assume cortetags will be in the same order as in last snapshot
        if current_snap < first_snap:
            # argsorted_fofm = coretags_fofm.argsort()
            insertions = np.searchsorted(coretags_fofm[argsorted_fofm], coretags_fofm_this)
            locations_this = argsorted_fofm[insertions]
        else:
            locations_this = np.arange(0, len(coretags_fofm))
        # get siblings as vectors in sorted-indices order 
        stime1 = time()
        first_siblings, next_siblings = get_sibling_vectors(coretags_fofm_this, index, inverse)
        stimes.append((time()-stime1)/60.)

        ptime1 = time()
        # loop over property matrices in coretree dict
        for p, v in coretrees.items():     
            print(p)
            #define/reorder property
            atime1 = time()
            prop_values = get_ordered_property(p, corecat, sorted_indices_this,
                                               snap_index, current_snap,
                                               first_siblings, next_siblings, foftags_fofm_this)
            atimes.append((time()-atime1)/60.)
            #fill row of property matrix with values for selected entries
            if p is 'Descendent' or p is 'FirstProgenitor' or p is 'NextProgenitor' or p is 'NextHaloInFOFGroup' or p is 'FirstHaloInFOFGroup':
                print(p)
                for i,location in enumerate(locations_this):
                    if snap_index==0: v[snap_index][location] = prop_values[i]
                    elif v[snap_index-1][location]!=-1: v[snap_index][location] = prop_values[i]
            else: v[snap_index][locations_this] = prop_values
        # fix first progenitors
        coretrees[FirstProgenitor] = fix_firstprogenitor_vector(coretrees[FirstProgenitor], snap_index)
        coretrees[FirstProgenitorOffset] = fix_firstprogenitor_vector(coretrees[FirstProgenitorOffset], snap_index)
    
    #print("\t4: add_snapshot_to_trees")
    return coretrees, atimes, stimes

def orphanize():
    descendent=-1


def get_sibling_vectors(coretags, index, inverse):
    """
    coretags: array sorted in decreasing tag, mass order
    index, inverse, counts: output of np.unique
    """
    #print("\n5: get_sibling_vectors")
    # find most massive core corresponding to parent halos
    first_sibling_coretags = coretags[index]
    # use inverse to reconstruct first-sibling array
    first_siblings = first_sibling_coretags[inverse]

    # initialize array and pointers for next siblings
    ncores = len(coretags)
    nparents = len(index)
    next_siblings = np.array([-1]*ncores) #null pointer
    last_sibling_loc = np.zeros(nparents - 1).astype(int)
    # shift coretag array by 1 argument to the left; last element is left as -1
    next_siblings[0:ncores-1] = coretags[1:ncores]
    # shift index locations by -1 to point to last sibling in each group
    last_sibling_loc[0:nparents-1] = index[1:nparents] - 1
    next_siblings[last_sibling_loc] = -1
    
    #print("\t5: get_sibling_vectors")
    return first_siblings, next_siblings


# fix first progenitors
def fix_firstprogenitor_vector(first_progenitor, row):
    #print("\n6: fix_firstprogenitor_vector")
    if row > 0:
        mask_this = (first_progenitor[row] == no_int) #no entry in this snap
        mask_prev = (first_progenitor[row-1] != no_int) #entry in previous snap
        mask = mask_this & mask_prev
        np.set_printoptions(edgeitems=50)
        #print("1",mask_this,"2",mask_prev,"3",mask)
        first_progenitor[row-1][mask] = -1 #overwrite
    
    #print("\t6: fix_firstprogenitor_vector")
    return first_progenitor

def get_ordered_property(p, corecat, sorted_indices_this, row, current_snap,
                         first_siblings, next_siblings, parent_tags):
    #print("\n7: get_ordered_property")
    
    ncores = len(corecat[coretag]) # = len(sorted_indices_this)
    orphans = np.where(np.round(corecat[coremass][sorted_indices_this]/particle_mass).astype(int)<100)[0]
    if Descendent in p:
        prop_values = np.array([row - 1]*ncores) #row = row of matrix array
        for i in tqdm(orphans):
            prop_values[i]=-999 
    elif FirstProgenitor in p:
        prop_values = np.array([row + 1]*ncores) if row != last_row else np.array([-1]*ncores)
        for i in tqdm(orphans):
            prop_values[i]=-999
    elif NextProgenitor in p:
        prop_values = np.array([-1]*ncores)
        for i in tqdm(orphans):
            prop_values[i]=-999
    elif FirstHaloInFOFGroup in p:
        prop_values = first_siblings
        for i in tqdm(orphans):
            prop_values[i]=-999
    elif NextHaloInFOFGroup in p:
        prop_values = next_siblings
        for i in tqdm(orphans):
            prop_values[i]=-999
    elif SnapNum in p:
        #prop_values = np.array([current_snap]*ncores)
        prop_values = np.array([last_row - row]*ncores) #L-galaxies needs consecutive integers
    elif 'FileNr' in p:
        prop_values = np.array([0]*ncores)
    elif ParentHaloTag in p:
        prop_values = parent_tags
    elif Len in p:
        prop_values = np.round(corecat[coremass][sorted_indices_this]/particle_mass).astype(int) #truncate
        mask = (prop_values == 0)
        prop_values[mask] = 1  #ensure non-zero values
    elif p in list(core_properties_int.keys()):
        if Zero in core_properties_int[p]:
            prop_values = np.zeros(ncores).astype(int)
        else:
            prop_values = corecat[core_properties_int[p]][sorted_indices_this] # reorder into sorted coretag order
    elif p in list(core_properties_float.keys()):
        if Zero in core_properties_float[p]:  #zero values
            prop_values = np.zeros(ncores)
        else:
            prop_values = corecat[core_properties_float[p]][sorted_indices_this] # reorder into sorted coretag order
            #if 'M_Crit' in p:
                #mask = (corecat['central'][sorted_indices_this]==0) # select non-centrals
                #prop_values[mask] = 0.   # set staellite masses to 0.
            if 'Spin' in p:
                prop_values /= corecat[infall_fof_mass][sorted_indices_this]
            if 'M_' in p:
                mask = (corecat['central'][sorted_indices_this]==0) # select non-centrals
                prop_values[mask] = 0.   # set staellite masses to 0.
                prop_values /= particle_mass
    else:
        print('Unknown property {}'.format(p))
    #print("Property =",p,":",prop_values)
    #print("\t7: get_ordered_property")
    return prop_values

def get_parent_boundary(foftags, loc, ncores, name='input', upper=True):
    #print("\n8: get_parent_boundary")
    # assumes foftags are in fofm order
    oldloc = loc
    foftag_this = foftags[loc-1] # parent of preceding index
    args =  np.where(foftags == foftag_this)[0] #locations of this parent
    if upper:
        loc = min(np.max(args) + 1, ncores) # maximum argument +1
    else:
        loc = min(np.min(args) - 1, 0)
    if loc != oldloc:
        print('Resetting {} to {} to line up with parent boundary'.format(name, loc))
    
    #print("\t8: get_parent_boundary")
    return loc

def write_outfile(outfile, coretrees, cores_to_write, vector=False, start=None, end=None,
                  column_counts=None):
    #print("\n9: write_outfile")
    """
    coretrees: trees to output
    cores_to_write: list of sorted cores
    vector: True if using vectorized code (requires start and end arguments)
    start: starting column of coretree matrix
    end: ending column of coretree matrix (= parent_boundary+1) 
    """
    wtime = time()
    if os.path.exists(outfile):
        os.remove(outfile)
        print("Removed old file {}".format(outfile))
    if len(cores_to_write) == 0:
        return
    hdfFile = h5py.File(outfile, 'w')
    hdfFile.create_group('treeData')
    hdfFile['treeData']['Ntrees'] = len(cores_to_write)
    # use Descendent to get lengths, masks etc. 
    if vector:
        if column_counts is None:
            print('Need column counts')
            return
        hdfFile['treeData']['TreeNHalos'] = column_counts
    hdfFile['treeData']['coretags'] = cores_to_write
    hdfFile['treeData']['totNHalos'] = np.sum(hdfFile['treeData']['TreeNHalos'])
    
    print('Number of core trees = {}, Total halos = {}'.format(hdfFile['treeData']['Ntrees'][()],
                                                          hdfFile['treeData']['totNHalos'][()]))
    
    # concatenate trees
    tgroup = hdfFile.create_group('trees')
    if vector:
        mtime = time()
        properties_list = list(coretrees.keys())
        mask = (coretrees[Descendent][:, start:end] != -999)  # mask for non-existent elements
        flat_mask = mask.flatten(order='F')
    else:
        properties_list = list(coretrees[cores_to_write[0]].keys())

    ptime = time()
    for p in properties_list:
        if vector:
            prop_array = coretrees[p][:, start:end].flatten(order='F') 
            tgroup[p] = prop_array[flat_mask]
            
        assert len(tgroup[p][()]) == hdfFile['treeData']['totNHalos'][()], \
            'Length mismatch for {}'.format(p,  len(tgroup[p][()]), hdfFile['treeData']['totNHalos'][()])

    #print("\t9: write_outfile")
    hdfFile.close()
                                      
def get_column_counts(coretree_matrix_transpose):
    #print("\n10: get_column_counts")
    #vectorize
    counts = np.zeros(len(coretree_matrix_transpose)).astype(int)
    for n, c in enumerate(coretree_matrix_transpose):
        mask = (c != -999)
        counts[n] = np.count_nonzero(mask)
    #print("\t10: get_column_counts")
    return counts

def replace_sibling_addresses(coretrees, col, column_counts, start, first_column_in_tree):
    for p in [FirstHaloInFOFGroup, NextHaloInFOFGroup]:
        cindx = col - start # correct for offset in column counts (evaluated for each file)
        mask = (coretrees[p][:, col][0:column_counts[cindx]] > 0)
        #print(col, column_counts[cindx], coretrees[p][:, col])
        if np.count_nonzero(mask) > 0:  # check if any entries are coretags
            unique_sibs, inverse = np.unique(coretrees[p][:, col][0: column_counts[cindx]], return_inverse=True)
            rows = np.arange(0, column_counts[cindx]) # vector of row indices
            #print(col, unique_sibs, inverse, rows)
            for ns, sib in enumerate(unique_sibs):
                if sib > 0:   #skip -1 entries
                    loc = np.where(coretrees[MostBoundID_Coretag][0]==sib)[0][0]
                    assert loc >= first_column_in_tree, "Sibling location is outside of tree"
                    offset = np.sum(column_counts[first_column_in_tree:loc])
                    sib_mask = (inverse==ns) & mask
                    coretrees[p + Offset][:, col][0:column_counts[cindx]][sib_mask] = offset + rows[sib_mask]

    return coretrees

def write_binary(outfile, coretrees, cores_to_write, foftags_to_write,
                 vector=True, start=None, end=None, column_counts=None):
    #print("\n11: write_binary")
    """
    coretrees: trees to output
    cores_to_write: list of sorted cores
    vector: True if using vectorized code (requires start and end arguments)
    start: starting column of coretree matrix
    end: ending column of coretree matrix (= parent_boundary+1) 
    column_counts: counts of cores in each tree to be written
    """

    wtime = time()
    rtimes = []
    if column_counts is None:
        print('Need column counts')
        return
    if len(cores_to_write) == 0:
        print('No cores to write')
        return
        
    with open(outfile, "wb") as fh:
        # L-galaxies requires trees organized by parent groups
        fof_groups, index, fof_counts = np.unique(foftags_to_write, return_index=True, return_counts=True)
        fof_group_order = index.argsort()   #reorder parent halos back to original order
        #get number of halos in each fof_group
        halos_per_tree = get_halo_counts(fof_counts[fof_group_order], column_counts)
        
        # write global information for file
        Ntrees = len(fof_groups)
        totNHalos = int(np.sum(column_counts))        
        values = [Ntrees, totNHalos] + halos_per_tree.tolist()
        header_this = header_format.format(Ntrees + 2)
        fh.write(pack(header_this, *values))
        
        if vector:
            assert np.array_equal(foftags_to_write, coretrees[ParentHaloTag][0, start:end]), "ParentHaloTags not in fofm order"
            validate_trees(column_counts, coretrees, start, end)
            #TODO validate parent tags and FirstHaloTags (order will be same but tags will not)
            
            # write struct for each tree (ie forest) (cores in same parent halo at z=0)
            column = start
            for fof_count in fof_counts[fof_group_order]:
                offset = 0   #reset offset for each tree
                for col, count in zip(np.arange(column, column + fof_count),
                                      column_counts[column-start:column-start+fof_count]):  #column counts are restarted for each file
                    coretrees = add_offsets(coretrees, col, offset, count)
                    rtime = time()
                    coretrees = replace_sibling_addresses(coretrees, col, column_counts,
                                                          start, column)
                    rtimes.append(time() - rtime)
                    for row in np.arange(count):
                        values = [coretrees[p].T[col][row] for p in struct_keys]
                        fh.write(pack(struct_format, *values))
                    offset = offset + count   #offset addresses by # of rows in column
                    
                column += fof_count    
    #print("\t11: write_binary")
    return

def get_halo_counts(counts, column_counts):
    #print("\n12: get_halo_counts")
    # sum up halos in each fof group
    
    halos_per_tree = np.zeros(len(counts)).astype(int)
    offset = 0
    for n, count in enumerate(counts):
        halos_per_tree[n] = np.sum(column_counts[offset:offset+count])
        offset += count
    
    #print("\t12: get_halo_counts")
    return halos_per_tree

def validate_trees(counts, coretrees, start, end):
    #print("\n13: validate_trees")
    # check coretag is unique and equal to counts
    for count, col in zip(counts, np.arange(start, end)):
        unique_tags = np.unique(coretrees[MostBoundID_Coretag].T[col][0:count])
        assert len(unique_tags)==1, 'Mutiple core tags {} found in column {}'.format(unique_tags,
                                                                                     col)
    #print("\t13: validate_trees")
    return

def main(argv):
    #print("\n0: main")
    #setup parameters
    if type(argv) == list:
        argv = dict(zip(np.arange(0, len(argv)), argv))
    print('Inputs:',argv)
    vector = True if argv.get(1, 'vector')=='vector' else False
    nfiles = int(argv.get(2, 3))   #number of files to write
    print_int = int(argv.get(3, 1000)) #for serial code
    ncore_min = int(argv.get(4, 0))    #for serial code
    ncore_max = int(argv.get(5, 67760)) #for serial code (number of cores)/100
    Nfiles = int(argv.get(6, 10000)) #total number of files for vector code
    name = argv.get(0, 'test')
    print('Outputs written to {}'.format(treedir))
    
    process = psutil.Process(os.getpid())
    
    corefiles = glob.glob(coredir+'/*')
    snapshots = sorted([int(os.path.basename(f).split('-')[1].split('.')[0]) for f in corefiles], reverse=True)  
    coretrees = {}
    delIndicesArr = []

    for n, s in enumerate(snapshots): #process in descending order
        print('\nProcessing snapshot {}'.format(s))
        stime = time() 
        corecat = get_core_snapshot(coredir, int(s))
        if corecat:
            corecat = add_coremass_column(corecat)
            if n == 0:
                sorted_coretags = np.sort(corecat[coretag]) #save for cleaning earlier snaps

                # sort corecat first by foftag and then coremass so siblings are grouped together
                indices_fofm = np.flip(np.lexsort(((corecat[coremass], corecat[foftag]))))
                coretags_fofm = corecat[coretag][indices_fofm] #order of coretags in coretrees
                foftags_fofm = corecat[foftag][indices_fofm]
                Ncores = len(coretags_fofm)
                argsorted_coretags_fofm = coretags_fofm.argsort()

                if vector: # initialize dict of matrices
                    print('Setting up ordered arrays and tree matrices')
                    ctime = time()
                    for i,p in enumerate(tqdm(float_properties)):
                        coretrees[p] = np.array([no_float]*Ncores*len(snapshots)).reshape(len(snapshots), Ncores)
                    for i,p in enumerate(tqdm(integer_properties)):
                        coretrees[p] = np.array([no_int]*Ncores*len(snapshots)).reshape(len(snapshots), Ncores)
                corecat = clean(corecat, sorted_coretags)
            else:
                # clean corecat to remove cores that disappeared before step 499
                corecat = clean(corecat, sorted_coretags)
            
            coretrees, atimes, stimes = add_snapshot_to_trees(coretrees, corecat, int(s), snap_index=n,
                                                              coretags_fofm=coretags_fofm,
                                                              argsorted_fofm=argsorted_coretags_fofm,
                                                              sorted_indices=indices_fofm,
                                                              print_int=print_int, ncore_max=ncore_max,
                                                              ncore_min=ncore_min, vector=vector)
                    
        mem = "Memory usage =  {0:.2f} GB"
        print(mem.format(process.memory_info().rss/1.e9))
    
    del corecat

    # TODO merging
    
    # output core trees
    numcores = Ncores if vector else len(list(coretrees.keys()))
    numfiles = Nfiles if vector else nfiles
    stride = int(np.ceil(numcores/numfiles))
    start = 0 if vector else ncore_min
    mode = '.vector' if vector else '.serial'
    print('Writing subsets of {} cores in {} files with stride {}'.format(numcores, nfiles, stride))
    
    for n in range(nfiles):
        fn = os.path.join(treedir, outfile_template.format(n)+mode)
        fn_bin = os.path.join(treedir, binfile_template.format(n)+mode)
        end = int(min(start + stride, numcores))
        # check that end is on a parent boundary
        end = get_parent_boundary(foftags_fofm, end, Ncores, name='end')
          
        cores_to_write = coretags_fofm[start:end]
        foftags_to_write = foftags_fofm[start:end]
        if len(cores_to_write) > 0:
            if vector:
                column_counts = get_column_counts(coretrees[Descendent].T[start:end]) #get counts for unfilled elements
            # write binary & hdf5
            write_outfile(fn, coretrees, cores_to_write, vector=vector, start=start, end=end,
                          column_counts=column_counts)
            write_binary(fn_bin, coretrees, cores_to_write, foftags_to_write,
                         vector=vector, start=start, end=end, column_counts=column_counts)
        start = end
    
    # testtree = build_core_trees.main({0:'test', 1:'serial', 2:1, 3:400, 4:246, 5:247})
    return coretrees
