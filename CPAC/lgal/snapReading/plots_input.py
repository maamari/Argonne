
# coding: utf-8

# In[ ]:


#PLOT OPTIONS
opt_stellar_mass_vs_halo_mass=1
opt_stellar_mass_function=1
opt_metals_vs_stellarmass=0
opt_BHBM=1
opt_SFRF=1
opt_gas_fraction=0
opt_HI_MF=1
opt_sfr_vs_stellar_mass=1
opt_ur_vs_r=0
opt_UVJ_colour=0
opt_redfraction_color_cut=0
    
opt_plot_MCMC_sample=1

#COSMOLOGIES & DARK MATTER SIMS
WMAP1=0
PLANCK=1
CATERPILLAR_PLANCK=0

DirName_MR = '/home/kmaamari/lgalaxies/lgal/output/'

Datadir = '/home/kmaamari/lgalaxies/lgal/AuxCode/Python/data/'
MCMCdir = '/Users/BrunoHenriques/Desktop/OneDrive/Workspace/GitHub_PR_Hen15/MCMC/'
MCMCSampledir = '/home/kmaamari/lgalaxies/lgal/output/'

prefix_this_model='Last Journey'
file_this_model='ThisWork'

do_previous_model1=1
file_previous_model1=Datadir+'Guo2013a_m05'
prefix_previous_model1='Guo2013a - WMAP7'
linestyle_previous_model1=':'

#do_previous_model2=1
#file_previous_model2=Datadir+'Henriques2013a'
#prefix_previous_model2='Henriques2013a - WMAP7'
#linestyle_previous_model2='--'

do_previous_model2=1
file_previous_model2=Datadir+'/Henriques2015a'
prefix_previous_model2='Henriques2015 - PLANCK1'
linestyle_previous_model2='--'


slope_red_fraction=[0.075,0.275, 0.3, 0.32,0.38]
offset_red_fraction=[1.85,1.213, 1.18,0.99,0.79]
minimum_y_red_fraction=[0.0,1.3,1.3,1.3,1.3]

RedshiftsToRead = [True,True,True,True,True,True,False]

CatalogType='snap'
#CatalogType='tree'

Hubble_h_WMAP1 = 0.732
Hubble_h_WMAP7 = 0.704


if WMAP1: 
    FullRedshiftList=[0.00,0.41,0.99,2.07,3.06,3.87] 
    FullSnapshotList=[63,50,41,32,27,24]  
    BoxSize_MR    = 500. #full MR 
    BoxSize_MRII  = 100. #full MRII      
    Hubble_h      = 0.73
    Omega_M       = 0.25 
    Omega_Lambda  = 0.75
    MaxTreeFiles  = 512
    

if PLANCK: 
    FullRedshiftList=[0.00,0.10,0.24,0.30,0.40,0.62,1.01] 
    FullSnapshotList=[99,95,90,88,85,79,70]  
    BoxSize_MR    = 256.* 0.960558 #full MR 
    BoxSize_MRII  = 256.* 0.960558 #full MRII      
    Hubble_h      = 0.71
    Omega_M       = 0.265
    Omega_Lambda  = 0.735
    MaxTreeFiles  = 999
    
if CATERPILLAR_PLANCK:
    FullRedshiftList=[0.00,0.10,0.40,1.00,2.01,2.98,3.96]
    FullSnapshotList=[255,231,181,125,81,60,47]
    BoxSize_MR    = 100.
    Hubble_h      = 0.673
    Omega_M       = 0.315
    Omega_Lambda  = 0.683
    MaxTreeFiles  = 8

    



