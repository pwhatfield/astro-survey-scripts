
# Create a range of relevant subsets of the data

import numpy as np
import pylab as pb
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import random
from astropy.cosmology import Planck13
pb.ion()
import matplotlib.pyplot as plt
from astropy.table import Table

#Intro params
axes_numbers_size=22
axes_labels_size=24
marker_size=0.5



upper_redshift_axis=4
upper_stellar_mass=12
lower_stellar_mass=8

version=1

script_directory='...'
data_directory='...'

mag_lim=np.genfromtxt(script_directory+'config_files/mag_lim_XMM.txt', delimiter=',')

# The redshift grid for the pdfs
zgrid_value=np.concatenate((np.arange(0,6,0.04),np.arange(6,9.1,0.1)))

###### IMPORT DATA

### COSMOS
## Import data
#print('Loading COSMOS Data')
#data_COS = Table.read(data_directory+'COS_Master_V4.fits')
#print('Finished')
#
#print('Loading COSMOS PDFs')
#pdfs_COS = np.genfromtxt(data_directory+'PDF/formatted_pdfs/COS_HB_pdfs',delimiter=',')
#print('Finished')
## Get number of sources
#num_sources_COS=np.shape(data_COS)[0]



## XMM
# Import data
print('Loading XMM-LSS Data')
data_XMM = Table.read(data_directory+'XMM_Master_V4.fits')
print('Finished')

print('Loading XMM-LSS PDFs')
pdfs_XMM = np.genfromtxt(data_directory+'PDF/formatted_pdfs/XMM_HB_pdfs',delimiter=',')
print('Finished')

# Get number of sources
num_sources_XMM=np.shape(data_XMM)[0]


# Load stellar masses, photometric redshifts and K band fluxes
#
## COS
#mass_galaxies_COS = data_COS['MASS_SE']
#redshifts_galaxies_COS = data_COS['Z_phot_best']
#K_flux_COS=data_COS['flux_Ks']
#K_mag_COS=data_COS['Ks']
#J_mag_COS=data_COS['J']
#g_mag_COS=data_COS['HSC-G_ssp']
#i_mag_COS=data_COS['HSC-I_ssp']
#RA_galaxies_COS = data_COS['RA']
#DEC_galaxies_COS = data_COS['DEC']

# XMM
mass_galaxies_XMM = data_XMM['MASS_SE']
redshifts_galaxies_XMM = data_XMM['Z_phot_best']
K_flux_XMM=data_XMM['flux_Ks']
K_mag_XMM=data_XMM['Ks']
J_mag_XMM=data_XMM['J']
g_mag_XMM=data_XMM['HSC-G']
i_mag_XMM=data_XMM['HSC-I']
RA_galaxies_XMM = data_XMM['RA']
DEC_galaxies_XMM = data_XMM['DEC']

###### 


# Defining functions for stellar locus

print('Define Stellar Locus Functions')

def f_stellar(x):
    
    if x<=0.4:
        return -0.58
    elif 0.4<x<=1.9:
        return -0.88 + 0.82*x -0.21*x*x
    else:
        return -0.08
    

def f_stellar_array(x):
    
    vals=np.zeros(len(x))
    
    for j in range(len(x)):
        vals[j]=f_stellar(x[j])
    
    return vals
        
    
    



def stellar_locus_exclusion(J,K,g,i):
    
    stellar_indicies=np.ones(len(J_mag_XMM), dtype=bool)
    
    for j in range(len(J)):
        if J[j]-K[j]<0.12 + f_stellar(g[j]-i[j]):
            stellar_indicies[j]=False
        else:
            stellar_indicies[j]=True
    
    return stellar_indicies



print('Create the Subsets')

z_bin_edges=np.genfromtxt(script_directory+'config_files/z_bins.txt', delimiter=',')
mass_bin_edges=np.genfromtxt(script_directory+'config_files/mass_bins.txt', delimiter=',')
bins_used=np.genfromtxt(script_directory+'config_files/bins_used.txt', delimiter=',')

for i in range(len(z_bin_edges)-1):
    for j in range(len(mass_bin_edges)-1):
        if bins_used[len(mass_bin_edges)-j-2,i]==1:
            
            # Save XMM data
            
            # Indicies
            bin_indicies=np.logical_and(np.logical_and(mass_bin_edges[j+1]>mass_galaxies_XMM,mass_galaxies_XMM>mass_bin_edges[j]),np.logical_and(z_bin_edges[i+1]>redshifts_galaxies_XMM,redshifts_galaxies_XMM>z_bin_edges[i]))
            stellar_indicies=stellar_locus_exclusion(J_mag_XMM,K_mag_XMM,g_mag_XMM,i_mag_XMM)
            mag_indicies=np.logical_and(np.logical_and(K_mag_XMM<mag_lim,g_mag_XMM<80),g_mag_XMM>-80)
            
            indicies=np.logical_and(np.logical_and(bin_indicies,stellar_indicies),mag_indicies)
            
            #Extract RA and DEC
            data=np.zeros([np.sum(indicies),2])
            data[:,0]=RA_galaxies_XMM[indicies]
            data[:,1]=DEC_galaxies_XMM[indicies]
            
            # Save
            np.savetxt(script_directory+'output_data/data_coordinates_XMM_redshift_'+str(i)+'_mass_'+str(j)+'.csv',data, delimiter=',')
            
            # Save XMM bin pdfs
            pdf=np.zeros(len(zgrid_value))
            pdf=np.sum(pdfs_XMM[indicies,:],axis=0)
            pdf=pdf/np.trapz(pdf, zgrid_value, axis=0)
            # Save
            np.savetxt(script_directory+'output_data/sample_redshift_pdfs_XMM_redshift_'+str(i)+'_mass_'+str(j)+'.csv',pdf, delimiter=',')
            


print('Finished')

g_i_values_plot_galaxies_XMM=g_mag_XMM[np.logical_and(mag_indicies,stellar_indicies)]-i_mag_XMM[np.logical_and(mag_indicies,stellar_indicies)]
J_K_values_plot_galaxies_XMM=J_mag_XMM[np.logical_and(mag_indicies,stellar_indicies)]-K_mag_XMM[np.logical_and(mag_indicies,stellar_indicies)]

np.savetxt(script_directory+'output_data/g_i_values_plot_galaxies_XMM.csv',g_i_values_plot_galaxies_XMM, delimiter=',')
np.savetxt(script_directory+'output_data/J_K_values_plot_galaxies_XMM.csv',J_K_values_plot_galaxies_XMM, delimiter=',')


g_i_values_plot_stellar_XMM=g_mag_XMM[np.logical_and(mag_indicies,1-stellar_indicies)]-i_mag_XMM[np.logical_and(mag_indicies,1-stellar_indicies)]
J_K_values_plot_stellar_XMM=J_mag_XMM[np.logical_and(mag_indicies,1-stellar_indicies)]-K_mag_XMM[np.logical_and(mag_indicies,1-stellar_indicies)]

np.savetxt(script_directory+'output_data/g_i_values_plot_stellar_XMM.csv',g_i_values_plot_stellar_XMM, delimiter=',')
np.savetxt(script_directory+'output_data/J_K_values_plot_stellar_XMM.csv',J_K_values_plot_stellar_XMM, delimiter=',')

# Save good indicies
np.savetxt(script_directory+'output_data/good_indicies_XMM.csv',np.logical_and(mag_indicies,stellar_indicies), delimiter=',')

