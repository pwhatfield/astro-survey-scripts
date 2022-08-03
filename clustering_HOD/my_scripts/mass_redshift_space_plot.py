
# Plot RA and DEC and also mass-redshift space

import numpy as np
import pylab as pb
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import random
from astropy.cosmology import Planck13
pb.ion()
import matplotlib.pyplot as plt
from astropy.table import Table

from matplotlib import rc
rc('text', usetex=True)
rc('font', family='sans-serif')





axes_numbers_size=22
axes_labels_size=24
marker_size=0.5

mag_lim=np.genfromtxt(script_directory+'config_files/mag_lim_XMM.txt', delimiter=',')

upper_redshift_axis=4
upper_stellar_mass=12
lower_stellar_mass=8

version=1

script_directory='...'
data_directory='...'


###### IMPORT DATA

good_indicies=np.genfromtxt(script_directory+'output_data/good_indicies_XMM.csv', delimiter=',')
good_indicies=np.array(good_indicies, dtype=bool)

## COSMOS
# Import data
print('Loading COSMOS Data')
data_COS = Table.read(data_directory+'COS_Master_V4.fits')
print('Finished')

# Get number of sources
num_sources_COS=np.shape(data_COS)[0]

RA_galaxies_COS = data_COS['RA']
DEC_galaxies_COS = data_COS['DEC']


## XMM
# Import data
print('Loading XMM-LSS Data')
data_XMM = Table.read(data_directory+'XMM_Master_V4.fits')
print('Finished')

# Get number of sources
num_sources_XMM=np.shape(data_XMM)[0]

RA_galaxies_XMM = data_XMM['RA'][good_indicies]
DEC_galaxies_XMM = data_XMM['DEC'][good_indicies]


###### PLOT RA AND DEC

print('Start Plotting RA and DEC')

RA_range=4
DEC_range=2.5

RA_min_COS=148
RA_max_COS=RA_min_COS+RA_range

DEC_min_COS=1.0
DEC_max_COS=DEC_min_COS+DEC_range

RA_min_XMM=33.5
RA_max_XMM=RA_min_XMM+RA_range

DEC_min_XMM=-6
DEC_max_XMM=DEC_min_XMM+DEC_range


# Plot COSMOS
fig=plt.figure(figsize=(12, 6),dpi=96)
#plt.figure()
ax1 = fig.add_subplot(1,2,1)
plt.scatter(RA_galaxies_COS,DEC_galaxies_COS,color=(1,0,0),s=marker_size,label="Galaxies")

pb.xlim(RA_min_COS,RA_max_COS)
pb.ylim(DEC_min_COS,DEC_max_COS)

pb.xlabel(r'RA', fontsize=axes_numbers_size)
pb.ylabel(r'DEC', fontsize=axes_numbers_size)

plt.gca().invert_xaxis()

plt.text(RA_min_COS+3.25, DEC_min_COS+2.2, r'COSMOS Field', fontsize=25)




# Plot XMM
ax2 = fig.add_subplot(1,2,2)
plt.scatter(RA_galaxies_XMM,DEC_galaxies_XMM,color=(1,0,0),s=marker_size)

###
data=np.genfromtxt(script_directory+'output_data/random_coordinates.csv', delimiter=',')
plt.scatter(data[:,0],data[:,1],color=(0,0,1),s=marker_size,label="Galaxies")
###

pb.xlim(RA_min_XMM,RA_max_XMM)
pb.ylim(DEC_min_XMM,DEC_max_XMM)

pb.xlabel(r'RA', fontsize=axes_numbers_size)
pb.ylabel(r'DEC', fontsize=axes_numbers_size)

plt.gca().invert_xaxis()

plt.text(RA_min_XMM+3.25, DEC_min_XMM+2.2, r'XMM-LSS Field', fontsize=25)


pb.savefig(script_directory+"output_images/RA_DEC_VIDEO_COS_XMM.pdf")

#plt.close()

print('Finished')





########### Plot Redshifts and Stellar Masses


# Load stellar masses, photometric redshifts and K band fluxes

# COS
mass_galaxies_COS = data_XMM['MASS_SE']
redshifts_galaxies_COS = data_XMM['Z_phot_best']
K_flux_COS=data_XMM['flux_Ks']
K_mag_COS=data_XMM['Ks']

# XMM
mass_galaxies_XMM = data_XMM['MASS_SE'][good_indicies]
redshifts_galaxies_XMM = data_XMM['Z_phot_best'][good_indicies]
K_flux_XMM=data_XMM['flux_Ks'][good_indicies]
K_mag_XMM=data_XMM['Ks'][good_indicies]

#######

print('Start Plotting Stellar Mass and Redshift')

fig = plt.figure(figsize=(7.5, 4.5),dpi=200)
ax1 = fig.add_subplot(1,1,1)
plt.gcf().subplots_adjust(top=0.8,bottom=0.25,left=0.2)

mass_lim = 0*K_mag_COS

# This is finding the mass it would have been if it had been at the detection threshold
for i in range(len(mass_galaxies_COS)):
    mass_lim[i]=10**(mass_galaxies_COS[i]+0.4*(K_mag_COS[i]-mag_lim))


# Introducing some scatter for plotting purposes
for i in range(len(redshifts_galaxies_COS)):
    redshifts_galaxies_COS[i]=redshifts_galaxies_COS[i]+random.uniform(-0.05,0.05)


#
p1=plt.scatter(redshifts_galaxies_COS,10**mass_galaxies_COS,color=(0,0,1),s=0.0012,label='_nolegend_') # Plot data
p2=plt.scatter(redshifts_galaxies_COS,mass_lim,color=(1,0,0),s=0.0012,label='_nolegend_')

print('Finished')

print('Calculate Completeness Threshold')

# This bit calculates the completeness threshold
numz_bins=20
lower_z=0
upper_z=upper_redshift_axis+0.5
bin_centres=np.arange(lower_z+0.5*(upper_z-lower_z)/numz_bins,upper_z+0.5*(upper_z-lower_z)/numz_bins,(upper_z-lower_z)/numz_bins)
mass_completeness_lim=np.zeros(shape=(numz_bins))
completeness=90

for j in range(numz_bins):
    indices = [i for i, x in enumerate(redshifts_galaxies_COS) if (upper_z-lower_z)*j/numz_bins<x <(upper_z-lower_z)*(j+1)/numz_bins]
    mass_completeness_lim[j]=np.percentile(mass_lim[indices],completeness)
    
completeness_function = spline(bin_centres,np.log10(mass_completeness_lim),k=3)

pb.plot(np.linspace(lower_z,upper_z,800),10**(completeness_function(np.linspace(lower_z,upper_z,800))+0.1),'g', label="$M_{lim}(K_{s})$")

pb.xscale("linear")
pb.yscale("log")
pb.xlim(0,upper_redshift_axis)
pb.xlabel(r'Redshift - z', fontsize=axes_numbers_size)
pb.ylim(10**lower_stellar_mass,10**upper_stellar_mass)
pb.ylabel(r'Stellar Mass - $M_{\star} / M_{\odot}$', fontsize=axes_numbers_size)


ax1.tick_params('both', length=8, width=1, which='major', labelsize=axes_numbers_size)
ax1.tick_params('both', length=5, width=0.5, which='minor')
ax1.set_yticks([10**7,10**8,10**9,10**10,10**11,10**12])

print('Finished')


print('Plot HAT2016 Bin Boundaries')

# HAT2016 Boundaries
z0=0.5;
z1=0.75;
z2=1;
z3=1.25;
z4=1.7;
z5=2.3;

# This is after upwards conversion
m1=10**9.35;
m2=10**9.6;
m3=10**9.85;
m4=10**10.1;
m5=10**10.35;
m6=10**10.6;
m7=10**10.85;

pb.plot([z0, z0], [m1,10**12], color='c', linestyle='-', linewidth=4)
pb.plot([z1, z1], [m1,10**12], color='c', linestyle='-', linewidth=4)
pb.plot([z2, z2], [m2,10**12], color='c', linestyle='-', linewidth=4)
pb.plot([z3, z3], [m3,10**12], color='c', linestyle='-', linewidth=4)
pb.plot([z4, z4], [m4,10**12], color='c', linestyle='-', linewidth=4)
#pb.plot([z5, z5], [m5,10**12], color='r', linestyle='-', linewidth=4)

pb.plot([z0, z1], [m1,m1], color='c', linestyle='-', linewidth=4)
pb.plot([z0, z2], [m2,m2], color='c', linestyle='-', linewidth=4)
pb.plot([z0, z3], [m3,m3], color='c', linestyle='-', linewidth=4)
pb.plot([z0, z4], [m4,m4], color='c', linestyle='-', linewidth=4)
pb.plot([z0, z4], [m5,m5], color='c', linestyle='-', linewidth=4)
pb.plot([z0, z4], [m6,m6], color='c', linestyle='-', linewidth=4)
pb.plot([z1, z3], [m7,m7], color='c', linestyle='-', linewidth=4)

print('Finished')


print('Plot New Bin Boundaries')

z_bin_edges=np.genfromtxt(script_directory+'config_files/z_bins.txt', delimiter=',')
mass_bin_edges=np.genfromtxt(script_directory+'config_files/mass_bins.txt', delimiter=',')
bins_used=np.genfromtxt(script_directory+'config_files/bins_used.txt', delimiter=',')

for i in range(len(z_bin_edges)-1):
    for j in range(len(mass_bin_edges)-1):
        if bins_used[len(mass_bin_edges)-j-2,i]==1:
            pb.plot([z_bin_edges[i], z_bin_edges[i]], [10**mass_bin_edges[j], 10**mass_bin_edges[j+1]], color='k', linestyle='-', linewidth=2)
            pb.plot([z_bin_edges[i+1], z_bin_edges[i+1]], [10**mass_bin_edges[j], 10**mass_bin_edges[j+1]], color='k', linestyle='-', linewidth=2)
            pb.plot([z_bin_edges[i], z_bin_edges[i+1]], [10**mass_bin_edges[j], 10**mass_bin_edges[j]], color='k', linestyle='-', linewidth=2)
            pb.plot([z_bin_edges[i], z_bin_edges[i+1]], [10**mass_bin_edges[j+1], 10**mass_bin_edges[j+1]], color='k', linestyle='-', linewidth=2)


print('Finished')




##### Do some cosmology for the top x-axis

print('Sort out the labels and lookback time')

handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[::-1], labels[::-1],prop={'size':16},loc=4,ncol=1,numpoints=1)

ax1.set_xlim(0, upper_redshift_axis)
ax2 = ax1.twiny()
ax1Xs = ax1.get_xticks()





z = np.linspace(0,upper_redshift_axis,100) # Array of redshifts
t_of_z=Planck13.lookback_time(z).value # Lookback times
z_of_t = spline(t_of_z,z,k=3) # call spline with inverted axes


ax2.set_xlim(0, upper_redshift_axis)
ax2.set_ylim(10**lower_stellar_mass, 10**upper_stellar_mass)
ax2.set_xticks([z_of_t(1), z_of_t(2), z_of_t(3),z_of_t(4),z_of_t(5),z_of_t(6),z_of_t(7),z_of_t(8),z_of_t(9),z_of_t(10),z_of_t(11),z_of_t(12)])
ax2.set_xticklabels(['1','2','3','4','5','6','7','8','9','10','11','12'])

pb.xlabel(r'Lookback Time (Gyr)', fontsize=axes_numbers_size)
ax2.tick_params('both', length=8, width=1, which='major', labelsize=axes_numbers_size/2)
ax2.tick_params('both', length=5, width=0.5, which='minor')

print('Finished')

pb.savefig(script_directory+"output_images/mass_redshift_VIDEO_COS_XMM.pdf")



# Save data

indicies=np.logical_and(np.logical_and(11>mass_galaxies_XMM,mass_galaxies_XMM>10),np.logical_and(1.5>redshifts_galaxies_XMM,redshifts_galaxies_XMM>1))

data=np.zeros([np.sum(indicies),2])
data[:,0]=RA_galaxies_XMM[indicies]
data[:,1]=DEC_galaxies_XMM[indicies]
np.savetxt(script_directory+'output_data/data_coordinates.csv',data, delimiter=',')

