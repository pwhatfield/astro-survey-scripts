
# Plot the pdfs and stellar locus from create_subsets.py


#import matplotlib
#matplotlib.use('Agg')
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


g_i_values_plot_galaxies_XMM=np.genfromtxt(script_directory+'output_data/g_i_values_plot_galaxies_XMM.csv', delimiter=',')
J_K_values_plot_galaxies_XMM=np.genfromtxt(script_directory+'output_data/J_K_values_plot_galaxies_XMM.csv', delimiter=',')

g_i_values_plot_stellar_XMM=np.genfromtxt(script_directory+'output_data/g_i_values_plot_stellar_XMM.csv', delimiter=',')
J_K_values_plot_stellar_XMM=np.genfromtxt(script_directory+'output_data/J_K_values_plot_stellar_XMM.csv',delimiter=',')

z_bin_edges=np.genfromtxt(script_directory+'config_files/z_bins.txt', delimiter=',')
mass_bin_edges=np.genfromtxt(script_directory+'config_files/mass_bins.txt', delimiter=',')
bins_used=np.genfromtxt(script_directory+'config_files/bins_used.txt', delimiter=',')


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
        
    
    



print('Start Plotting Colour Space')

fig = plt.figure(figsize=(7.5, 4.5),dpi=200)
ax1 = fig.add_subplot(1,1,1)
plt.gcf().subplots_adjust(top=0.8,bottom=0.25,left=0.2)

plt.scatter(g_i_values_plot_galaxies_XMM,J_K_values_plot_galaxies_XMM,color=(1,0,0),s=0.00001,label='_nolegend_')
plt.scatter(g_i_values_plot_stellar_XMM,J_K_values_plot_stellar_XMM,color=(0,0,1),s=0.00001,label='_nolegend_')

x_values=np.arange(-1,3.5,0.1)
y_values=f_stellar_array(x_values)
pb.plot(x_values, y_values, color='k', linestyle='--', linewidth=2,label='Stellar Locus')


plt.scatter(1000,1000,color=(1,0,0),s=0.5,label='Galaxies')
plt.scatter(1000,1000,color=(0,0,1),s=0.5,label='Stars')


pb.xscale("linear")
pb.yscale("linear")
pb.xlim(0,3.5)
pb.xlabel(r'g-i', fontsize=axes_numbers_size)
pb.ylim(-1,2)
pb.ylabel(r'J-Ks', fontsize=axes_numbers_size)

handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[::-1], labels[::-1],prop={'size':9})

pb.savefig(script_directory+"output_images/stellar_locus.pdf")




print('Start Plotting PDFs')

y_scale_max=5
y_scale_min=-1

fig = plt.figure(figsize=(7.5, 4.5),dpi=200)
ax1 = fig.add_subplot(1,1,1)
plt.gcf().subplots_adjust(top=0.8,bottom=0.25,left=0.2)

for i in range(len(z_bin_edges)-1):
    pdf=np.genfromtxt(script_directory+'output_data/sample_redshift_pdfs_XMM_redshift_'+str(i)+'_mass_'+str(len(mass_bin_edges)-2)+'.csv',delimiter=',')
    pdf=pdf/np.trapz(pdf, zgrid_value, axis=0)
    pb.plot(zgrid_value,pdf,color=(1.0-((len(z_bin_edges)-i-1)/(len(z_bin_edges)-1)),0,((len(z_bin_edges)-i-1)/(len(z_bin_edges)-1))))
    pb.plot([z_bin_edges[i], z_bin_edges[i]], [y_scale_min,y_scale_max], color='k', linestyle='--', linewidth=2)
    pb.plot([z_bin_edges[i+1], z_bin_edges[i+1]], [y_scale_min,y_scale_max], color='k', linestyle='--', linewidth=2)



pb.xscale("linear")
pb.yscale("linear")
pb.xlim(0,5)
pb.xlabel(r'$z$', fontsize=axes_numbers_size)
pb.ylim(y_scale_min,y_scale_max)
pb.ylabel(r'$P(z)$', fontsize=axes_numbers_size)


print('Finished')

# Save good indicies
pb.savefig(script_directory+"output_images/redshift_distributions.pdf")
