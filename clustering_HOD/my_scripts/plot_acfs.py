# Plot the ACFs from acf_calculator.py

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pylab as pb
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import random
from astropy.cosmology import Planck13
#pb.ion()
import matplotlib.pyplot as plt
from astropy.table import Table

from matplotlib import rc

rc('text', usetex=True)
rc('font', family='sans-serif')


script_directory='...'
data_directory='...'

# Angular scales probed
ang_scale_min_log=-3.0
ang_scale_max_log=-0.5
nbins = 10

axis_name_label_size=10
major_tick_size=5

n_rows=3
n_columns=3

##
fig = plt.figure(figsize=(7.5, 4.5),dpi=200)
#ax1 = fig.add_subplot(1,1,1)
plt.gcf().subplots_adjust(top=0.9,bottom=0.1,left=0.1,right=0.9,wspace = 0.2,hspace = 0.3)




###### 


print('Calculate ACF for all bins')

z_bin_edges=np.genfromtxt(script_directory+'config_files/z_bins.txt', delimiter=',')
mass_bin_edges=np.genfromtxt(script_directory+'config_files/mass_bins.txt', delimiter=',')
bins_used=np.genfromtxt(script_directory+'config_files/bins_used.txt', delimiter=',')




for i in range(len(z_bin_edges)-1):

    ax1 = fig.add_subplot(n_rows,n_columns,i+1)

    for j in range(len(mass_bin_edges)-1):
        
        
        if bins_used[len(mass_bin_edges)-j-2,i]==1:
            
            
            
            print('zbin'+str(i)+' massbin'+str(j))
            
            # Load ACFs
            print('Load ACFs')
            data=np.genfromtxt(script_directory+'output_data/acf_XMM_redshift_'+str(i)+'_mass_'+str(j)+'.csv', delimiter=',')
            bin_centres=data[:,0]
            wtheta=data[:,1]
            err_wtheta=data[:,2]
            
            print('Loaded')
            
            #pb.plot(bin_centres,wtheta,color=(1.0-((len(mass_bin_edges)-j-1)/(len(mass_bin_edges)-1)),0,((len(mass_bin_edges)-j-1)/(len(mass_bin_edges)-1))),linewidth=1.0,linestyle='-') # Plot data
            plt.errorbar(bin_centres,wtheta, err_wtheta,color=(1.0-((len(mass_bin_edges)-j-1)/(len(mass_bin_edges)-1)),0,((len(mass_bin_edges)-j-1)/(len(mass_bin_edges)-1))),label='$'+str(mass_bin_edges[j])+'  < \log_{10} (M_{\star}/M_{\odot}) <'+str(mass_bin_edges[j+1])+'$')
            
            print('Plotted')
    
    z_range_label=r'$'+str(z_bin_edges[i])+'<z<'+str(z_bin_edges[i+1])+'$'
    plt.text(10**(-1.8), 10**(0), z_range_label, fontsize=axis_name_label_size*0.6)
    pb.xscale("log")
    pb.yscale("log")
    pb.xlim(10**ang_scale_min_log,10**ang_scale_max_log)
    pb.ylim(10**(-1.5),10**(0.5))
    
    # Tick size
    ax1.tick_params(axis='both', which='major', labelsize=major_tick_size)
    #ax1.tick_params(axis='both', which='minor', labelsize=8)
    
    # Add x label just for bottom row
    if i>=(n_columns-1)*n_rows:
        pb.xlabel(r'$\theta \mathrm{ (^{\circ})}$', fontsize=axis_name_label_size)
    
    # Add y label just for left column
    if 0==np.mod(i,n_columns):
        pb.ylabel(r'$w(\theta)$', fontsize=axis_name_label_size)
    
    # Add legend
    if i==0:
        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels,fontsize=6, loc='upper center',ncol=3)



print('Finished')



pb.savefig(script_directory+"output_images/acfs_for_all_redshifts.png")
