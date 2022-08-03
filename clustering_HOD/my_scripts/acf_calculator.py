# Load packages
from os.path import dirname, abspath, join as pjoin
import numpy as np
import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.io import read_catalog
from Corrfunc.utils import convert_3d_counts_to_cf
import matplotlib.pyplot as plt
import pylab as pb
pb.ion()

axes_numbers_size=22
axes_labels_size=24
marker_size=0.5

upper_redshift_axis=4
upper_stellar_mass=12
lower_stellar_mass=8

# Angular scales probed
ang_scale_min_log=-3.0
ang_scale_max_log=-0.5
nbins = 10

number_bootstrap=10

# Location of data
script_directory='...'
mask_directory='...'

def angular_corr_gal_calculator(RA_data,DEC_data,RA_random, DEC_random, theta_min_log, theta_max_log, theta_num,RR_counts,nthreads):
    
    N = len(RA_data)
    rand_N = len(RA_random)

    # Setup the bins
    bins = np.logspace(theta_min_log, theta_max_log, theta_num + 1) # note the +1 to nbins
    bin_centres= np.logspace(theta_min_log +0.5*(theta_max_log-theta_min_log)/theta_num, theta_max_log -0.5*(theta_max_log-theta_min_log)/theta_num, theta_num)
    
    print('Calculate DD')
    # Auto pair counts in DD
    autocorr=1
    DD_counts = DDtheta_mocks(autocorr, nthreads, bins,RA_data, DEC_data)
    
    print('Calculate DR')
    # Cross pair counts in DR
    autocorr=0
    DR_counts = DDtheta_mocks(autocorr, nthreads, bins,RA_data, DEC_data,RA2=RA_random, DEC2=DEC_random)
    
#    print('Calculate RR')
#    # Auto pairs counts in RR
#    autocorr=1
#    RR_counts = DDtheta_mocks(autocorr, nthreads, bins,RA_random, DEC_random)
    
    print('Calculate wtheta')
    # All the pair counts are done, get the angular correlation function
    wtheta = convert_3d_counts_to_cf(N, N, rand_N, rand_N, DD_counts, DR_counts,DR_counts, RR_counts)
    
    return bin_centres,wtheta

def angular_corr_gal_calculator_with_uncertainty(RA_data,DEC_data,RA_random, DEC_random, theta_min_log, theta_max_log, theta_num,n_bootstrap):
    
    N = len(RA_data)
    rand_N = len(RA_random)
    
    # Number of threads to use
    nthreads = 1
    # Setup the bins
    bins = np.logspace(theta_min_log, theta_max_log, theta_num + 1) # note the +1 to nbins
    bin_centres= np.logspace(theta_min_log +0.5*(theta_max_log-theta_min_log)/theta_num, theta_max_log -0.5*(theta_max_log-theta_min_log)/theta_num, theta_num)

    print('Calculate RR')
    # Auto pairs counts in RR
    autocorr=1
    RR_counts = DDtheta_mocks(autocorr, nthreads, bins,RA_random, DEC_random)
    
    wtheta_array=np.zeros([n_bootstrap,nbins])
    
    for i in range(n_bootstrap):
        
        RA_resampled=RA[np.random.randint(0,N,size=N)]
        DEC_resampled=DEC[np.random.randint(0,N,size=N)]
        
        # Calculate ACF
        bin_centres,wtheta_array[i,:]=angular_corr_gal_calculator(RA_resampled,DEC_resampled,rand_RA, rand_DEC, ang_scale_min_log, ang_scale_max_log, nbins,RR_counts,nthreads)
    
    wtheta=np.percentile(wtheta_array,50,axis=0)
    wtheta_upper=np.percentile(wtheta_array,84,axis=0)
    err_wtheta=wtheta_upper-wtheta
    
    return bin_centres,wtheta,err_wtheta




# Load the RA and DEC of the random points
data=np.genfromtxt(script_directory+'output_data/random_coordinates.csv', delimiter=',')
rand_RA, rand_DEC = data.T





###### 


print('Calculate ACF for all bins')

z_bin_edges=np.genfromtxt(script_directory+'config_files/z_bins.txt', delimiter=',')
mass_bin_edges=np.genfromtxt(script_directory+'config_files/mass_bins.txt', delimiter=',')
bins_used=np.genfromtxt(script_directory+'config_files/bins_used.txt', delimiter=',')



for i in range(len(z_bin_edges)-1):
    for j in range(len(mass_bin_edges)-1):
        if bins_used[len(mass_bin_edges)-j-2,i]==1:
            
            print('zbin'+str(i)+' massbin'+str(j))
            
            # Load RA and DEC
            data=np.genfromtxt(script_directory+'output_data/data_coordinates_XMM_redshift_'+str(i)+'_mass_'+str(j)+'.csv',delimiter=',')
            RA=data[:,0]
            DEC=data[:,1]
            
            print('Loaded')
            
            # Calculate ACF
            bin_centres,wtheta,err_wtheta=angular_corr_gal_calculator_with_uncertainty(RA,DEC,rand_RA, rand_DEC, ang_scale_min_log, ang_scale_max_log, nbins,number_bootstrap)
            
            # Save correlation function
            data=np.zeros([len(bin_centres),3])
            data[:,0]=bin_centres
            data[:,1]=wtheta
            data[:,2]=err_wtheta
            np.savetxt(script_directory+'output_data/acf_XMM_redshift_'+str(i)+'_mass_'+str(j)+'.csv',data, delimiter=',')
            


print('Finished')




