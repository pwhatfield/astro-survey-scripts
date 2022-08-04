# Do the hierarchical combination as per https://arxiv.org/pdf/1709.09183.pdf (see also https://arxiv.org/pdf/1712.04476.pdf with some script from Ken Duncan

# Import stuff
import numpy as np
from scipy.stats import norm

from scipy.interpolate import InterpolatedUnivariateSpline as spline

print('Starting')

# Open files
ML_pdfs=np.loadtxt('.../COS_ml_pdfs',delimiter=',')
template_pdfs=np.loadtxt('.../COS_template_pdfs',delimiter=',')

# Blank array for combined pdfs
combined_pdfs=0*template_pdfs

# Number of galaxies
n_gal=995049#1674689


# Define HB parameters
beta_value=1
zgrid_value=np.concatenate((np.arange(0,1,0.01),np.arange(1,9.04,0.04)))
fbad_min_value=0
fbad_max_value=0.05
nbad_value=30

# Number of redshift points
n_z=np.shape(zgrid_value)
n_z=n_z[0]

# High resolution for spline
zgrid_high_resolution=(np.arange(0,9,0.0001))


# Create prior
weights=np.genfromtxt('.../omega_GCSL.txt')
spec_z=np.genfromtxt('.../Y_train.txt',delimiter=',')

z_prior=0*zgrid_value

for i in range(len(spec_z)):
    z_prior=z_prior+weights[i]*norm.pdf(zgrid_value, spec_z[i], 0.1)

z_prior /= np.trapz(z_prior, zgrid_value, axis=0)

# Flat prior
flat_prior=0*zgrid_value
flat_prior[:]=1.0/9.0
flat_prior /= np.trapz(flat_prior, zgrid_value, axis=0)

# Combine informative prior and flat prior
z_prior=0.999*z_prior+0.001*flat_prior
z_prior /= np.trapz(z_prior, zgrid_value, axis=0)

# bin_width
bin_width=0*zgrid_value
bin_width[0:100]=0.01# This possibly needs changing
bin_width[100:len(zgrid_value)]=0.04


# Make to array
pzbad_all=np.vstack((z_prior, z_prior))


# Create weights
AGN_chi2=np.loadtxt('.../COS_agn_chi2',delimiter=',')
gal_chi2=np.loadtxt('.../COS_galaxy_chi2',delimiter=',')

# Get template ch2
chi2=np.minimum(gal_chi2,AGN_chi2)
weights_template=np.exp(-0.5*chi2) # Template weight

# Get the different types of uncertainty from GPz
sigma=np.loadtxt('.../sigma.txt')
nu=np.loadtxt('.../nu.txt')
beta_i=np.loadtxt('.../beta_i.txt')
gamma=np.loadtxt('.../gamma.txt')

#
weights_ML=beta_i/(beta_i+nu)
weights_ML=weights_ML[0:n_gal] # ML weight


# Define Hierarchical Bayes function

def HBpz(pzarr, zgrid, pzbad, beta, fbad_min, fbad_max, nbad,specific_ML_weight,specific_template_weight):
            
    fbad_range = np.linspace(fbad_min, fbad_max, nbad) # nbad values to integrate over in Equation 5
    
    pzarr_fbad = np.zeros((n_z, len(fbad_range))) # P(z,f_bad)
    
    for f, fbad in enumerate(fbad_range):
        #print(fbad) # Print what loop on
        
        pzb = (fbad*pzbad) + (pzarr*(1-fbad)) # Equation 3
        #pzarr_fbad[:, f] = np.exp(np.sum(np.log(pzb)/beta, axis=0)) # Equation 4
        pzarr_fbad[:, f] = np.exp(  np.log(pzb[0,:])*specific_ML_weight + np.log(pzb[1,:])*specific_template_weight      ) # Generalised Equation 4 with the beta weights that account for the ML and template reliability
        
    
    # Equation 5
    pzarr_hb = np.trapz(pzarr_fbad, fbad_range, axis=1)
    pzarr_hb /= np.trapz(pzarr_hb, zgrid, axis=0)
    return pzarr_hb



# Set up empty arrays
template_uncertainty=np.zeros(n_gal)
template_best=np.zeros(n_gal)
template_mean=np.zeros(n_gal)
template_median=np.zeros(n_gal)
combination_uncertainty=np.zeros(n_gal)
combination_best=np.zeros(n_gal)
combination_mean=np.zeros(n_gal)
combination_median=np.zeros(n_gal)

template_high_res_best=np.zeros(n_gal)
combination_high_res_best=np.zeros(n_gal)

template_low_res_best=np.zeros(n_gal)
combination_low_res_best=np.zeros(n_gal)



for i in range(n_gal):
    
    print(i)
    
    index=i
    
    specific_ML_pdf=ML_pdfs[index,:] # Find ML pdf
    specific_ML_pdf /= np.trapz(specific_ML_pdf, zgrid_value, axis=0) #Normalise
    specific_template_pdf=template_pdfs[index,:] # Find template pdf
    specific_template_pdf /= np.trapz(specific_template_pdf, zgrid_value, axis=0) # Normalise
    
    # Array of pdfs
    pzarr_specific=np.vstack((specific_ML_pdf, specific_template_pdf))
    
    #Find combination
    specific_combined_pdf=HBpz(pzarr_specific,zgrid_value,pzbad_all,beta_value,fbad_min_value,fbad_max_value,nbad_value,weights_ML[i],weights_template[i])
    
    # Save resulting combination in array
    combined_pdfs[i,:]=specific_combined_pdf
    
    
    
    ## Template point values
    
    # Find mean and mode of array
    index_of_max=np.argmax(specific_template_pdf)
    template_best[i]=zgrid_value[index_of_max]
    template_mean[i]=np.sum(zgrid_value*specific_template_pdf/flat_prior)/np.sum(specific_template_pdf/flat_prior)
    
    # Pick spline best (find mode of spline of pdf which is better)
    template_pdf= spline(zgrid_value,specific_template_pdf,k=3)
    high_res_redshifts=template_pdf(zgrid_high_resolution)
    index_of_max=np.argmax(high_res_redshifts)
    template_high_res_best[i]=zgrid_high_resolution[index_of_max]
    template_low_res_best[i]=template_best[i]
    template_best[i]=template_high_res_best[i]
    
    # Find median
    elements = zgrid_value
    probabilities = specific_template_pdf*bin_width
    probabilities=probabilities/np.sum(probabilities)
    x=np.random.choice(elements, 1000, p=probabilities)
    template_median[i]=np.median(x)    

    
    ## Combination point values
    
    # Find mean and mode of array
    index_of_max=np.argmax(specific_combined_pdf)
    combination_best[i]=zgrid_value[index_of_max]
    combination_mean[i]=np.sum(zgrid_value*specific_combined_pdf/flat_prior)/np.sum(specific_combined_pdf/flat_prior)
    
    # Pick spline best (find mode of spline of pdf which is better)
    template_pdf= spline(zgrid_value,specific_combined_pdf,k=3)
    high_res_redshifts=template_pdf(zgrid_high_resolution)
    index_of_max=np.argmax(high_res_redshifts)
    combination_high_res_best[i]=zgrid_high_resolution[index_of_max]
    combination_low_res_best[i]=combination_best[i]
    combination_best[i]=combination_high_res_best[i]
    
    # Find median
    elements = zgrid_value
    probabilities = specific_combined_pdf*bin_width
    probabilities=probabilities/np.sum(probabilities)
    x=np.random.choice(elements, 1000, p=probabilities)
    combination_median[i]=np.median(x)    
    
    
    
    # Point values of uncertainty for template and combination (standard deviation of pdf, even though both can be multi-modal)
    template_uncertainty[i]=((np.sum((zgrid_value**2)*specific_template_pdf/flat_prior))/np.sum(specific_template_pdf/flat_prior)-(template_mean[i]**2))**0.5
    combination_uncertainty[i]=((np.sum((zgrid_value**2)*specific_combined_pdf/flat_prior))/np.sum(specific_combined_pdf/flat_prior)-(combination_mean[i]**2))**0.5
    
    
    
    
    
# ML point values
ML_photoz=np.genfromtxt('.../GPz_prediction_data_COSMOS.txt',delimiter=',')
ML_uncertainty=ML_photoz[:,5]
ML_best=ML_photoz[:,4]


# Save values
np.savetxt('.../COS_HB_pdfs',combined_pdfs, delimiter=',')

np.savetxt('.../COS_ML_uncertainty',ML_uncertainty, delimiter=',')
np.savetxt('.../COS_ML_best',ML_best, delimiter=',')
np.savetxt('.../COS_template_uncertainty',template_uncertainty, delimiter=',')
np.savetxt('.../COS_template_best',template_best, delimiter=',')
np.savetxt('.../COS_combination_uncertainty',combination_uncertainty, delimiter=',')
np.savetxt('.../COS_combination_best',combination_best, delimiter=',')

np.savetxt('.../COS_template_median',template_median, delimiter=',')
np.savetxt('.../COS_ML_median',ML_median, delimiter=',')
np.savetxt('.../COS_combination_median',combination_median, delimiter=',')

