# This is an old code to use MCMC to fit a HOD model to ACF data


# Import everything needed
import halomod as hm
import numpy as np
import pylab as pb
pb.ion()
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import cosmolopy.distance as cd
import cosmolopy.parameters as cp
from scipy.integrate import simps
import triangle
import random
import time
import emcee
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='sans-serif')

##################################################################


# Read in data [columns = theta, w, dw, ....] 
massbin=6 # Mass data to fit
zbin=4 # Redshift data to fit

little_h=0.6711 #0.71

mass_median=[1,2,3,4,5,6,7]

z_median=[0.62,0.88,1.12,1.48,1.96] # Medians of the redshifts
volumes=[1211788,1821023,2317136,5023644,7561368] # Calculated with volume_calc.py
x_low=[1808,2400,3200,3800,4700] # Ranges to do deprojection on
x_high=[2750,3400,4030,4850,5700] # Ranges to do deprojection on

MCMC_num_steps=1000 #1000 # Number of MCMC steps

alpha_max=2.5 # Max alpha value for prior
alpha_min=0.5  #0.2 # Min alpha value for prior
M_min_min=10 # Min M_min value for prior
M_min_max=15 # Max M_min value for prior
M_1_max=17 # Max M_1 value for prior
M_0_min=8 # Min M_0 value for prior
sigma_min=0.00001 # Min sigma value for prior
sigma_max=0.6 # Max sigma value for prior

mask_percent=0.97 # Fraction of field not covered by mask

version=1 # MCMC version
version_acf=1 # ACF versions

print 'version '+str(version)

# Number of galaxies before being modified for mask
num_gal_sample=np.zeros(shape=(5,7)) # Redshift first index, mass second
num_gal_sample[0,0]=6535
num_gal_sample[0,1]=5061
num_gal_sample[0,2]=3877
num_gal_sample[0,3]=2847
num_gal_sample[0,4]=1884
num_gal_sample[0,5]=1022
num_gal_sample[0,6]=321

num_gal_sample[1,1]=9791
num_gal_sample[1,2]=7365
num_gal_sample[1,3]=5453
num_gal_sample[1,4]=3824
num_gal_sample[1,5]=2330
num_gal_sample[1,6]=1023

num_gal_sample[2,2]=7512
num_gal_sample[2,3]=5529
num_gal_sample[2,4]=3892
num_gal_sample[2,5]=2412
num_gal_sample[2,6]=1064

num_gal_sample[3,3]=10800
num_gal_sample[3,4]=5875
num_gal_sample[3,5]=2542
num_gal_sample[3,6]=752

num_gal_sample[4,4]=5165
num_gal_sample[4,5]=1331

# Cosmic variance estimates from Trenti
cv_errors=np.zeros(shape=(5,7)) # Redshift first index, mass second
cv_errors[0,0]=582
cv_errors[0,1]=465
cv_errors[0,2]=367
cv_errors[0,3]=282
cv_errors[0,4]=198
cv_errors[0,5]=119
cv_errors[0,6]=48

cv_errors[1,1]=774
cv_errors[1,2]=605
cv_errors[1,3]=467
cv_errors[1,4]=345
cv_errors[1,5]=235
cv_errors[1,6]=120

cv_errors[2,2]=598
cv_errors[2,3]=462
cv_errors[2,4]=344
cv_errors[2,5]=231
cv_errors[2,6]=120

cv_errors[3,3]=662
cv_errors[3,4]=404
cv_errors[3,5]=205
cv_errors[3,6]=79

cv_errors[4,4]=343
cv_errors[4,5]=117

# Calculate density
n_data_gal_den=num_gal_sample[zbin-1,massbin-1]/(volumes[zbin-1]*mask_percent)
cv_n_error=cv_errors[zbin-1,massbin-1]/volumes[zbin-1]


# Print steps
print 'z='+str(zbin)+'   m='+str(massbin) + '   Number of steps =' + str(MCMC_num_steps)

# Load RR values
RR= np.genfromtxt('.../TPCF3_python_RR.dat')
theta_RR= np.genfromtxt('.../TPCF3_python_theta.dat')

# Load acf
d_long = np.genfromtxt(".../TPCF3_python_data_function_zbin"+str(zbin)+"_massbin"+str(massbin)+"_v"+str(version_acf)+".dat")
d=d_long[0:300:10] # Pick every tenth value
wh = (d[:,0] < -1) & (d[:,0] > -3) # where theta<0.6, for fitting

# Do everything on restricted set
x = d[wh,0]    # theta
y = d[wh,1]    # w(theta)
yerr_log_std=d[wh,2] # error on w(theta)
yerr_std=(10**(y+yerr_log_std))-(10**(y)) # Move out of log space
yerr=(10**(y+yerr_log))-(10**(y)) # Move out of log space
x=10**x # Move out of log space
y=10**y # Move out of log space

# Do everything on whole set for plotting purposes
x_long = d_long[:,0]    # theta
y_long = d_long[:,1]    # w(theta)
yerr_log_long_std=d_long[:,2] # error on w(theta)
yerr_long_std=(10**(y_long+yerr_log_long_std))-(10**(y_long)) # Move out of log space
yerr_long=(10**(y_long+yerr_log_long))-(10**(y_long)) # Move out of log space
x_long=10**x_long # Move out of log space
y_long=10**y_long # Move out of log space


# Define own acf model directly

# NB h.r comes out in Mpc/h , so must be corrected. cd.comoving_distance comes out in just Mpc, as does angular correction etc

def integral(x,u,theta):
    xi = spline((1.0/little_h)*h.r,h.corr_gal, k=3) # Spatial correlation function to be integrated
    integral=2*f(x)*f(x)*xi(np.sqrt((u**2)+((x*theta)**2))) # Spatial correlation function to be integrated
    return integral

def angular_corr_gal_PH(h,f, theta_min, theta_max, theta_num,x_min,x_max):
    # Confirm angular range to calculate makes sense
    if theta_min <= 0 or theta_min >= theta_max:
        raise ValueError("theta_min must be > 0 and < theta_max")
    if theta_max >= 180.0:
        raise ValueError("theta_max must be < pi")
    if x_min <= 0 or x_min >= x_max:
        raise ValueError("x_min must be >0 and < x_max")
    
    theta = 10 ** np.linspace(np.log10(theta_min), np.log10(theta_max), theta_num) # Set up theta values
    w = np.zeros_like(theta) # Set up empty matrix for acf
    
    # Range for dummy variable integral
    umin = -4
    umax = np.log10(120)
    u = np.logspace(umin, umax, 100)
    
    # Range for comoving distance variable integral
    x = np.logspace(np.log10(x_min), np.log10(x_max), 10)
    
    integral_along_x=0*u
    
    # Do the integral
    for i, th in enumerate(theta):
            
            for j, u_val in enumerate(u):
                integral_along_x[j]=simps(integral(x,u_val,th),x)
            
            w[i] = simps(integral_along_x,u)
    
    return w, theta


# Read in N(z) and create function
nz = np.loadtxt(".../TPCF3_python_data_redshift_zbin"+str(zbin)+"_massbin"+str(massbin)+"_v"+str(version_acf)+".dat")
if nz.ndim == 2: # Normal for z values and number values
    s = spline(nz[:,0],nz[:,1],k=3) # Create spline
elif nz.ndim ==1:
    print 'Potential error in redshift input'
    hist,edges = np.histogram(nz,bins=100,density=True) #Normalized histogram with 100 bins
    # Create a function s(z)
    centres = (edges[1:] + edges[:-1])/2 # Get mid-points of bins
    s = spline(centres,hist,k=3)


 # Set up z/ inversion
cosmo = cp.WMAP7_ML(flat=True)  # Set up the WMAP7 cosmology
z = np.linspace(0,6,100)
r_of_z = cd.comoving_distance(z,**cosmo) # get the comoving distances
z_of_r = spline(r_of_z,z,k=3) # call spline with inverted axes

r_of_z_spline = spline(z,r_of_z,k=3)
dr_dz_spline=r_of_z_spline.derivative()
# now can define final f(r) 
def f(r):
    z = z_of_r(r)  # This is the redshift corresponding to the radius
     # This changes the variables (including the derrivative), and also 
    return np.maximum(s(z)/dr_dz_spline(z),0*s(z))/simps(np.maximum(s(z_of_r(np.linspace(0,8000,800)))/dr_dz_spline(z_of_r(np.linspace(0,8000,800))),0*s(z_of_r(np.linspace(0,8000,800)))),np.linspace(0,8000,800))

###############################################################


# Define a halo model
h = hm.HaloModel(rmin=1e-3,rmax=320.0,rnum=50,z=z_median[zbin-1],hod_model="Zheng05")
#h.update(omegav=0.73)
#h.update(omegab=0.045)
#h.update(omegam=0.27)
#h.update(n=0.96)
#h.update(sigma_8=0.8)
#h.update(h=0.71)
#h.update(cm_relation='bullock_rescaled') # Set concentration relation
#h.update(mf_fit='Jenkins') # Set hmf

# Define log likelihood, prior and probabilities for model fitting
def lnlike(theta, x, y, yerr,n_data_gal_den,cv_n_error):
    M_min,M_1,alpha,sig_logm,M_0  = theta # HOD parameters
    h.update(hod_params={"M_min":M_min, "M_1":M_1, "alpha":alpha, "sig_logm":sig_logm, "M_0":M_0}) # Update HOD model with new parameters
    w,xm = angular_corr_gal_PH(h,f,theta_min=0.001*np.pi/180,theta_max=1.0*np.pi/180.0,theta_num=30,x_min=x_low[zbin-1],x_max=x_high[zbin-1]) # Calculate acf model at range of values
    s = spline(xm,w,k=3) # Crate spline of the acf
    model = s(x*np.pi/180.0) # Define model w(theta) over same theta values as data
    C=np.sum(RR*s(theta_RR*np.pi/180.0))/np.sum(RR) # Calculate an estimate of the integral constraint
    # Calculate the log-liklihood
    inv_sigma2 = 1.0/(yerr**2)
    prov_output=-0.5*(np.sum(((y-(model-C))**2)*inv_sigma2 - np.log(inv_sigma2)))
    prov_output=prov_output+(-0.5*((((h.mean_gal_den*(little_h**3))-n_data_gal_den)/cv_n_error)**2))-np.log(cv_n_error) # Add the number counts agreement
    return prov_output # Log likelihood

def lnprior(theta):
    M_min,M_1,alpha,sig_logm,M_0 = theta # HOD parameters
    # Flat priors over some range...
    if M_min_min < M_min < M_min_max and M_min < M_1 < M_1_max and alpha_min < alpha < alpha_max and sigma_min < sig_logm < sigma_max  and M_0_min < M_0 < M_1: 
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr,n_data_gal_den,cv_n_error):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr,n_data_gal_den,cv_n_error)





ndim, nwalkers = 5, 48 # No. of parameters to be fitted for, and no. of MCMC chains (must be >2*ndim, I think...)

# Set up initial positions of walkers randomly in the allowed parameter space
pos=[] # Empty list

for i in range(nwalkers):
    rand1=random.uniform(M_min_min, M_min_max)
    rand2=random.uniform(rand1, M_1_max)
    rand3=random.uniform(alpha_min, alpha_max)
    rand4=random.uniform(sigma_min, sigma_max)
    rand5=random.uniform(M_min_min, rand2)
    pos.append([rand1,rand2,rand3,rand4,rand5])
    
print 'done'

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y,yerr,n_data_gal_den,cv_n_error)) # Define MCMC sampler

t1 = time.time() # Time the fitting process
sampler.run_mcmc(pos,MCMC_num_steps) # Run MCMC fit
t2 = time.time()
print str(np.floor((t2-t1)/60))+' minutes' # Prints how long the MCMC fitting took

samples = sampler.chain.reshape((-1, ndim)) # Combine chains into one array (consider burn-in)
lnprob = sampler.lnprobability.reshape(-1) # Probability values for each point in chain

############ Now plot stuff

## Triangle plot

plt.figure()
triangle.corner(samples, labels=[r'$M_\mathrm{min}$', r'$M_1$', r'$\alpha$',r'$\sigma$',r'$M_0$'],quantiles=[0.16, 0.5, 0.84]) # Triangle plot of parameters
plt.title(r'$z='+str(z_median[zbin-1])+'$ , $m='+str(mass_median[massbin-1])+'$', fontsize=15)

#Save triangle plot
pb.savefig(".../hod_triangle_z_" + str(zbin) + "_m_" + str(massbin) + "_numstep_" + str(MCMC_num_steps) +'_v'+str(version)+ ".png")
print 'Masses at this stage still in M/h'

## ACF plot and saving derrived parameters

sample_rate=10 # How much of the full sample to use in these calculations

model_omega_distrib=np.zeros(shape=(30,nwalkers*MCMC_num_steps/sample_rate)) # Calculate acf for each parameter sampling
bias_effective_distrib=np.zeros(shape=(1,nwalkers*MCMC_num_steps/sample_rate)) # Bias for distribution
mass_effective_distrib=np.zeros(shape=(1,nwalkers*MCMC_num_steps/sample_rate)) # Effective halo mass for distribution
satellite_fraction_distrib=np.zeros(shape=(1,nwalkers*MCMC_num_steps/sample_rate)) # Satellite fraction for distribution
r0_distrib=np.zeros(shape=(1,nwalkers*MCMC_num_steps/sample_rate)) # r0 for distribution
number_density_distrib=np.zeros(shape=(1,nwalkers*MCMC_num_steps/sample_rate)) # Galaxy comoving density for distribution


for i in range(nwalkers*MCMC_num_steps/sample_rate): # Go through the 
    h.update(hod_params={"M_min":samples[i*sample_rate,0],"M_1":samples[i*sample_rate,1],"alpha":samples[i*sample_rate,2],"sig_logm":samples[i*sample_rate,3],"M_0":samples[i*sample_rate,4]}) # Update from the samples
    # theta is in radians
    w,theta = angular_corr_gal_PH(h,f,theta_min=0.001*np.pi/180,theta_max=1.*np.pi/180.0,theta_num=30,x_min=x_low[zbin-1],x_max=x_high[zbin-1]) # Calculate acf
    # Save all the derrived parameters
    model_omega_distrib[:,i]=w
    bias_effective_distrib[0,i]=h.bias_effective
    satellite_fraction_distrib[0,i]=h.satellite_fraction
    mass_effective_distrib[0,i]=h.mass_effective-np.log10(little_h) #This turns it from M/h to just M
    number_density_distrib[0,i]=h.mean_gal_den*(little_h**3) #This turns it from h^3/Mpc^3 to just Mpc^-3
    
    r = h.r*(1.0/little_h) #This turns it from Mpc/h to just Mpc
    xi = h.corr_gal
    r=r[0:39]
    xi=xi[0:39]
    r=r[::-1]
    xi=xi[::-1]
    r_of_xi = spline(np.log10(xi),np.log10(r),k=3)
    r0_distrib[0,i]=10**r_of_xi(0)

# Print the derrived parameters
print 'Satellite Fraction'
print np.percentile(satellite_fraction_distrib[0,:],16), np.percentile(satellite_fraction_distrib[0,:],50), np.percentile(satellite_fraction_distrib[0,:],84)
print 'Bias'
print np.percentile(bias_effective_distrib[0,:],16), np.percentile(bias_effective_distrib[0,:],50), np.percentile(bias_effective_distrib[0,:],84)
print 'Effective Halo Mass'
print np.percentile(mass_effective_distrib[0,:],16), np.percentile(mass_effective_distrib[0,:],50), np.percentile(mass_effective_distrib[0,:],84)
print 'r0'
print np.percentile(r0_distrib[0,:],16), np.percentile(r0_distrib[0,:],50), np.percentile(r0_distrib[0,:],84)
print 'Derrived Number Density'
print np.percentile(number_density_distrib[0,:],16), np.percentile(number_density_distrib[0,:],50), np.percentile(number_density_distrib[0,:],84)
print 'Observed Number Density'
print n_data_gal_den


# Create objects for holding the models
model_best=np.zeros(shape=(30,1))
model_lower=np.zeros(shape=(30,1))
model_upper=np.zeros(shape=(30,1))

# Find percentiles of acf
for i in range(30):
    model_best[i,0]=np.percentile(model_omega_distrib[i,:],50)
    model_lower[i,0]=np.percentile(model_omega_distrib[i,:],16)
    model_upper[i,0]=np.percentile(model_omega_distrib[i,:],84)

# Plot data
plt.figure()
pb.plot(x_long,y_long,color=(1.0-((7.0-massbin)/6.0),0,((7.0-massbin)/6.0)),linewidth=2.0) # Plot data
pb.plot(x_long,y_long+yerr_long,color=(1.0-((7.0-massbin)/6.0),0,((7.0-massbin)/6.0)),linewidth=1.0) # Plot data
pb.plot(x_long,y_long-yerr_long,color=(1.0-((7.0-massbin)/6.0),0,((7.0-massbin)/6.0)),linewidth=1.0) # Plot data
#plt.errorbar(x,y,yerr=yerr,linestyle="None",color='k')


# Plot model fits
xm=theta # theta and xm in radians

# Calculate integral constraint for models

g = spline(xm,model_best,k=3)
C=np.sum(RR*g(theta_RR*np.pi/180.0))/np.sum(RR)
pb.plot(xm*180/np.pi,model_best-C,'k',linewidth=2.0) # Plot HOD model c.f.

g = spline(xm,model_lower,k=3)
C=np.sum(RR*g(theta_RR*np.pi/180.0))/np.sum(RR)
pb.plot(xm*180/np.pi,model_lower-C,'k',linewidth=1.0) # Plot HOD model c.f.

g = spline(xm,model_upper,k=3)
C=np.sum(RR*g(theta_RR*np.pi/180.0))/np.sum(RR)
pb.plot(xm*180/np.pi,model_upper-C,'k',linewidth=1.0) # Plot HOD model c.f.

# Set up plot etc
pb.xscale("log")
pb.yscale("log")
pb.xlim(0.001,0.1)
pb.ylim(0.01,10**(0.2))
pb.legend(loc='lower left',frameon=False, fontsize=11)
pb.xlabel(r'$\theta \mathrm{ (degrees)}$')
pb.ylabel(r'$w(\theta)')
plt.text(10**(-2), 10**(-0.4), r'$z='+str(z_median[zbin-1])+'$ , $m='+str(mass_median[massbin-1])+'$', fontsize=15)

# Save plot
pb.savefig(".../hod_plot_z_" + str(zbin) + "_m_" + str(massbin) + "_numstep_" + str(MCMC_num_steps) +'_v'+str(version)+ ".png")

# Save data
with open(".../hod_samples_z_" + str(zbin) + "_m_" + str(massbin) + "_numstep_" + str(MCMC_num_steps) +'_v'+str(version)+ ".dat", "w") as out_file:
    out_string = "zbin, massbin, MCMC_num_steps, alpha_max, alpha_min, M_min_min, M_min_max, M_1_max, M_0_min, sigma_min,sigma_max, nwalkers"
    out_string += "\n"
    out_string += str(zbin) + ',' + str(massbin) + ',' + str(MCMC_num_steps) + ','+ str(alpha_max) + ','+str(alpha_min) + ','+ str(M_min_min) + ','+','+ str(M_min_max) +','+ str(M_1_max)+','+str(M_0_min)+','+str(sigma_min)+','+str(sigma_max)+','+str(nwalkers)
    out_string += "\n"
    out_file.write(out_string)
    for i in range(nwalkers*MCMC_num_steps):
        out_string = ""
        out_string += str(samples[i,0]-np.log10(little_h))
        out_string += ','
        out_string += str(samples[i,1]-np.log10(little_h))
        out_string += ','
        out_string += str(samples[i,2])
        out_string += ','
        out_string += str(samples[i,3])
        out_string += ','
        out_string += str(samples[i,4]-np.log10(little_h))
        out_string += "\n"
        out_file.write(out_string)
        

with open(".../hod_samples_deduced_params_z_" + str(zbin) + "_m_" + str(massbin) + "_numstep_" + str(MCMC_num_steps) +'_v'+str(version)+ ".dat", "w") as out_file:
    out_string = "zbin, massbin, MCMC_num_steps, alpha_max, alpha_min, M_min_min, M_min_max, M_1_max, M_0_min, sigma_min,sigma_max, nwalkers"
    out_string += "\n"
    out_string += "bias, satellite fraction, effective halo mass, galaxy comoving density, r0"
    out_string += "\n"
    out_file.write(out_string)
    for i in range(nwalkers*MCMC_num_steps/sample_rate):
        out_string = ""
        out_string += str(bias_effective_distrib[0,i])
        out_string += ','
        out_string += str(satellite_fraction_distrib[0,i])
        out_string += ','
        out_string += str(mass_effective_distrib[0,i])
        out_string += ','
        out_string += str(number_density_distrib[0,i])
        out_string += ','
        out_string += str(r0_distrib[0,i])
        out_string += "\n"
        out_file.write(out_string)

# Save data
with open(".../acf_quantiles_samples_z_" + str(zbin) + "_m_" + str(massbin) + "_numstep_" + str(MCMC_num_steps) +'_v'+str(version)+ ".dat", "w") as out_file:
    out_string = "model best"
    out_string += "\n"
    out_string += str(model_best)
    out_string += "\n"
    out_string += str(model_upper)
    out_string += "\n"
    out_string += str(model_lower)
    out_string += "\n"
    out_file.write(out_string)