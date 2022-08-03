# Create the random points in the field

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from astropy.coordinates import SkyCoord
import random
import numpy as np
import pylab as pb
pb.ion()

axes_numbers_size=22
axes_labels_size=24
marker_size=0.5


num_randoms=1000000

RA_random=np.zeros(num_randoms)
DEC_random=np.zeros(num_randoms)

script_directory='...'
mask_directory='...'



RA_min=33.5
RA_max=37.5

DEC_min=-6
DEC_max=-3.5

## XMM1

mask_name_XMM1='MASK_XMM1_UDEEP_VISTA_v2.fits'
mask_name_XMM2='MASK_XMM2_DEEP_VISTA_v2.fits'
mask_name_XMM3='MASK_XMM3_DEEP_VISTA_v2.fits'


mask_name_SCZ='MASK_XMM1_UDEEP_ch12.fits'


# Plot mask XMM1
hdul1 = fits.open(mask_directory+mask_name_XMM1)
#plt.figure()
#plt.imshow(hdul1[0].data, interpolation='nearest')
fn_1 = get_pkg_data_filename(mask_directory+mask_name_XMM1)
f_1 = fits.open(fn_1)
w_1 = WCS(f_1[0].header)

# Plot mask XMM2
hdul2 = fits.open(mask_directory+mask_name_XMM2)
#plt.figure()
#plt.imshow(hdul2[0].data, interpolation='nearest')
fn_2 = get_pkg_data_filename(mask_directory+mask_name_XMM2)
f_2 = fits.open(fn_2)
w_2 = WCS(f_2[0].header)

# Plot mask XMM1
hdul3 = fits.open(mask_directory+mask_name_XMM3)
#plt.figure()
#plt.imshow(hdul3[0].data, interpolation='nearest')
fn_3 = get_pkg_data_filename(mask_directory+mask_name_XMM3)
f_3 = fits.open(fn_3)
w_3 = WCS(f_3[0].header)


# Taking the HSC masks from /mnt/zfsusers/pwhatfield/vardygroupshare/data/masks/HSC_circle_xmm.reg
ra_circle_1=(2+27.0/60.0+19.86/(60.0*60.0) )*(360.0/24.0)
dec_circle_1=-4-27.0/60.0-47.42/(60.0*60.0)
radius_1=0.736111111111111*u.degree

ra_circle_2=(2+18/60.0+08.827/(60.0*60.0)) *(360.0/24.0)
dec_circle_2=-4-49.0/60.0-09.607/(60.0*60.0)
radius_2=0.81315055*u.degree

ra_circle_3=(2+23.0/60.0+19.029/(60.0*60.0))*(360.0/24.0)
dec_circle_3=-5-19.0/60.0-52.339/(60.0*60.0)
radius_3=0.7344997222*u.degree

ra_circle_4=(2+22.0/60.0+42.795/(60.0*60.0))*(360.0/24.0)
dec_circle_4=-4-3.0/60.0-44.429/(60.0*60.0)
radius_4=0.7777777777777778*u.degree

#circle(02:27:19.86,-04:27:47.42,2650")
#circle(2:18:08.827,-4:49:09.607,2927.342")
#circle(2:23:19.029,-5:19:52.339,2644.199")
#circle(2:22:42.795,-4:03:44.429,2800.000")

coordinates_circle_1=SkyCoord(ra=ra_circle_1*u.degree, dec=dec_circle_1*u.degree, frame='icrs')
coordinates_circle_2=SkyCoord(ra=ra_circle_2*u.degree, dec=dec_circle_2*u.degree, frame='icrs')
coordinates_circle_3=SkyCoord(ra=ra_circle_3*u.degree, dec=dec_circle_3*u.degree, frame='icrs')
coordinates_circle_4=SkyCoord(ra=ra_circle_4*u.degree, dec=dec_circle_4*u.degree, frame='icrs')


print('Starting')

for i in range(num_randoms):
    
    inside_mask=0 # Counter to confirm if outside mask
    
    #print(str(i))
    
    while inside_mask==0:
        
        # Create random coordinates in general patch of sky (accounting for curvature of sky)
        RA_random_trial=RA_min+random.uniform(0, 1)*(RA_max-RA_min)
        DEC_random_trial=(360.0/(2*np.pi))*np.arcsin(np.sin(DEC_max*(2*np.pi/360.0))+random.uniform(0, 1)*(np.sin(DEC_min*(2*np.pi/360))-np.sin(DEC_max*(2*np.pi/360))))
        
        coordinates=SkyCoord(ra=RA_random_trial*u.degree, dec=DEC_random_trial*u.degree, frame='icrs')

        ## Convert Coordinates in Pixels
        check_if_inside_all=0
        check_if_inside_1=0
        check_if_inside_2=0
        check_if_inside_3=0
        
        check_if_inside_circle_1=0
        check_if_inside_circle_2=0
        check_if_inside_circle_3=0
        check_if_inside_circle_4=0
        inside_any_circle=0
        

        #XMM1
        pixels=w_1.world_to_pixel(coordinates)
        x_pixel=int(np.round(w_1.world_to_pixel(coordinates)[1]))
        y_pixel=int(np.round(w_1.world_to_pixel(coordinates)[0]))
        

        if x_pixel>=0 and np.shape(hdul1[0].data)[0]>x_pixel and y_pixel>=0 and np.shape(hdul1[0].data)[1]>y_pixel:
            check_if_inside_1=hdul1[0].data[x_pixel,y_pixel]
        
        #XMM2
        pixels=w_2.world_to_pixel(coordinates)
        x_pixel=int(np.round(w_2.world_to_pixel(coordinates)[1]))
        y_pixel=int(np.round(w_2.world_to_pixel(coordinates)[0]))
        

        if x_pixel>=0 and np.shape(hdul2[0].data)[0]>x_pixel and y_pixel>=0 and np.shape(hdul2[0].data)[1]>y_pixel:
            check_if_inside_2=hdul2[0].data[x_pixel,y_pixel]
        
        #XMM3
        pixels=w_3.world_to_pixel(coordinates)
        x_pixel=int(np.round(w_3.world_to_pixel(coordinates)[1]))
        y_pixel=int(np.round(w_3.world_to_pixel(coordinates)[0]))
        

        if x_pixel>=0 and np.shape(hdul3[0].data)[0]>x_pixel and y_pixel>=0 and np.shape(hdul3[0].data)[1]>y_pixel:
            check_if_inside_3=hdul3[0].data[x_pixel,y_pixel]
            
        #HSC
        if coordinates.separation(coordinates_circle_1)<radius_1:
            check_if_inside_circle_1=1
        
        if coordinates.separation(coordinates_circle_2)<radius_2:
            check_if_inside_circle_2=1
        
        if coordinates.separation(coordinates_circle_3)<radius_3:
            check_if_inside_circle_3=1
        
        if coordinates.separation(coordinates_circle_4)<radius_4:
            check_if_inside_circle_4=1
        
        inside_any_circle=check_if_inside_circle_1+check_if_inside_circle_2+check_if_inside_circle_3+check_if_inside_circle_4

        
        
        #Combine all VISTA
        check_if_inside_all=check_if_inside_1+check_if_inside_2+check_if_inside_3
        
        
        if check_if_inside_all>0 and inside_any_circle>0:
            inside_mask=1
        
    
    RA_random[i]=RA_random_trial
    DEC_random[i]=DEC_random_trial
        




# Plot Random Points
fig=plt.figure(figsize=(12, 6),dpi=96)
#plt.figure()
ax1 = fig.add_subplot(1,1,1)
plt.scatter(RA_random,DEC_random,color=(1,0,0),s=marker_size,label="Galaxies")
plt.gca().invert_xaxis()

# Save data
data=np.zeros([num_randoms,2])
data[:,0]=RA_random
data[:,1]=DEC_random
np.savetxt(script_directory+'output_data/random_coordinates.csv',data, delimiter=',')

