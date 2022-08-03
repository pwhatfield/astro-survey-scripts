# Calculate comoving volume corresponding to a survey and redshift range

import numpy as np
import pylab as pb
pb.ion()

import cosmolopy.distance as cd
import cosmolopy.parameters as cp


survey_area_deg2=12
z_upper=0.6
z_lower=0.3
number_gal=3118

cosmo = cp.WMAP7_ML(flat=True)

volume=(cd.comoving_volume(z_upper, **cosmo)-cd.comoving_volume(z_lower, **cosmo))*(survey_area_deg2/41253.0)
ng=number_gal/volume
