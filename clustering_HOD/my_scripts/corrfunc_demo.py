# A demo to check Corrfunc working, see https://github.com/manodeep/Corrfunc

from __future__ import print_function
import os.path as path
import numpy as np
import Corrfunc
from Corrfunc.theory import wp

# Setup the problem for wp
boxsize = 500.0
pimax = 40.0
nthreads = 4

# Create a fake data-set.
Npts = 100000
x = np.float32(np.random.random(Npts))
y = np.float32(np.random.random(Npts))
z = np.float32(np.random.random(Npts))
x *= boxsize
y *= boxsize
z *= boxsize

# Setup the bins
rmin = 0.1
rmax = 20.0
nbins = 20

# Create the bins
rbins = np.logspace(np.log10(0.1), np.log10(rmax), nbins + 1)

# Call wp
wp_results = wp(boxsize, pimax, nthreads, rbins, x, y, z, verbose=True, output_rpavg=True)

# Print the results
print("#############################################################################")
print("##       rmin           rmax            rpavg             wp            npairs")
print("#############################################################################")
print(wp_results)