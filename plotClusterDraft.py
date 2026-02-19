# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 11:04:18 2026

@author: mares
"""

import numpy as np
import pandas as pd
import os
from astropy.io import fits
import matplotlib.pyplot as plt

g1 = 'ABELL0370'
g2 = 'ABELL2744'
g3 = 'MACS0257'
g4 = 'MACS0416NE'
g5 = 'MACS0416S'
g6 = 'MACS0940'

#LAMBDA

# set global values
workingDir = r'C:\Users\mares\HDD\Notes\Learning\Tertiary\2026\2026_Summer\Clusters'
clusterName = g2

# set wd
os.chdir(rf'{workingDir}\{clusterName}')

# grabs the names of all galaxies in output folder
outputList = [galaxy for galaxy in os.listdir('ngistOutput') if os.path.isdir(os.path.join('ngistOutput', galaxy))]

# read profoundsources, extract ellipticit
clusterData = pd.read_csv(f'{clusterName}_profoundsources_tidy.csv')
axrat_rt = np.array(clusterData.axrat_rt)
R50_rt = np.array(clusterData.R50_rt)

rotationList = np.zeros_like(range(len(clusterData)))
thresholdList = []
lambdaList = []
ellipticity = (1 - axrat_rt)

# filter out those not spatially resolved - i didnt apply it cuz theyre all good...
spatiallyResolved = np.array(clusterData.semimin_rt > 0.65)

for i in range(len(outputList)):
    # reads output of spatial binning ngist module                                  Columns: ID Spaxel ID | BIN_ID Bin ID | X Y Pixel coordinates of spaxels | XBIN YBIN Pixel coordinates of bins | FLUX Flux in the spaxel | SNR Signal-to-noise ratio in the spaxel | SNRBIN Signal-to-noise ratio in the bin | NSPAX Number of spaxels in the bin
    spatialFile = fits.open(f'ngistOutput/{outputList[i]}/{outputList[i]}_table.fits')
    spatialData = spatialFile[1].data
    
    # reads output of stellar kinematics ngist module                               Columns: V SIGMA H3 H4 etc Stellar kinematics. The number of columns here will correspond to value of the MOM keyword. For example, MOM: 4 will return V SIGMA H3 H4, while MOM:6 will return V SIGMA H3 H4 H5 H6 | ERR_* Errors on the stellar kinematics from MC-simulations | FORM_ERR_* Formal errors on the stellar kinematics | REDDENING Extinction calculated by PPXF
    kinFile = fits.open(f'ngistOutput/{outputList[i]}/{outputList[i]}_kin.fits')
    kinData = kinFile[1].data
    
    # find indexes for all unmasked bins (pos) from binning table
    goodBinMask = np.where(spatialData.BIN_ID >= 0)[0]
    # abs neg (masked) and pos bins to match masked and unmasked spaxels
    allSpaxelBins = np.abs(spatialData.BIN_ID)
    # grab index of unmasked spaxels
    goodBinIndex = allSpaxelBins[goodBinMask]
    
    V = kinData['V'][goodBinIndex]
    
    v_min = np.nanmin(V)
    v_max = np.nanmax(V)
    v_shift = (v_min + v_max) / 2
    V = V - v_shift
    
    SIGMA = kinData['SIGMA'][goodBinIndex]
    FLUX = spatialData['FLUX'][goodBinMask][goodBinIndex]
    
    # mask relevant X,Y coords
    X = spatialData['X'][goodBinMask]
    Y = spatialData['Y'][goodBinMask]
    pixelSize = spatialFile[0].header['PIXSIZE']
    # calculate direct distance from centre pixel and convert to dist in arcsec 
    R = ( np.sqrt(X**2 + Y**2) ) * pixelSize # MUSE 0.1765................arcsec
    
    fluxMask = R < R50_rt[i]
    
    FLUX = FLUX[fluxMask]
    R = R[fluxMask]
    SIGMA = SIGMA[fluxMask]
    V = V[fluxMask]
    
    # calculate lambda
    lambdaR = (np.sum(FLUX*R*np.abs(V)) / np.sum(FLUX*R*np.sqrt(np.square(V) + np.square(SIGMA))) )
    lambdaList.append(lambdaR)
                          #still add effective half radius bit 
                          
    thresholdList.append((0.31) * np.sqrt(ellipticity[i]))
    rotationList[i] = lambdaR > (0.31) * np.sqrt(ellipticity[i])
    
    print(f'for {i}, pixels masked = {np.sum(fluxMask) - len(fluxMask)}')
    
rotationList = rotationList.astype(bool)
# FR = TRUE
# SR = FALSE

workingDir = r'C:\Users\mares\HDD\Notes\Learning\Tertiary\2026\2026_Summer\Clusters'
plotCluster = True
evenAxes = False
legendLocation = 'lower right'
legendSize = 15
plotSize = 8,15

# set working directory to cluster folder
os.chdir(rf'{workingDir}\{clusterName}')

if evenAxes == True:
    axes = "equal"
if evenAxes == False: 
    axes = "auto"

# read in data
clusterRot = pd.read_csv(f"{clusterName}_kinematicStates.csv") # rotation state
clusterTidy = pd.read_csv(f"{clusterName}_profoundsources_tidy.csv") # tidy data
clusterAll = pd.read_csv(f"{clusterName}_profoundsources.csv") # all Z = 0.3
clusterAll = clusterAll[(clusterAll['Z'] > 0.2) & (clusterAll['Z'] < 0.4)]

# get coords for analysed galaxies and others in cluster
tidyRA = np.array(clusterTidy.RAcen_rt) / 206264.806247
tidyDEC = np.array(clusterTidy.Deccen_rt) / 206264.806247
allRA = np.array(clusterAll.RAcen_rt) / 206264.806247
allDEC = np.array(clusterAll.Deccen_rt) / 206264.806247

# Find the index of max flux for centre of coord system
maxFluxI = np.argmax(clusterAll.flux_rt)
originRA = allRA[maxFluxI]
originDEC = allDEC[maxFluxI]

# init array for loop to find x,y coords for plotting, for analysed and not
arrayTidyRA = np.zeros_like(tidyRA)
arrayTidyDEC = np.zeros_like(tidyDEC)
arrayAllRA = np.zeros_like(allRA)
arrayAllDEC = np.zeros_like(allDEC)

for i in range(len(tidyRA)):
    X = (tidyRA[i] - originRA) * np.cos(originDEC)   # RA
    Y = tidyDEC[i] - originDEC                       # DEC
    arrayTidyRA[i] = -X                              # mirrored on the x-axis
    arrayTidyDEC[i] = Y

for i in range(len(allRA)):
    X = (allRA[i] - originRA) * np.cos(originDEC)   # RA
    Y = allDEC[i] - originDEC                       # DEC
    arrayAllRA[i] = -X                              # mirrored on the x-axis
    arrayAllDEC[i] = Y

tidyPlot = np.column_stack((arrayTidyRA, arrayTidyDEC, clusterRot[['ROTATION']]))
rot = tidyPlot[clusterRot.ROTATION == 'rotating']
nonrot = tidyPlot[clusterRot.ROTATION == 'nonrotating']
SR = tidyPlot[rotationList == False]
FR = tidyPlot[rotationList == True]

# set working directory to cluster folder
os.chdir(f'{workingDir}')

#plt.style.use('dark_background')
#plt.style.use('default') 

# plt.figure(figsize = (plotSize))
# plt.scatter(arrayAllRA, arrayAllDEC, c='gray', label='Other Galaxies') 
# plt.scatter(rot[:,0], rot[:,1], c='red', label='Rotating') 
# plt.scatter(nonrot[:,0], nonrot[:,1], label='Non-Rotating') 
# plt.xlabel('X')
# plt.ylabel('Y')
# #plt.grid(True, alpha=0.2)
# plt.title(f'{clusterName} Cluster Kinematics', fontsize = 20)
# plt.axis(axes)
# plt.legend(loc=legendLocation, fontsize = legendSize)
# plt.savefig(f'{clusterName}_kinematicStates_visual_N.jpg', dpi=300, bbox_inches='tight')
# plt.show()

F_SR = (len(SR)/(len(SR)+len(FR)))*100

plt.style.use('dark_background')
plt.figure()
plt.scatter(arrayAllRA, arrayAllDEC, c='#575757', label='NA') 
plt.scatter(FR[:,0], FR[:,1], c='#f74f4f', label='FR') 
plt.scatter(SR[:,0], SR[:,1], c='#4fcaf7', label='SR') 
plt.xlabel('X')
plt.ylabel('Y')
#plt.grid(True, alpha=0.2)
plt.title(f'{clusterName} Cluster Kinematics', fontsize = 20)
plt.axis(axes)
plt.legend(loc=legendLocation, fontsize = legendSize)
plt.savefig(f'{clusterName}_kinematicStates_lambda_N.jpg', dpi=300, bbox_inches='tight')
plt.show()





# plot with density as grid behind

plt.style.use('dark_background')
plt.figure()

# create histogram of density in grid of binsize. calculates the edges of
# the plot based on the plotsize&bin amount (12 = 12+1)
densityData, xEnd, yEnd = np.histogram2d(arrayAllRA, arrayAllDEC, bins=20)

heatmap = plt.pcolormesh(xEnd, yEnd, densityData.T, cmap = 'gray', alpha = 0.7)
plt.colorbar(heatmap, label = 'Cluster Density (count of galaxies in 20x20 grid)')

# creates plot similar to b4
plt.scatter(arrayAllRA, arrayAllDEC, c='#575757', label='NA') 
plt.scatter(FR[:,0], FR[:,1], c = '#f74f4f', label = 'FR') 
plt.scatter(SR[:,0], SR[:,1], c = '#4fcaf7', label = 'SR') 
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'{clusterName} Cluster Kinematics', fontsize=20)
plt.axis(axes)
plt.legend(loc=legendLocation, fontsize=legendSize)
#plt.savefig(f'{clusterName}_density_grid.jpg', dpi=300, bbox_inches='tight')
#plt.show()


# same thing but with the blended histogram grid and no unanalysed galaxies

plt.figure()
im = plt.imshow(densityData.T,
                cmap='gray',
                extent = [xEnd[0], xEnd[-1], yEnd[0], yEnd[-1]],
                origin ='lower',
                interpolation = 'bilinear',
                alpha = 0.7)
plt.colorbar(im, label='Cluster Density (count of galaxies in 20x20 grid)')

# plt.scatter(arrayAllRA, arrayAllDEC, c='#575757', label='NA') 
plt.scatter(FR[:,0], FR[:,1], c='#f74f4f', label='FR') 
plt.scatter(SR[:,0], SR[:,1], c='#4fcaf7', label='SR') 
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'{clusterName} Cluster Kinematics', fontsize=20)
plt.axis(axes)
plt.legend(loc=legendLocation, fontsize=legendSize)
#plt.savefig(f'{clusterName}_density_blurred.jpg', dpi=300, bbox_inches='tight')
#plt.show()

# density of dots estimates kernel
# probability density function


import seaborn as sns
import matplotlib.pyplot as plt

plt.style.use('default')
plt.figure(figsize = (8,6))

kdeObject = sns.kdeplot(x=arrayAllRA, y=arrayAllDEC, fill=True, bw_adjust=0.8, 
                       cmap='gray_r', alpha=0.6)
plt.scatter(arrayAllRA, arrayAllDEC, c='black', label='NA', s=6) 
plt.scatter(FR[:,0], FR[:,1], c='#f74f4f', label='FR', s=6)
plt.scatter(SR[:,0], SR[:,1], c='#4fcaf7', label='SR', s=6)

cbar = plt.colorbar(kdeObject.collections[0], ax=plt.gca())
cbar.set_label('Density Histogram Scale')

plt.xlabel('X'); plt.ylabel('Y')
plt.title(f'{clusterName}')
plt.axis(axes)
plt.legend(loc=legendLocation, fontsize=8)
plt.savefig(f'{clusterName}_kinematics_density.jpg', dpi=300, bbox_inches='tight')
plt.show()

lambdaList = np.array(lambdaList)
thresholdList = np.array(thresholdList)

thresholdList = np.sort(thresholdList)
ellipticitySorted = np.sort(ellipticity)

plt.figure(figsize = (10,10))
plt.plot(ellipticitySorted, thresholdList, zorder=1)
plt.scatter(x=ellipticity[rotationList == True], y=lambdaList[rotationList == True],zorder=2)
plt.scatter(x=ellipticity[rotationList == False], y=lambdaList[rotationList == False],zorder=2)
plt.ylabel('lambdaR')
plt.xlabel('ellipticity')
plt.savefig(f'{clusterName}_lambda_ellip.jpg', dpi=300, bbox_inches='tight')
plt.show()



