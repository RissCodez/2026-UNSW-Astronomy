# import libraries
import numpy as np
import pandas as pd
import os
from astropy.io import fits

# first download onedspec and profound sources to Galaxy folder.
# creates masterConfig, ngistInput and ngistOutput folders.
# this function will output the tidied profound sources csv.
arcsec = 0.6
minSNR = 30

# set wd and working cluster
clusterName = 'ABELL0370'
workingDir = r'C:\Users\mares\HDD\Notes\Learning\Tertiary\2026\2026_Summer\Clusters'

def limitClusterSNR(clusterName, workingDir, arcsec, minSNR):
    
    os.chdir(rf'{workingDir}\{clusterName}') # set WD to cluster folder
    
    # setup for later  
    folderNames = ['ngistInput',     # creates dir. for minicubes and masks
                   'ngistOutput',    # creates dir. ngist output files
                   'masterConfig']   # creates dir. to store ngist config files

    for folder in folderNames:
        try: 
            os.mkdir(folder)
        except FileExistsError:
            print(f"{folder} already exists")
        
    # import downloaded profoundsources CSV as df
    clusterData = pd.read_csv(f"{clusterName}_profoundsources.csv")
    
    # extract redshift and redshift confidence values
    redshift = clusterData[['Z']].to_numpy()
    redshiftProb = clusterData[['Z1_PROB']].to_numpy()
    
    # create bool mask for galaxies within MAGPI redshift range (0.2-0.4 inclusive)
    redshiftMask = (redshift >= 0.2) & (redshift <= 0.4)
    redshiftProbMask = redshiftProb >= 0.8
    redshiftBestMask = (redshiftMask & redshiftProbMask)[:,0]
    
    # apply bool mask to df
    redshift0_3 = clusterData[redshiftBestMask]
    
    # find onedspec folder
    clusterNames0_3 = redshift0_3[['MAGPIID']].to_numpy()[:,0]
    
    # init array for loop checking good S/N finishes
    signalNoiseArray = np.zeros(len(clusterNames0_3))
    signalNoiseArrayMax = np.zeros(len(clusterNames0_3))
    signalNoiseArrayMin = np.zeros(len(clusterNames0_3))
    i = 0
    
    # grab only required galaxies from onedspec folder and open FITS file data  
    for galaxy in clusterNames0_3:
        filepath = rf'onedspec\MAGPI{galaxy}_1dspec_{arcsec}arcsec.fits'            
        file = fits.open(filepath)
        DATA = fits.getdata(filepath, 1)
        STAT = fits.getdata(filepath, 2)
        
        # create np array with fits file data as columns
        galaxyData = np.column_stack((DATA, STAT))
        
        # create array of all wavelengths (WL)                                      
        startWL = file[1].header['CRVAL1']                                         
        lengthWL = file[1].header['NAXIS1']
        incrementWL = file[1].header['CDELT1']
        endWL = startWL + lengthWL * incrementWL
        WL = np.arange(startWL, endWL, incrementWL)
        
        # choose 'clean' wavelength range and mask range using bool mask
        maskWL = np.logical_and(WL >= 6000, WL <= 6200)
        galaxyDataClean = galaxyData[maskWL]
        
        # Calculate S/N ratio per Angstrom(0.1nm), DATA = S + N, STAT = N^2
        noise = np.sqrt(galaxyDataClean[:, 1])
        signal = galaxyDataClean[:, 0] - noise
        signalNoiseData = signal/noise
        signalNoiseArray[i] = np.median(signalNoiseData)
        signalNoiseArrayMax[i] = np.max(signalNoiseData)
        signalNoiseArrayMin[i] = np.min(signalNoiseData)
        
        # set new index for loop
        i += 1
    
    # create bool mask of S/N >= minSNR on S/N medians array
    signalNoiseMask = signalNoiseArray >= minSNR
    goodGalaxies = redshift0_3[signalNoiseMask]
    goodGalaxies['SN_MEDIAN'] = signalNoiseArray[signalNoiseMask]
    goodGalaxies['SN_MAX'] = signalNoiseArrayMax[signalNoiseMask]
    goodGalaxies['SN_MIN'] = signalNoiseArrayMin[signalNoiseMask]
    
    # export the filtered data to a new CSV
    goodGalaxies.to_csv(f'{clusterName}_profoundsources_tidy.csv')
    
    originalClusterAmount = len(clusterData)
    goodClusterAmount = len(goodGalaxies)
    return f'{originalClusterAmount} input. {goodClusterAmount} output.'

limitClusterSNR(clusterName, workingDir, arcsec, minSNR)
