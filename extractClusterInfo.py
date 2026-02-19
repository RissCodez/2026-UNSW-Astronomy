# -*- coding: utf-8 -*-

# Import Libraries
import os
import numpy as np
import pandas as pd
from astropy.io import fits

# Set Global Values
clusterName = 'ABELL2744'
workingDir = r'C:\Users\mares\HDD\Notes\Learning\Tertiary\2026\2026_Summer\Clusters'

def extractClusterInfo(clusterName, workingDir):
    
    # set working directory to cluster folder
    os.chdir(rf'{workingDir}\{clusterName}')
    
    inputData = [f for f in os.listdir('ngistInput') 
                 if f.endswith('_minicube.fits')]
    
    coordArray = np.zeros((len(inputData), 2))
    
    for i, minicube in enumerate(inputData):
        file = fits.open(f'ngistInput/{minicube}')
        
        XCEN = file[1].header['NAXIS1'] / 2 + 0.5
        YCEN = file[1].header['NAXIS2'] / 2 + 0.5
        
        coordArray[i, 0] = XCEN  # Column 0 = X center
        coordArray[i, 1] = YCEN  # Column 1 = Y center
        
        file.close()
        print(f"{minicube}: ({XCEN:.1f}, {YCEN:.1f})")
    
    print(f"\ncoordArray shape: {coordArray.shape}")
    print("First 5 rows:")
    print(coordArray[:5])
    
    clusterDataTidy = pd.read_csv(f"{clusterName}_profoundsources_tidy.csv")
    
    galaxyInfo = clusterDataTidy[['MAGPIID']]
    galaxyInfo['Z'] = clusterDataTidy[['Z']]
    galaxyInfo['XCEN'] = coordArray[:, 0]
    galaxyInfo['YCEN'] = coordArray[:, 1]
    
    galaxyInfo.to_csv(f'{clusterName}_masterConfigValues.csv')
    
extractClusterInfo(clusterName, workingDir)
