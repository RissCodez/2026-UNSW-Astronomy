# Import Libraries
import os
from astropy.io import fits

# Set Global Values
clusterName = 'ABELL2744'
workingDir = r'C:\Users\mares\HDD\Notes\Learning\Tertiary\2026\2026_Summer\Clusters'

def extractGalaxyMasks(clusterName, workingDir):
    
    os.chdir(rf'{workingDir}\{clusterName}')
    
    inputData = [f for f in os.listdir('ngistInput')  # created in first script
                 if f.endswith('_minicube.fits')]
    
    for minicube in inputData:
        galaxy = minicube.replace('_minicube.fits', '')
        mask = os.path.join('ngistInput', f"{galaxy}_mask.fits")
        
        # Open FITS file
        cube_path = os.path.join('ngistInput', minicube)
        with fits.open(cube_path) as extension:
            maskData = extension[7].data
        
        # Save mask
        maskFile = fits.PrimaryHDU(maskData)
        maskFile.header['PIXSIZE'] = 0.2
        maskFile.writeto(mask, overwrite=True)
        
        print(f"{mask}")
        
extractGalaxyMasks(clusterName, workingDir)
