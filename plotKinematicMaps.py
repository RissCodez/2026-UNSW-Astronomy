# Import Libraries
import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt

# Set Global Values
clusterName = 'ABELL2744'
workingDir = r'C:\Users\mares\HDD\Notes\Learning\Tertiary\2026\2026_Summer\Clusters'

def plotKinematicMaps(clusterName, workingDir):
    
    # set working directory to cluster folder
    os.chdir(rf'{workingDir}\{clusterName}')
    
    # grabs the names of all galaxies
    outputList = [d for d in os.listdir('ngistOutput') if os.path.isdir(os.path.join('ngistOutput', d))]
    
    failedFiles = []
    
    for galaxy in outputList:
        
        if os.path.exists(f"ngistOutput/{galaxy}/{galaxy}_table.fits") == False or os.path.exists(f"ngistOutput/{galaxy}/{galaxy}_kin.fits") == False:
            failedFiles.append(galaxy)
            continue
    
        # read fits files with V
        tableFile = f"ngistOutput/{galaxy}/{galaxy}_table.fits"
        kinFile = f"ngistOutput/{galaxy}/{galaxy}_kin.fits"
        
        # Read table and get pixel size
        tableHDU = fits.open(tableFile)
        table = tableHDU[1].data
        pixelSize = tableHDU[0].header["PIXSIZE"]
        tableHDU.close()
        
        # Read kin stuff
        kinHDU = fits.open(kinFile)
        kinResults = kinHDU[1].data
        kinHDU.close()
        
        # Get bin ids + convert
        _, idxConvertShortToLong = np.unique(np.abs(table.BIN_ID), return_inverse=True)
        
        # Get velocity values for Voronoi bins only (BIN_ID >= 0)
        idx_voronoi = np.where(table.BIN_ID >= 0)[0]
        velocities = kinResults['V'][idxConvertShortToLong[idx_voronoi]]
        
        v_min = np.nanmin(velocities)
        v_max = np.nanmax(velocities)
        
        # v_range_sides = (abs(v_min) + abs(v_max)) / 2
        # v_adjustment = abs(v_min) - v_range_sides
        v_shift = (v_min + v_max) / 2
        velocities_shifted = velocities - v_shift
        
        # Create pixel image
        xmin = np.nanmin(table.X[idx_voronoi]) - 1
        xmax = np.nanmax(table.X[idx_voronoi]) + 1
        ymin = np.nanmin(table.Y[idx_voronoi]) - 1
        ymax = np.nanmax(table.Y[idx_voronoi]) + 1
        
        npixels_x = int(np.round((xmax - xmin) / pixelSize) + 1)
        npixels_y = int(np.round((ymax - ymin) / pixelSize) + 1)
        
        i = np.array(np.round((table.X[idx_voronoi] - xmin) / pixelSize), dtype=int)
        j = np.array(np.round((table.Y[idx_voronoi] - ymin) / pixelSize), dtype=int)
        
        image = np.full((npixels_x, npixels_y), np.nan)
        image[i, j] = velocities_shifted
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 10), dpi=150)
    
        im = ax.imshow(
            np.rot90(image),
            cmap='RdBu_r',  # Red-Blue diverging colormap good for velocities
            vmin=np.nanpercentile(velocities_shifted, 5),
            vmax=np.nanpercentile(velocities_shifted, 95),
            interpolation="none",
            extent=[
                xmin - pixelSize / 2,
                xmax + pixelSize / 2,
                ymin - pixelSize / 2,
                ymax + pixelSize / 2,
            ]
        )
        
        # Colourbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("L.O.S. Stellar Velocity V [km/s]", fontsize=12)
        cbar.ax.tick_params(labelsize=10)
        
        # Labels and title
        ax.set_title("V (Voronoi Binned)", fontsize=14, pad=20)
        ax.set_xlabel("X [arcsec]", fontsize=12)
        ax.set_ylabel("Y [arcsec]", fontsize=12)
        
        plt.tight_layout()
        
        # Save plot
        output_file = rf"ngistOutput\{galaxy}_V.png"
        if os.path.exists(output_file):
            print(f"Skipping {output_file} (already exists)")
            plt.close()
            continue
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"Velocity map saved as: {output_file}")
    
    print(failedFiles)
    return failedFiles

plotKinematicMaps(clusterName, workingDir)
