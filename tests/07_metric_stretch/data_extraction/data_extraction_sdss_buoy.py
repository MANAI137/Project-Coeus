"""
PROJECT COEUS: WAYPOINT & BUOY OPERATIONALIZATION
-------------------------------------------------
Author: Miguel Antonio Navarro

DESCRIPTION:
This script performs a 3D Spatial Audit of the SDSS spectroscopic 
catalog to identify 'Isolated Waypoints' (Buoys). Using a cKDTree 
algorithm, it filters for galaxies with zero neighbors within a 
5 Mpc radius. 

These isolated objects serve as pure probes of the manifold's 
metric stretch, as their motion is minimally influenced by local 
gravitational clustering. By flagging these waypoints, we can 
test the 'stiffness' of the metric grain independent of the 
large-scale structure.
"""
import astropy.io.fits as fits
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree

# 1. Configuration
file_path = r'D:\specObj-dr17.fits'
output_csv = 'sdss_waypoint_map.csv'
ISOLATION_RADIUS_MPC = 5.0  # The 5 Mpc "Buoy" Criterion
H0 = 70.0                   # Hubble constant for distance proxy

# Metric Epicenter (CMB Dipole Coordinates)
RA_DIPOLE = 168.0
DEC_DIPOLE = -7.0

def angular_separation(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    a = np.sin((dec2-dec1)/2)**2 + np.cos(dec1)*np.cos(dec2)*np.sin((ra2-ra1)/2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(a)))

def to_cartesian(ra, dec, z):
    """Convert spherical coordinates + redshift to 3D Cartesian Mpc."""
    dist = (z * 3e5) / H0  # Simple linear distance proxy for local isolation
    ra_rad, dec_rad = np.radians(ra), np.radians(dec)
    x = dist * np.cos(dec_rad) * np.cos(ra_rad)
    y = dist * np.cos(dec_rad) * np.sin(ra_rad)
    z_coord = dist * np.sin(dec_rad)
    return np.column_stack((x, y, z_coord))

print(f"Opening {file_path}...")

with fits.open(file_path, memmap=True) as hdul:
    data = hdul[1].data
    
    # Auto-detect columns
    ra_key = next((k for k in data.columns.names if k.upper() in ['RA', 'PLUG_RA']), 'RA')
    dec_key = next((k for k in data.columns.names if k.upper() in ['DEC', 'PLUG_DEC']), 'DEC')
    
    print("Applying initial filters (Galaxy, ZWARNING, Z > 0.01)...")
    mask = (data['CLASS'] == 'GALAXY') & (data['ZWARNING'] == 0) & (data['Z'] > 0.01)
    
    ra = data[ra_key][mask]
    dec = data[dec_key][mask]
    z = data['Z'][mask]

# 2. 3D Spatial Audit (Operationalizing 'Isolation')
print(f"Building 3D Cartesian map for {len(z)} objects...")
coords = to_cartesian(ra, dec, z)
tree = cKDTree(coords)

print(f"Running {ISOLATION_RADIUS_MPC} Mpc neighbor search...")
# count_neighbors finds how many other galaxies are within the radius
# we subtract 1 because the search always finds the galaxy itself
neighbor_count = tree.query_ball_point(coords, r=ISOLATION_RADIUS_MPC, return_sorted=False)
neighbor_count = np.array([len(x) - 1 for x in neighbor_count])

# 3. Flagging the Samples
is_isolated = (neighbor_count == 0)
print(f"Audit Results: Found {np.sum(is_isolated)} Isolated Buoys out of {len(z)} total galaxies.")

# 4. Final Processing
dist_from_epicenter = angular_separation(ra, dec, RA_DIPOLE, DEC_DIPOLE)

df = pd.DataFrame({
    'ra': ra,
    'dec': dec,
    'z': z,
    'angle_to_epicenter': dist_from_epicenter,
    'neighbors_5mpc': neighbor_count,
    'is_buoy': is_isolated.astype(int)  # 1 for isolated, 0 for clustered
})

print(f"Saving {output_csv}...")
df.to_csv(output_csv, index=False)
print("Done! The 'Isolated Waypoint' flag is now operationalized.")