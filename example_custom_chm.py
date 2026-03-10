"""
Advanced Example: Using Custom CHM Data

This example demonstrates how to use the PATH_RT_T model with your own
Canopy Height Model (CHM) data from a GeoTIFF file.
"""

import numpy as np
import matplotlib.pyplot as plt
import PATH_RT_T as MyModel
import tifffile

# ============================================================================
# Configuration
# ============================================================================

# Path to your CHM file (GeoTIFF format)
CHM_FILE = 'path/to/your/chm.tif'  # Update this path!

# Terrain parameters
SLOPE = 25      # Manually specified slope (degrees)
ASPECT = 135    # Manually specified aspect (degrees)
FAVD = 0.4      # Foliage Area Volume Density

# Solar geometry
SZA = 35        # Solar zenith angle (degrees)
SAA = 180       # Solar azimuth angle (degrees)

# Optical properties (Red band example)
RHO_L = 0.10    # Leaf reflectance (Red band)
TAU_L = 0.075   # Leaf transmittance (Red band)
SOIL_R = 0.065  # Soil reflectance (Red band)
SKY_RATIO = 0.1 # Proportion of diffuse radiation

# Leaf angle distribution
LEAF_CLASS = 6  # Spherical distribution
HOTSPOT = 0.02

# ============================================================================
# Load and Inspect CHM Data
# ============================================================================

print("Loading CHM data...")
try:
    with tifffile.TiffFile(CHM_FILE) as tif:
        chm = tif.asarray()
        
        # Get metadata
        tags = tif.pages[0].tags
        if 'ModelPixelScaleTag' in tags:
            pixel_scale = tags['ModelPixelScaleTag'].value
            resolution = pixel_scale[0]
            print(f"  Resolution: {resolution} m")
        else:
            resolution = None
            print("  Resolution: Unknown (not in GeoTIFF tags)")
        
    print(f"  Shape: {chm.shape}")
    print(f"  Height range: {np.nanmin(chm):.2f} - {np.nanmax(chm):.2f} m")
    print(f"  Mean height: {np.nanmean(chm[chm > 0]):.2f} m")
    
except FileNotFoundError:
    print(f"ERROR: CHM file not found: {CHM_FILE}")
    print("\nCreating a simulated CHM for demonstration...")
    
    # Create a simulated CHM
    chm = np.zeros((200, 200), dtype=np.float32)
    np.random.seed(42)
    
    # Add some random trees
    for i in range(50):
        r, c = np.random.randint(0, 200, 2)
        h = np.random.uniform(10, 20)
        radius = np.random.uniform(3, 6)
        
        y, x = np.ogrid[:200, :200]
        mask = ((y - r)**2 + (x - c)**2) <= (radius / 0.5)**2
        chm[mask] = np.maximum(chm[mask], h)
    
    CHM_FILE = 'simulated_chm.tif'
    tifffile.imwrite(CHM_FILE, chm)
    resolution = 0.5
    print(f"  Simulated CHM created and saved to: {CHM_FILE}")

# ============================================================================
# Create Terrain with Custom CHM
# ============================================================================

print("\nCreating terrain from CHM...")
terrain = MyModel.TerrainGeometry(
    mode='custom',
    tif_path=CHM_FILE,
    FAVD=FAVD,
    slope=SLOPE,
    aspect=ASPECT,
    margin=0  # No margin pixels to exclude
)

print(f"Terrain created:")
print(f"  Size: {terrain.rows} x {terrain.cols} pixels")
print(f"  Resolution: {terrain.res:.2f} m")
print(f"  LAI: {terrain.LAI:.2f}")
print(f"  Mean canopy height: {terrain.mean_h:.2f} m")
print(f"  Max canopy height: {terrain.max_h:.2f} m")
print(f"  Vegetation cover: {terrain.veg_cover:.1%}")
print(f"  Canopy continuity index: {terrain.c_index:.3f}")

# ============================================================================
# Calculate BRF for Different View Angles
# ============================================================================

print(f"\nCalculating BRF for multiple view angles...")

# Define view angle grid
vza_array = np.arange(0, 61, 10)  # 0, 10, 20, ..., 60 degrees
vaa_array = [0, 90, 180, 270]      # North, East, South, West

results = {}

for vaa in vaa_array:
    brf_list = []
    for vza in vza_array:
        try:
            geo_comps = terrain.get_fast_geometry(
                terrain.sun_mask, SZA, SAA, vza, vaa
            )
            
            brf = MyModel.PATH_RT_Terrain(
                terrain, TAU_L, RHO_L, SOIL_R,
                geo_comps, SKY_RATIO,
                SZA, SAA, vza, vaa,
                leaf_class=LEAF_CLASS,
                Hotspot=HOTSPOT
            )
            brf_list.append(brf)
        except:
            brf_list.append(np.nan)
    
    results[vaa] = np.array(brf_list)
    print(f"  VAA={vaa:3d}°: BRF range = {np.nanmin(brf_list):.4f} - {np.nanmax(brf_list):.4f}")

# ============================================================================
# Visualize Results
# ============================================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: CHM visualization
ax = axes[0]
im = ax.imshow(terrain.chm_top, cmap='YlGn', interpolation='nearest')
ax.set_title('Canopy Height Model (CHM)', fontsize=12, fontweight='bold')
ax.set_xlabel('X (pixels)')
ax.set_ylabel('Y (pixels)')
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Height (m)', fontsize=11)

# Plot 2: BRF vs View Zenith Angle
ax = axes[1]
colors = ['blue', 'green', 'red', 'orange']
labels = ['North (0°)', 'East (90°)', 'South (180°)', 'West (270°)']

for i, (vaa, color, label) in enumerate(zip(vaa_array, colors, labels)):
    ax.plot(vza_array, results[vaa], 'o-', linewidth=2, markersize=7,
            color=color, label=label)

ax.set_xlabel('View Zenith Angle (degrees)', fontsize=11)
ax.set_ylabel('BRF', fontsize=11)
ax.set_title(f'BRF vs View Angle\n(SZA={SZA}°, SAA={SAA}°, Red Band)', 
             fontsize=12, fontweight='bold')
ax.legend(fontsize=10, loc='best')
ax.grid(True, alpha=0.3)

plt.tight_layout()
output_file = 'example_custom_chm.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nPlot saved: {output_file}")

plt.show()

print("\nAdvanced example completed successfully!")
