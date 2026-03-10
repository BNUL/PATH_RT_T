"""
Basic Example: Simulating BRF for a Homogeneous Canopy on Sloped Terrain

This example demonstrates how to use the PATH_RT_T model to simulate
bidirectional reflectance factor (BRF) for a simple homogeneous forest
canopy on sloped terrain.
"""

import numpy as np
import matplotlib.pyplot as plt
import PATH_RT_T as MyModel

# ============================================================================
# Configuration
# ============================================================================

# Terrain parameters
SLOPE = 30      # Terrain slope (degrees)
ASPECT = 180    # Terrain aspect (degrees, 0=North, 90=East, 180=South, 270=West)

# Canopy structure parameters
CANOPY_HEIGHT = 15  # Tree height (m)
BRANCH_HEIGHT = 5   # Branch-free height (m)
FAVD = 0.4         # Foliage Area Volume Density (LAI / canopy height)

# Solar geometry
SZA = 30        # Solar zenith angle (degrees)
SAA = 180       # Solar azimuth angle (degrees, 0=North)

# Optical properties (Near-Infrared band example)
RHO_L = 0.51    # Leaf reflectance
TAU_L = 0.47    # Leaf transmittance
SOIL_R = 0.1    # Soil reflectance
SKY_RATIO = 0.09  # Proportion of diffuse radiation

# Leaf angle distribution
LEAF_CLASS = 6  # 6 = Spherical distribution
HOTSPOT = 0.02  # Hotspot parameter (related to leaf size)

# ============================================================================
# Create Terrain
# ============================================================================

print("Creating terrain...")
terrain = MyModel.TerrainGeometry(
    mode='homogeneous',
    homogeneous_h=CANOPY_HEIGHT,
    homogeneous_base=BRANCH_HEIGHT,
    FAVD=FAVD,
    slope=SLOPE,
    aspect=ASPECT,
    scale=100,      # Scene size (m)
    res=0.5         # Spatial resolution (m)
)

print(f"Terrain created:")
print(f"  LAI: {terrain.LAI:.2f}")
print(f"  Mean canopy height: {terrain.mean_h:.2f} m")
print(f"  Sky View Factor (SVF): {terrain.svf:.3f}")

# ============================================================================
# Calculate BRF for a Single View Angle
# ============================================================================

VZA = 20    # Viewing zenith angle (degrees)
VAA = 0     # Viewing azimuth angle (degrees)

print(f"\nCalculating BRF...")
print(f"  Solar geometry: SZA={SZA}°, SAA={SAA}°")
print(f"  View geometry: VZA={VZA}°, VAA={VAA}°")

# Get geometric components
geo_comps = terrain.get_fast_geometry(
    terrain.sun_mask, SZA, SAA, VZA, VAA
)

# Calculate BRF
brf = MyModel.PATH_RT_Terrain(
    terrain, TAU_L, RHO_L, SOIL_R,
    geo_comps, SKY_RATIO,
    SZA, SAA, VZA, VAA,
    leaf_class=LEAF_CLASS,
    Hotspot=HOTSPOT
)

print(f"\nResults:")
print(f"  BRF: {brf:.4f}")
print(f"  Geometric components (Kc, Kt, Kg, Kz, Kgz_s): {geo_comps}")

# ============================================================================
# Calculate BRF for Multiple View Angles (Principal Plane)
# ============================================================================

print(f"\nCalculating BRF in principal plane...")

# View zenith angles from -75° to +75° (negative = backward direction)
vza_range = np.arange(-75, 76, 5)
brf_values = []

for vza in vza_range:
    if vza < 0:
        # Backward scattering direction
        vza_actual = abs(vza)
        vaa_actual = SAA  # Same side as sun
    else:
        # Forward scattering direction
        vza_actual = vza
        vaa_actual = SAA + 180  # Opposite side from sun

    try:
        geo_comps = terrain.get_fast_geometry(
            terrain.sun_mask, SZA, SAA, vza_actual, vaa_actual
        )
        
        brf = MyModel.PATH_RT_Terrain(
            terrain, TAU_L, RHO_L, SOIL_R,
            geo_comps, SKY_RATIO,
            SZA, SAA, vza_actual, vaa_actual,
            leaf_class=LEAF_CLASS,
            Hotspot=HOTSPOT
        )
        brf_values.append(brf)
    except:
        brf_values.append(np.nan)

brf_values = np.array(brf_values)

# ============================================================================
# Plot Results
# ============================================================================

plt.figure(figsize=(10, 6))

plt.plot(vza_range, brf_values, 'o-', linewidth=2, markersize=6, 
         color='darkgreen', label='NIR Band')

# Mark hotspot
hotspot_idx = np.argmin(np.abs(vza_range - (-SZA)))
plt.plot(vza_range[hotspot_idx], brf_values[hotspot_idx], '*', 
         markersize=15, color='red', label='Hotspot')

plt.axvline(x=-SZA, color='orange', linestyle='--', alpha=0.5, 
            label=f'Solar Position (SZA={SZA}°)')
plt.axvline(x=0, color='gray', linestyle=':', alpha=0.5)

plt.xlabel('View Zenith Angle (degrees)\n[Negative=Backward, Positive=Forward]', 
           fontsize=12)
plt.ylabel('BRF', fontsize=12)
plt.title(f'Bidirectional Reflectance Factor - Principal Plane\n' +
          f'Slope={SLOPE}°, Aspect={ASPECT}°, SZA={SZA}°', 
          fontsize=13, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.legend(fontsize=11)
plt.tight_layout()

output_file = 'example_brf_principal_plane.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nPlot saved: {output_file}")

plt.show()

print("\nExample completed successfully!")
