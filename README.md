# PATH_RT_T: CHM-based Radiative Transfer Model for Sloped Terrain

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18937065.svg)](https://doi.org/10.5281/zenodo.18937065)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

PATH_RT_T is a Python implementation of a path-length based radiative transfer model that uses Canopy Height Models (CHM) to simulate vegetation reflectance over sloped terrain. The model explicitly incorporates the 3D structure of forest canopies from CHM data and accounts for terrain effects on radiative transfer processes.

### Key Features

- **3D Canopy Structure**: Uses real or simulated CHM data to represent forest canopy architecture
- **Terrain Effects**: Accounts for slope and aspect effects on incident and reflected radiation
- **Hotspot Modeling**: Implements Kuusk's hotspot factor for accurate near-specular reflectance
- **Multiple Leaf Angle Distributions**: Supports 7 different leaf orientation types (planophile, erectophile, spherical, Campbell's ellipsoidal, etc.)
- **Flexible Input**: Works with custom CHM data or generates homogeneous canopy scenarios
- **Fast Computation**: Optimized ray-tracing algorithm for efficient gap probability calculation

## Model Components

### PATH_RT_T.py
The main radiative transfer model containing:
- `TerrainGeometry`: Class for terrain and canopy geometry management
- `PATH_RT_Terrain`: Main function computing bidirectional reflectance factor (BRF)
- `calc_sun_position`: High-precision solar position calculation
- `Ross_G_function`: Projection functions for various leaf angle distributions
- `area_scatter_phase_function`: Leaf scattering phase function
- `calc_hotspot_factor`: Hotspot effect calculation

### Polar_BRDF_Analysis.py
Script for analyzing BRDF patterns over different slopes and aspects:
- Generates polar plots of BRDF for various terrain configurations
- Visualizes how terrain orientation affects reflectance anisotropy
- Useful for understanding terrain effects on satellite observations

### BRDF_Matrix_Analysis.py
Additional analysis tools for BRDF visualization and sensitivity studies.

## Installation

### Requirements

- Python 3.7+
- NumPy >= 1.19.0
- Matplotlib >= 3.3.0
- SciPy >= 1.5.0
- tifffile >= 2020.9.3

### Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/PATH_RT_T.git
cd PATH_RT_T

# Install dependencies
pip install -r requirements.txt
```

##使用 Usage

### Basic Example

```python
import numpy as np
import PATH_RT_T as MyModel

# Create a homogeneous canopy over sloped terrain
terrain = MyModel.TerrainGeometry(
    mode='homogeneous',
    homogeneous_h=15,      # Canopy height (m)
    homogeneous_base=5,    # Branch-free height (m)
    FAVD=0.4,              # Foliage Area Volume Density
    slope=30,              # Terrain slope (degrees)
    aspect=180,            # Terrain aspect (degrees, 0=North)
    scale=100,             # Scene size (m)
    res=0.5                # Resolution (m)
)

# Define observation geometry
SZA, SAA = 30, 180  # Solar zenith and azimuth angles
VZA, VAA = 20, 0    # Viewing zenith and azimuth angles

# Get geometric components
geo_comps =terrain.get_fast_geometry(
    terrain.sun_mask, SZA, SAA, VZA, VAA
)

# Define optical properties (NIR band example)
tau_l = 0.47    # Leaf transmittance
rho_l = 0.51    # Leaf reflectance
soil_r = 0.1    # Soil reflectance
sky_ratio = 0.09  # Proportion of diffuse radiation

# Calculate BRF
brf = MyModel.PATH_RT_Terrain(
    terrain, tau_l, rho_l, soil_r,
    geo_comps, sky_ratio,
    SZA, SAA, VZA, VAA,
    leaf_class=6,  # Spherical leaf distribution
    Hotspot=0.02   # Hotspot parameter
)

print(f"Calculated BRF: {brf:.4}")
```

### Using Custom CHM Data

```python
# Load your own CHM from GeoTIFF
terrain = MyModel.TerrainGeometry(
    mode='custom',
    tif_path='path/to/your/chm.tif',
    FAVD=0.4,
    slope=25,
    aspect=135
)

# Rest of the code is the same...
```

### Polar BRDF Analysis

```python
from Polar_BRDF_Analysis import plot_polar_brdf_grid

# Analyze BRDF patterns for multiple slopes and aspects
slopes = [0, 15, 30, 45]
aspects = [0, 90, 180, 270]

fig, data = plot_polar_brdf_grid(
    slopes=slopes,
    aspects=aspects,
    sza=30,
    saa=0,
    FAVD=0.25,
    rho_l=0.55,    # NIR reflectance
    tau_l=0.405,   # NIR transmittance
    soil_r=0.211   # NIR soil reflectance
)

fig.savefig('polar_brdf_analysis.png', dpi=300)
```

## Model Theory

The PATH_RT_T model combines:

1. **Geometric-Optical Modeling**: Four-component decomposition (sunlit/shaded canopy and ground)
2. **Radiative Transfer**: Multiple scattering within canopy layers
3. **Terrain Effects**: Correction for sloped surfaces and sky view factor
4. **3D Structure**: Path-length distributions from real CHM data

### Four-Component Model

The model separates reflectance into:
- **Kc**: Sunlit canopy
- **Kt**: Shaded canopy
- **Kg**: Sunlit ground
- **Kz**: Shaded ground

### Gap Probability

Gap probabilities are calculated using:
```
P_gap = exp(-G × FAVD × Path_length)
```
where:
- G: Projection function (depends on leaf angle distribution)
- FAVD: Foliage Area Volume Density (LAI / canopy height)
- Path_length: Ray-tracing distance through canopy

## Citation

If you use this model in your research, please cite:

```bibtex
@software{PATH_RT_T,
  author = {Weihua Li},
  title = {PATH\_RT\_T: CHM-based Radiative Transfer Model for Sloped Terrain},
  year = {2026},
  version = {1.0.1},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.18937065},
  url = {https://github.com/BNUL/PATH_RT_T}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

For questions or issues, please open an issue on GitHub or contact [weihuahk@hku.hk](mailto:weihuahk@hku.hk).

## Acknowledgments

- Based on radiative transfer theory from Ross (1981), Kuusk (1991), and others

## References

- Li, W., Yan, G., Mu, X., Tong, Y., Zhou, K., & Xie, D. (2024). Modeling the hotspot effect for vegetation canopies based on path length distribution. *Remote Sensing of Environment*, 303, 113985. https://doi.org/10.1016/j.rse.2023.113985
- Li, W., Yan, G., Geng, J., Guo, Y., Xie, T., Mu, X., Xie, D., Roujean, J.-L., Zhou, G., & Gastellu-Etchegorry, J.-P. (2025). A model based on spectral invariant theory for correcting topographic effects on vegetation canopy reflectance. *Remote Sensing of Environment*, 322, 114695. https://doi.org/10.1016/j.rse.2025.114695
- Li, W., & Mu, X. (2021). Using fractal dimension to correct clumping effect in leaf area index measurement by digital cover photography. *Agricultural and Forest Meteorology*, 311, 108695. https://doi.org/10.1016/j.agrformet.2021.108695
- Gastellu-Etchegorry, J.-P., Li, W., Yan, G., Cao, B., Kallel, A., Hedman, J., Wang, Y., Zhen, Z., Yin, T., & Lauret, N. (2026). DART 3D radiative transfer modeling applied to RAMI forests. Part 1: Assessing canopy structure effects on directional TIR emissivity. *Journal of Remote Sensing*, 0, 0738. https://doi.org/10.34133/remotesensing.0738
- Li, W., Gastellu-Etchegorry, J.-P., Yin, T., Lauret, N., & Yan, G. (2026). DART 3D radiative transfer modeling applied to RAMI forests. Part 2: LiDAR waveform simulation and canopy structure analysis. *Journal of Remote Sensing*, 0, 0737. https://doi.org/10.34133/remotesensing.0737
- Ross, J. (1981). The radiation regime and architecture of plant stands. The Hague: W. Junk.
- Kuusk, A. (1991). The hot-spot effect in plant canopy reflectance. In R.B. Myneni, & J. Ross (Eds.), *Photon-Vegetation Interactions, Application in Optical Remote Sensing and Plant Ecology* (pp. 139-159). New York: Springer Verlag.
- Li, X.W., Strahler, A.H., & Woodcock, C.E. (1995). A Hybrid Geometric Optical-Radiative Transfer Approach for Modeling Albedo and Directional Reflectance of Discontinuous Canopies. *IEEE Transactions on Geoscience and Remote Sensing*, 33, 466-480. https://doi.org/10.1109/TGRS.1995.8746028
- Chen, J.M., & Leblanc, S.G. (1997). A four-scale bidirectional reflectance model based on canopy architecture. *IEEE Transactions on Geoscience and Remote Sensing*, 35, 1316-1337. https://doi.org/10.1109/36.628798
- Li, X.W., & Strahler, A.H. (1988). Modeling the Gap Probability of a Discontinuous Vegetation Canopy. *IEEE Transactions on Geoscience and Remote Sensing*, 26, 161-170. https://doi.org/10.1109/36.3017

