import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import shift, binary_fill_holes
import warnings
import tifffile
from datetime import datetime, timedelta
warnings.filterwarnings("ignore")

# ==============================================================================
# Part 1: functions
# ==============================================================================

def cosd(theta): return np.cos(np.deg2rad(theta))
def sind(theta): return np.sin(np.deg2rad(theta))
def tand(theta): return np.tan(np.deg2rad(theta))
def acosd(val): return np.degrees(np.arccos(np.clip(val, -1, 1)))

def Campbell_g(theta, ALA):
    """Campbell (1990)"""
    theta_rad = np.deg2rad(theta)
    if ALA <= 0:
        return 0.5
    alpha = np.deg2rad(ALA)
    chi = np.power((alpha/9.65), -(1/1.65)) - 3
    tmp = np.sqrt(chi**2 * np.sin(theta_rad)**2 + np.cos(theta_rad)**2)
    tmp = np.maximum(np.abs(tmp), 1e-6)
    Lambda = chi + 1.744 * np.power((chi + 1.182), -0.733)
    g_val = (2 * (chi**3) * np.sin(theta_rad)) / (Lambda *
                                                  (np.cos(theta_rad)**2 + (chi**2) * np.sin(theta_rad)**2)**2)
    if np.ndim(g_val) == 0:
        if np.isnan(g_val):
            g_val = 0.5
    else:
        g_val[np.isnan(g_val)] = 0.5
    return g_val
def get_gFun(iorien, theta_L):
    """
    Calculate leaf angle distribution function g_L(theta_L)
    
    Parameters:
        iorien: Leaf normal distribution type (1~6)
            1 - planophile     (a=1,  b=2)
            2 - erectophile    (a=-1, b=-2)
            3 - plagiophile    (a=-1, b=4)
            4 - extremophile   (a=1,  b=4)
            5 - uniform        (a=0,  b=arbitrary)
            6 - spherical      g(θ) = sin(θ)
        theta_L: Leaf inclination angle (radians), can be scalar or numpy array
    Returns:
        g_L: Corresponding leaf angle distribution function value
    """
    
    if iorien == 6:  # spherical
        g_L = np.sin(theta_L)
        
    else:
        # 根据类型设定 a, b 参数
        if iorien == 1:    # planophile
            a, b = 1, 2
        elif iorien == 2:  # erectophile
            a, b = -1, -2   # 注意：很多文献这里用 b=2，此处保留原MATLAB逻辑
        elif iorien == 3:  # plagiophile
            a, b = -1, 4
        elif iorien == 4:  # extremophile
            a, b = 1, 4
        elif iorien == 5:  # uniform
            a, b = 0, 0     # b实际上不起作用
        else:
            raise ValueError("iorien must be an integer between 1 and 6")
            
        g_L = (2 / np.pi) * (1 + a * np.cos(b * theta_L))
    
    return g_L

def get_APF(a, b, rol, taul, cteta, phi, cteta1, phi1):
    """
    Calculate Area Projection Factor (APF) for bidirectional reflection/transmission
    Numerical integration method consistent with original MATLAB code
    
    Parameters:
        a, b         : Leaf angle distribution parameters (for 1 + a*cos(b*theta_L))
        rol          : Reflection coefficient
        taul         : Transmission coefficient
        cteta, phi   : Cosine of zenith angle & azimuth angle (degrees) for incident direction
        cteta1, phi1 : Cosine of zenith angle & azimuth angle (degrees) for viewing/scattering direction
    
    Returns:
        get_APF      : Calculated area projection factor value
    """
    # Preprocess angles
    cteta1 = cteta1 + 180          # Important: 180-degree reversal
    
    ccteta  = np.cos(np.deg2rad(cteta))     # cos(θ)
    ccteta1 = np.cos(np.deg2rad(cteta1))    # cos(θ₁)
    
    steta   = np.sin(np.deg2rad(cteta))     # sin(θ)
    steta1  = np.sin(np.deg2rad(cteta1))    # sin(θ₁)
    
    phi_rad  = np.deg2rad(phi)
    phi1_rad = np.deg2rad(phi1)
    
    kpi3 = 1.0 / (np.pi ** 3)
    
    # Integration grid (consistent with original code)
    n = 30
    m = 4 * n
    h_theta = 0.5 * np.pi / n
    h_fi    = 2.0 * np.pi / m
    
    theta_i = 0.5 * h_theta      # Start from midpoint
    fi_1    = 0.5 * h_fi
    
    integral = 0.0
    
    for i in range(n):
        fi_j = fi_1
        xx = 0.0
        c_i = np.cos(theta_i)
        s_i = np.sin(theta_i)
        
        for j in range(m):
            # Cosine of angles between two directions and leaf normal
            yy = ccteta  * c_i + steta  * s_i * np.cos(phi_rad  - fi_j)
            zz = ccteta1 * c_i + steta1 * s_i * np.cos(phi1_rad - fi_j)
            
            # Product
            zz *= yy
            
            # Choose reflection or transmission based on sign
            if zz <= 0:
                xx += rol * abs(zz)
            else:
                xx += taul * zz
            
            fi_j += h_fi
        
        xx *= h_fi
        yy = 1.0 + a * np.cos(b * theta_i)     # Note: This is a simplified form
        integral += yy * xx
        
        theta_i += h_theta
    
    integral *= h_theta
    get_APF = kpi3 * integral
    
    return get_APF
def Ross_G_function(theta, phi, leaf_orientation_class, ALA):
    if leaf_orientation_class == 7:
        return Campbell_g(theta, ALA)
    elif leaf_orientation_class == 6:
        return 0.5
    else:
        return get_G(leaf_orientation_class, theta, phi)
def get_G(iorien, theta, fi):
    """
    Calculate integral: G(Omega) = (1/(2*pi)) ∫∫ g_L(θ_L) |cos(Ω, Ω_L)| dθ_L dφ_L
    Parameters:
        iorien: Leaf normal distribution type (1~6), same as get_gFun
        theta: Solar/viewing zenith angle (degrees)
        fi: Solar/viewing azimuth angle (degrees)
    Returns:
        G_Fun: Projection function G(Omega)
    Note: Uses simple numerical integration with Simpson's rule or trapezoidal integration
          (original code uses rectangular integration)
    """
    # Convert to radians
    theta = np.deg2rad(theta)
    fi = np.deg2rad(fi)
    
    # Integration grid setup (consistent with original MATLAB)
    n = 30
    m = 4 * n
    h_theta = 0.5 * np.pi / n
    h_fi = 2 * np.pi / m
    
    theta_i = 0.5 * h_theta
    fi_1 = 0.5 * h_fi
    
    G_Fun = 0.0
    
    for i in range(n):
        fi_j = fi_1
        c_i = np.cos(theta_i)
        s_i = np.sin(theta_i)
        xx = 0.0
        
        for j in range(m):
            # |cos(Ω · Ω_L)| = |cosθ cosθ_L + sinθ sinθ_L cos(φ - φ_L)|
            yy = np.cos(theta) * c_i + np.sin(theta) * s_i * np.cos(fi - fi_j)
            xx += np.abs(yy)
            fi_j += h_fi
        
        xx *= h_fi
        g_val = get_gFun(iorien, theta_i)  # Call previously defined get_gFun
        G_Fun += g_val * xx
        theta_i += h_theta
    
    G_Fun *= h_theta
    G_Fun /= (2 * np.pi)
    
    return G_Fun
def area_scatter_phase_function(rho_l, tau_l, theta_0, phi_0, theta, phi, leaf_class, ALA=57.5):
    a = 0
    b = 0
    if leaf_class == 6:
        sza_r, vza_r = np.deg2rad(theta_0+180), np.deg2rad(theta)
        diff_phi = np.deg2rad(phi - phi_0)
        cos_beta = np.cos(sza_r)*np.cos(vza_r) + np.sin(sza_r) * \
            np.sin(vza_r)*np.cos(diff_phi)
        # cos_beta = np.clip(cos_beta, -1, 1)
        beta = np.arccos(cos_beta)
        omega = rho_l + tau_l
        func = (np.sin(beta) - beta * cos_beta) / np.pi
        Gamma_val = omega * func / 3.0 + (tau_l * cos_beta) / 3.0
        return np.abs(Gamma_val)
    elif leaf_class == 5: # uniform
        a = 0   
    elif leaf_class == 4: # extremophile
        a = 1
        b = 4
    elif leaf_class == 3: # plagiophile
        a = -1
        b = 4
    elif leaf_class == 2: # erectophile
        a = -1
        b = 2   
    elif leaf_class == 1: # planophile
        a = 1
        b = 2
    Gamma_val = get_APF(a,b,rho_l,tau_l,theta_0,phi_0,theta,phi) * np.pi
    return Gamma_val
def calc_hotspot_factor(par, SZA, SAA, VZA, VAA, Gs, Gv, uL, z, ls, lv):
    """
    Calculate Hot Spot Factor (Chs) in Kuusk model
    
    Parameters:
        par   : Average leaf size (bl, leaf size parameter)
        SZA   : Solar zenith angle (degrees)
        SAA   : Solar azimuth angle (degrees)
        VZA   : Viewing zenith angle (degrees)
        VAA   : Viewing azimuth angle (degrees)
        Gs    : G function value for solar direction (projection function)
        Gv    : G function value for viewing direction
        uL    : Average leaf area density (LAI / canopy height) or FAVD
        z     : Relative height (from canopy top downward, same unit as total canopy height)
        ls    : Relative optical path per unit canopy depth in solar direction (mu0)
        lv    : Relative optical path per unit canopy depth in viewing direction (muv)
    
    Returns:
        Chs   : Hot spot factor value at height z
    """
    # Convert angles to radians
    SZA = np.deg2rad(SZA)
    SAA = np.deg2rad(SAA)
    VZA = np.deg2rad(VZA)
    VAA = np.deg2rad(VAA)
    
    mu0 = 1/ls
    muv = 1/lv

    # f1 component
    f1 = np.sqrt(Gs * Gv / (mu0 * muv))
    
    # Calculate relative azimuth angle phi and normalize to [0, 360)
    phi = VAA - SAA
    if phi < 0:
        phi += 2 * np.pi
    if phi > 2 * np.pi:
        phi -= 2 * np.pi
    
    # Calculate cos(gamma) - cosine of angle between two directions
    cosgamma = (np.cos(SZA) * np.cos(VZA) +
                np.sin(SZA) * np.sin(VZA) * np.cos(VAA - SAA))
    
    # Calculate delta (key parameter in Kuusk hotspot model)
    delta = np.sqrt(1/np.cos(SZA)**2 + 1/np.cos(VZA)**2 -
                    2 * cosgamma / (np.cos(SZA) * np.cos(VZA)))
    
    # Prevent division by zero or very small values
    delta = max(delta, 1e-5)   # Original code uses 0.00001
    
    # Calculate f2
    f2 = (uL * par / delta) * (1 - np.exp(-z * delta / par))
    
    # Final hotspot factor
    Chs = np.exp(f1 * f2)
    
    return Chs


def calc_gap_prob(path_array, favd, G):
    if path_array.size == 0:
        return 1.0
    return np.mean(np.exp(-G * favd * path_array))


def calc_sun_position(lat, lon, dt, timezone=8, pressure=1010, temp=10):
    
    # --- 1. Time processing (convert to UTC) ---
    if dt.tzinfo is None:
        # If Naive Time, subtract timezone to get UTC
        dt_utc = dt - timedelta(hours=timezone)
    else:
        # If Aware Time, convert directly to UTC
        dt_utc = dt.astimezone(timezone.utc)

    # Calculate Julian Day (Julian Day)
    def get_julian_day(d):
        Y, M, D = d.year, d.month, d.day
        h, m, s = d.hour, d.minute, d.second + d.microsecond / 1e6
        if M <= 2:
            Y -= 1
            M += 12
        A = int(Y / 100)
        B = 2 - A + int(A / 4)
        JD = int(365.25 * (Y + 4716)) + int(30.6001 * (M + 1)) + D + B - 1524.5
        JD += (h + m / 60.0 + s / 3600.0) / 24.0
        return JD

    jd = get_julian_day(dt_utc)
    t = (jd - 2451545.0) / 36525.0  # Julian century from J2000

    # Solar mean longitude (L0)
    L0 = 280.46646 + 36000.76983 * t + 0.0003032 * t**2
    L0 = L0 % 360
    
    # Solar mean anomaly (M)
    M = 357.52911 + 35999.05029 * t - 0.0001537 * t**2
    M_rad = np.radians(M)
    
    # Orbital eccentricity (e)
    e = 0.016708634 - 0.000042037 * t - 0.0000001267 * t**2
    
    # Solar equation of center (C)
    C = (1.914602 - 0.004817 * t - 0.000014 * t**2) * np.sin(M_rad) + \
        (0.019993 - 0.000101 * t) * np.sin(2 * M_rad) + \
        0.000289 * np.sin(3 * M_rad)
    
    # Solar true longitude (True Longitude, theta)
    theta = L0 + C
    
    omega = 125.04 - 1934.136 * t
    lambda_sun = theta - 0.00569 - 0.00478 * np.sin(np.radians(omega))
    lambda_rad = np.radians(lambda_sun)
    
    # Mean obliquity of ecliptic (Mean Obliquity, epsilon0)
    epsilon0 = 23 + 26/60.0 + 21.448/3600.0 - 46.8150/3600.0 * t - \
               0.00059/3600.0 * t**2 + 0.001813/3600.0 * t**3
    
    epsilon = epsilon0 + 0.00256 * np.cos(np.radians(omega))
    eps_rad = np.radians(epsilon)
    

    alpha_rad = np.arctan2(np.cos(eps_rad) * np.sin(lambda_rad), np.cos(lambda_rad))
    # Declination (delta)
    delta_rad = np.arcsin(np.sin(eps_rad) * np.sin(lambda_rad))
    

    GMST0 = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + \
            0.000387933 * t**2 - (t**3) / 38710000
    GMST0 = GMST0 % 360
    # Apparent sidereal time needs nutation correction (omitted here for small errors <0.001 degrees)
    GAST = GMST0 
    
    # --- 5. Convert to horizontal coordinate system (Az, El) ---
    # Hour angle (H)
    H = (GAST + lon) - np.degrees(alpha_rad)
    H = H % 360
    if H > 180: H -= 360
    H_rad = np.radians(H)
    
    lat_rad = np.radians(lat)
    
    # Elevation angle (Elevation)
    sin_el = np.sin(lat_rad) * np.sin(delta_rad) + np.cos(lat_rad) * np.cos(delta_rad) * np.cos(H_rad)
    el_rad = np.arcsin(np.clip(sin_el, -1, 1))
    el_deg = np.degrees(el_rad)
    
    # Geometric zenith angle (Geometric SZA)
    sza_geo = 90 - el_deg
    
    # --- 6. Solar azimuth angle (Azimuth) ---
    # North=0, East=90 standard
    # Use atan2 to ensure correct quadrant
    y = -np.cos(delta_rad) * np.sin(H_rad)
    x = np.sin(delta_rad) * np.cos(lat_rad) - np.cos(delta_rad) * np.cos(H_rad) * np.sin(lat_rad)
    
    az_rad = np.arctan2(y, x)
    saa = np.degrees(az_rad)
    if saa < 0: saa += 360
    
    if el_deg > -1: 
        refraction = 1.02 / np.tan(np.radians(el_deg + 10.3 / (el_deg + 5.11)))
        # Temperature and pressure correction
        refraction = refraction * (pressure / 1010.0) * (283.0 / (273.0 + temp))
        # Here refraction is in arc minutes, needs to be converted to degrees
        refraction_deg = refraction / 60.0
    else:
        refraction_deg = 0
        
    el_apparent = el_deg + refraction_deg
    sza_apparent = 90 - el_apparent
    
    return sza_apparent, saa
def utm_to_latlon_approx(easting, northing, zone_number=50, northern=True):
    """
    UTM 转经纬度近似算法 (针对塞罕坝区域优化)
    """
    k0 = 0.9996
    a = 6378137
    f = 1 / 298.257223563

    x = easting - 500000
    y = northing

    m = y / k0
    mu = m / (a * (1 - f / 4 - 3 * f ** 2 / 64 - 5 * f ** 3 / 256))

    e1 = (1 - np.sqrt(1 - f)) / (1 + np.sqrt(1 - f))
    j1 = (3 * e1 / 2 - 27 * e1 ** 3 / 32)
    j2 = (21 * e1 ** 2 / 16 - 55 * e1 ** 4 / 32)
    j3 = (151 * e1 ** 3 / 96)

    fp = mu + j1 * np.sin(2 * mu) + j2 * np.sin(4 * mu) + j3 * np.sin(6 * mu)

    # Calculate Lat/Lon (Simplified)
    # 对于简单的太阳角度计算，这种精度足够
    lat_rad = fp

    c_meridian = (zone_number * 6 - 183)
    lon = c_meridian + np.degrees(x / (a * k0 * np.cos(lat_rad)))
    lat = np.degrees(lat_rad)

    return lat, lon
# ==============================================================================
# Part 2: terrain parameters
# ==============================================================================

class TerrainGeometry:
    def __init__(self, scale=100.0, res=0.5, slope=30.0, aspect=180.0,
                 mode='forest', tree_count=2500.0, margin=0,
                 SZA=30.0, SAA=180.0, FAVD=0.4, 
                 custom_chm=None, custom_trunk=None, tif_path=None,
                 homogeneous_h=15.0, homogeneous_base=5.0):

        self.res = res
        self.slope = slope
        self.aspect = aspect
        self.rows = int(scale / res)
        self.cols = int(scale / res)
        self.margin = margin
        self.mode = mode


        # print(
        #     f"Generating Terrain ({self.rows}x{self.cols}), Mode: {mode}, Margin={self.margin}")

        self.dem, self.grad_x, self.grad_y = self._generate_slope_dem(
            slope, aspect)

        if mode == 'custom':
            # print("Using Custom CHM input.")
            # print(f"Loading CHM from: {tif_path}")
            with tifffile.TiffFile(tif_path) as tif:
                self.chm_top = tif.asarray()
                self.chm_top[self.chm_top < 0] = 0.0  # negetative to zero    
                tags = tif.pages[0].tags
                self.res = tags['ModelPixelScaleTag'].value[0] if 'ModelPixelScaleTag' in tags else 0.1
                self.lat, self.lon = 42.4, 117.2 
                if 'ModelTiepointTag' in tags:
                    tp = tags['ModelTiepointTag'].value
                    if tp[4] > 1000000:
                        self.lat, self.lon = utm_to_latlon_approx(tp[3], tp[4])

        elif mode == 'homogeneous':
            # print(
            #     f"Generating Homogeneous Canopy (H={homogeneous_h}m, Base={homogeneous_base}m).")
            self.chm_top = np.full(
                (self.rows, self.cols), homogeneous_h, dtype=np.float32)
            self.trunk_mask = np.zeros((self.rows, self.cols), dtype=bool)

        else:  # Default: 'forest'
            print(f"Generating Random Forest (Count={tree_count})...")
            self.chm_top, self.chm_bot, self.trunk_mask = self._generate_forest(
                tree_count)
        self.chm_top = np.nan_to_num(self.chm_top, nan=0.0)
        self.rows, self.cols = self.chm_top.shape
        self.roi_mask = binary_fill_holes(self.chm_top > 0) 
        self.chm_bot = self.chm_top * homogeneous_base / homogeneous_h
            
        valid_chm = self.chm_top[self.roi_mask]    
        valid_chm_bot = self.chm_bot[self.roi_mask]
        self.max_h = np.max(valid_chm) if valid_chm.size > 0 else 1.0
        self.mean_h = np.mean(valid_chm - valid_chm_bot) if valid_chm.size > 0 else 0.0
        self.veg_cover = np.sum(valid_chm > 0) / valid_chm.size if valid_chm.size > 0 else 0.0            
            #  (Continuity Index)
        self.c_index = np.mean(valid_chm) / self.max_h if self.max_h > 0 else 0
        # self.c_index = 0
        self.dem, self.grad_x, self.grad_y = self._generate_slope_dem(slope, aspect)
        self.svf = (1 + cosd(slope)) / 2.0
        self.scene_max_z = np.max(self.dem + self.chm_top) + 10.0
        self.hemi_paths = None
        self.path_normal, self.normal_vector = self.get_slope_normal_params()
        self.h_normal, self.normal = self.path_normal, self.normal_vector
        self.sun_mask = self.calc_sun_params(SZA, SAA)
        self.path_s = self.path_angles(SZA, SAA)
        self.FAVD = FAVD
        self.LAI = self.mean_h * FAVD
        # print(f"Terrain LAI: {self.LAI:.2f}, FAVD: {FAVD:.2f}")
        self.iD_val = self.precompute_diffuse_iD(FAVD)
    def _generate_slope_dem(self, slope_deg, aspect_deg):
        y, x = np.mgrid[0:self.rows, 0:self.cols]
        y, x = y * self.res, x * self.res
        s_rad, a_rad = np.radians(slope_deg), np.radians(aspect_deg)

        # Calculate terrain gradient (m/m)
        grad_x = -np.tan(s_rad) * np.sin(a_rad)
        grad_y = np.tan(s_rad) * np.cos(a_rad)

        base = max(self.rows, self.cols) * self.res * np.tan(s_rad) + 10
        dem = base + (grad_x * x + grad_y * y)

        # Return dem array and gradient constants
        return dem, grad_x, grad_y

    def _generate_forest(self, count, max_h=15, min_h=5, crown_r=3):
        # ... (keep unchanged) ...
        chm_top = np.zeros((self.rows, self.cols), dtype=np.float32)
        chm_bot = np.zeros_like(chm_top)
        trunk_mask = np.zeros((self.rows, self.cols), dtype=bool)

        np.random.seed(42)
        rows_idx = np.random.randint(0, self.rows, count)
        cols_idx = np.random.randint(0, self.cols, count)
        heights = np.random.uniform(min_h, max_h, count)
        r_px = crown_r / self.res

        y_g, x_g = np.ogrid[:self.rows, :self.cols]
        sorted_idx = np.argsort(heights)

        for i in sorted_idx:
            r, c, h = rows_idx[i], cols_idx[i], heights[i]
            h_base = h * 0.33
            crown_h = h - h_base
            dist_sq = (y_g - r)**2 + (x_g - c)**2
            norm_dist_sq = dist_sq / (r_px**2 + 1e-6)
            in_radius = dist_sq <= r_px**2
            ellipsoid_z = np.zeros_like(dist_sq, dtype=np.float32)
            valid_mask = in_radius & (norm_dist_sq <= 1)

            if np.any(valid_mask):
                ellipsoid_z[valid_mask] = h_base + crown_h * \
                    np.sqrt(1.0 - norm_dist_sq[valid_mask])

            update = valid_mask & (ellipsoid_z > chm_top)
            chm_top[update] = ellipsoid_z[update]
            chm_bot[update] = h_base
            trunk_r = max(1, r_px * 0.15)
            trunk_mask[update] = (dist_sq[update] <= trunk_r**2)

        return chm_top, chm_bot, trunk_mask

    def get_slope_normal_params(self):
        nx = sind(self.slope) * sind(self.aspect)
        ny = -sind(self.slope) * cosd(self.aspect)
        nz = cosd(self.slope)
        paths_normal = self.path_angles(self.slope, self.aspect)
        paths_normal_nogaps = np.mean(paths_normal[paths_normal > 0])
        return paths_normal_nogaps, (nx, ny, nz)

    @staticmethod
    def calc_cos_local(zenith, azimuth, normal_vec):
        nx, ny, nz = normal_vec
        z_rad, a_rad = np.radians(zenith), np.radians(azimuth)
        lx = np.sin(z_rad) * np.sin(a_rad)
        ly = -np.sin(z_rad) * np.cos(a_rad)
        lz = np.cos(z_rad)
        cos_loc = lx*nx + ly*ny + lz*nz
        return max(0.001, cos_loc)

    def _get_step_params(self, azimuth, zenith):
        az_rad = np.radians(azimuth)
        zen_rad = np.radians(zenith)
        if zen_rad < 1e-4:
            return 0.0, 0.0, self.res, self.res
        else:
            vec_x = np.sin(az_rad)
            vec_y = -np.cos(az_rad)
            vec_z = 1.0 / np.tan(zen_rad)
            scale = 1.0 / max(abs(vec_x), abs(vec_y))
            dx, dy, dz = vec_x*scale, vec_y*scale, vec_z*scale*self.res
            step_len = np.sqrt((dx*self.res)**2 + (dy*self.res)**2 + dz**2)

            # --- 优化：限制步长 ---
            # 原来的限制只防止过大 (2.0 * res)
            # 现在为了适应低矮冠层，当 step_len 导致 dz 超过 0.5m 时，强制缩小步长
            # 这样保证在薄冠层（比如1m厚）里至少采样2次
            max_dz = 0.5  # 允许的最大垂直步长 (米)

            # 检查当前 dz 大小
            current_dz = dz
            if current_dz > max_dz:
                factor = max_dz / current_dz
                dx *= factor
                dy *= factor
                dz *= factor
                step_len *= factor

            # 保留原有的基于 res 的限制 (防止水平跨度太大)
            if step_len > 2.0 * self.res:
                factor = 2.0 * self.res / step_len
                dx *= factor
                dy *= factor
                dz *= factor
                step_len = 2.0 * self.res

            return dx, dy, dz, step_len

    def trace_surface_mask(self, azimuth, zenith):
        ground_height = 1
        """Used to determine if surface is illuminated (binary)"""
        if zenith >= 90: 
            return np.zeros_like(self.dem, dtype=bool)

        if zenith < 1e-3: 
            return self.chm_top < ground_height

        dx, dy, dz, _ = self._get_step_params(azimuth, zenith)
        

        max_search_height = np.nanmax(self.chm_top)
        max_steps = int(max_search_height / dz) + 10 

        max_steps = min(max_steps, 500)

        curr_z = self.dem + ground_height
        
        lit_mask = np.ones_like(self.dem, dtype=bool)
        
        if self.mode == 'homogeneous':
            return lit_mask
        else:

            pad_val = max_search_height 
            
            for k in range(1, max_steps):
                sx, sy = dx * k, dy * k 
                curr_z += dz
                
                shift_res_x, shift_res_y = sx * self.res, sy * self.res
                delta_z_terrain = (self.grad_x * shift_res_x) + (self.grad_y * shift_res_y)
                s_dem = self.dem + delta_z_terrain
                
                s_top = shift(self.chm_top, (-sy, -sx), order=0, cval = pad_val)
                s_bot = shift(self.chm_bot, (-sy, -sx), order=0, cval = 0)
                
                hit_crown = (s_top > 0) & (curr_z > (s_dem + s_bot)) & (curr_z < (s_dem + s_top - 1))
                
                # Combined occlusion
                is_blocked = hit_crown
                
                # Update mask: Once blocked, permanently False
                lit_mask[is_blocked] = False
                
                if not np.any(lit_mask):
                    break
                    
            return lit_mask

    def precompute_diffuse_iD(self, FAVD, leaf_class=6, ALA=57.5):
        za = np.arange(0, 86, 5) 
        za = np.append(za, 89)
        iv = np.ones_like(za, dtype=float)
        for i in range(len(iv)):
            G = Ross_G_function(za[i], 0, leaf_class, ALA)
            pg_corr = calc_gap_prob(self.path_angles_normal(za[i]), FAVD, G)
            iv[i] = (1-pg_corr) * sind(2*za[i])
        iD = np.trapz(iv, np.deg2rad(za))
        return iD

    def calc_sun_params(self, sza, saa):
        return self.trace_surface_mask(saa, sza)
    
    def get_fast_geometry(self, sun_params, sza, saa, vza, vaa):
        if self.mode == 'homogeneous':
            geo_comps = (1, 0, 0, 0, 0)
            return geo_comps
        else:
            sun_mask = sun_params
            view_mask = self.trace_surface_mask(vaa, vza)
            
            m = self.margin
            roi = (slice(m, -m), slice(m, -m)
                ) if m > 0 else (slice(None), slice(None))
            
            sm_roi = sun_mask[roi]
            vm_roi = view_mask[roi]
            veg_roi = (self.chm_top[roi] > 0)

            total_pixels = sm_roi.size
            if total_pixels == 0:
                return (0, 0, 0, 0), np.array([0.]), np.array([0.])

            kgz_s = np.sum(sun_mask[roi]) / total_pixels
            kgz_v = np.sum(view_mask[roi]) / total_pixels
            
            # --- Calculate four components ---
            # Kg (Sunlit Ground): 
            kg = np.sum(view_mask[roi] & sun_mask[roi]) / total_pixels
            
            # Kc (Sunlit Vegetation):

            kct = 1.0 - kgz_v
        
            delta_deg = cosd(sza)*cosd(vza) + sind(sza)*sind(vza)*cosd(vaa-saa)
            phi = acosd(delta_deg)
            
            
            c = self.c_index
            # print(f"Canopy Continuity Index c: {c:.3f}")
            delta_mod = cosd(phi * (1 - np.sin(np.pi * c / 2.0)))
        
            # Kc 
            kc = 0.5 * (1 + delta_mod) * kct
                
            # Kt (Shaded Vegetation) = Total Visible Veg - Sunlit Veg
            kt = kct - kc
            if kt < 0: kt = 0 # 
            
            # Kz (Shaded Ground) = Total Visible Ground - Sunlit Ground
            kz = kgz_v - kg
            if kz < 0: kz = 0
            
            geo_comps = (kc, kt, kg, kz, kgz_s)
            
            # print(f"  P_gap Sun: {kgz_s:.3f}, View: {kgz_v:.3f}, View zenith: {vza}, RAA : {vaa-saa}, KC: {kc:.3f}, KT: {kt:.3f}, KG: {kg:.3f}, KZ: {kz:.3f}")

            return geo_comps
    
    def path_angles(self, vza, vaa):
        denom = cosd(vza) * (1 + tand(self.slope) * cosd(vaa - self.aspect) * np.abs(tand(vza)))
        path_v = (self.chm_top - self.chm_bot) / denom
        path_nogap = path_v[path_v > 0]
        return path_nogap   
    def path_angles_normal(self, vza):
        denom = cosd(vza)
        path_v = (self.path_normal) / denom
        path_nogap = path_v[path_v > 0]
        return path_nogap     

# ==============================================================================
# Part 3: path_RT_T
# ==============================================================================


def PATH_RT_Terrain(terrain, tau_l, rho_l, soil_r,
                    geo_comps, sky_ratio,
                    sza, saa, vza, vaa, branchFactor = 0.9,
                    leaf_class= 6, ALA=57.5, Hotspot=0.02):
    
    LAI = terrain.LAI
    FAVD = terrain.FAVD
    iD_value = terrain.iD_val
    path_s = terrain.path_s
    path_v = terrain.path_angles(vza, vaa)
    SLOPE = terrain.slope
    h_eff_slope = terrain.h_normal
    normal_vec = terrain.normal_vector
    svf = terrain.svf
    (kc, kt, kg, kz, kgz_s) = geo_comps
    omega = rho_l + tau_l   
    Gs = Ross_G_function(sza, saa, leaf_class, ALA)
    Gv = Ross_G_function(vza, vaa, leaf_class, ALA)
    # 1. Geometry
    cos_s_local = TerrainGeometry.calc_cos_local(sza, saa, normal_vec)
    cos_v_local = TerrainGeometry.calc_cos_local(vza, vaa, normal_vec)
    # 2. Gap Probabilities
    ps_dir_all = calc_gap_prob(path_s, FAVD, Gs) * (1 - kgz_s) + kgz_s
    pv_dir_all = calc_gap_prob(path_v, FAVD, Gv) * (1 - kg - kz) + kg + kz

    paths_s_withincrown = np.mean(path_s) 
    paths_v_withincrown = np.mean(path_v)

    # paths_v_withincrown = min(paths_v_withincrown, 1000)

    # inv_ls_loc = 1.0 / cos_s_local
    # inv_lv_loc = 1.0 / cos_v_local
    inv_ls_loc = paths_s_withincrown / h_eff_slope
    inv_lv_loc = paths_v_withincrown / h_eff_slope
    # 3. Diffuse Light
    iD = iD_value

    # 4. Veg Single Scattering
    gamma_val = area_scatter_phase_function(
        rho_l, tau_l, sza, saa, vza, vaa, leaf_class, ALA)

    H_eff = h_eff_slope
    z_steps = 100

    i_arr = np.arange(1, z_steps+1)
    z_arr = (H_eff / z_steps) * (i_arr)

    pz_s_arr = np.exp(-Gs * FAVD * z_arr * inv_ls_loc)
    pz_v_arr = np.exp(-Gv * FAVD * z_arr * inv_lv_loc)

    chs_arr = calc_hotspot_factor(Hotspot, sza, saa, vza, vaa, Gs, Gv,
                                  FAVD, z_arr, inv_ls_loc, inv_lv_loc)

    p_bi_arr = pz_s_arr * pz_v_arr * chs_arr
    p_bi_kt_arr = pz_s_arr * pz_v_arr

    int_sun = np.trapz(p_bi_arr, z_arr)
    # int_sun = np.sum(p_bi_arr) * (H_eff / z_steps)
    int_shade = np.sum(p_bi_kt_arr) * (H_eff / z_steps)
    # brf_veg_sun = kc * (gamma_val * int_sun * FAVD / cosd(sza)) / cosd(vza)
    # brf_veg_shade = kt * (gamma_val * int_shade * FAVD * 0.1 / cosd(sza))/ cosd(vza)
    brf_veg_sun = kc * (gamma_val * int_sun * FAVD / cosd(sza)) * inv_lv_loc
    brf_veg_shade = kt * (gamma_val * int_shade * FAVD *
                          0.1 / cosd(sza))*inv_lv_loc
    # 5. Soil BRF
    ps_dir = pz_s_arr[-1]
    pv_dir = pz_v_arr[-1]
    chs = chs_arr[-1]
    intensity_ratio = cos_s_local / cosd(sza)
    brf_soil = ((kg + (kc + kt) * (ps_dir * pv_dir * chs)) *
                soil_r + (kz * ps_dir * soil_r)) * intensity_ratio

    # 6. Multiple Scattering
    i_0, i_v = 1 - ps_dir_all, 1 - pv_dir_all
    
    i_0 = sky_ratio * iD * svf + (1 - sky_ratio) * i_0
    
    LAI = LAI * cosd(SLOPE)
    denom = 1 - (1 - iD/(LAI+1e-6)) * omega + 1e-6
    p_recol = 1 - iD/(LAI+1e-6)
    esc_v = i_v / (2*LAI+1e-6)
    esc_h = iD / (2*LAI+1e-6)

    BRF_vm = i_0 * (omega**2) * p_recol * esc_v / denom 
    Tdn = 1 - i_0 + i_0 * omega * esc_h / denom
    Tup = 1 - i_v + iD * omega * esc_v / denom
    Rdn = iD * omega * esc_h / denom
    BRF_vs = (soil_r / (1 - soil_r * Rdn)) * Tdn * \
        Tup - (1 - i_0) * soil_r * (1 - i_v) 
        

        
    if cos_v_local < 0.09:
        return np.nan
    else:
        return (1 - sky_ratio) * (brf_soil + brf_veg_sun + brf_veg_shade) + (BRF_vm + BRF_vs) * \
            cos_s_local/cosd(sza) * branchFactor

# ==============================================================================
# Main
# ==============================================================================

if __name__ == "__main__":
    print("Initializing...")
