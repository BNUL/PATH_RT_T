import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import shift, binary_fill_holes, label
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
    计算叶片角度分布函数 g_L(theta_L)

    参数:
        iorien: 叶片法向分布类型 (1~6)
            1 - planophile     (a=1,  b=2)
            2 - erectophile    (a=-1, b=-2)
            3 - plagiophile    (a=-1, b=4)
            4 - extremophile   (a=1,  b=4)
            5 - uniform        (a=0,  b=任意)
            6 - spherical      g(θ) = sin(θ)
        theta_L: 叶倾角 (弧度)，可以是标量或 numpy 数组
    返回:
        g_L: 对应的叶片角度分布函数值
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
            raise ValueError("iorien 必须是 1~6 之间的整数")

        g_L = (2 / np.pi) * (1 + a * np.cos(b * theta_L))

    return g_L

def get_APF(a, b, rol, taul, cteta, phi, cteta1, phi1):
    """
    计算双向反射/透射相关的面积投影因子 (APF)
    数值积分方法与原 MATLAB 代码一致

    参数:
        a, b         : 叶片角度分布参数 (用于 1 + a*cos(b*theta_L))
        rol          : 反射率 (reflection coefficient)
        taul         : 透射率 (transmission coefficient)
        cteta, phi   : 入射方向的天顶角余弦 & 方位角 (度)
        cteta1, phi1 : 观测/散射方向的天顶角余弦 & 方位角 (度)

    返回:
        get_APF      : 计算得到的面积投影因子值
    """
    # 预处理角度
    cteta1 = cteta1 + 180          # 重要：180度反向处理

    ccteta  = np.cos(np.deg2rad(cteta))     # cos(θ)
    ccteta1 = np.cos(np.deg2rad(cteta1))    # cos(θ₁)

    steta   = np.sin(np.deg2rad(cteta))     # sin(θ)
    steta1  = np.sin(np.deg2rad(cteta1))    # sin(θ₁)

    phi_rad  = np.deg2rad(phi)
    phi1_rad = np.deg2rad(phi1)

    kpi3 = 1.0 / (np.pi ** 3)

    # 积分网格（与原代码一致）
    n = 30
    m = 4 * n
    h_theta = 0.5 * np.pi / n
    h_fi    = 2.0 * np.pi / m

    theta_i = 0.5 * h_theta      # 从中点开始
    fi_1    = 0.5 * h_fi

    integral = 0.0

    for i in range(n):
        fi_j = fi_1
        xx = 0.0
        c_i = np.cos(theta_i)
        s_i = np.sin(theta_i)

        for j in range(m):
            # 两个方向与叶片法线的夹角余弦
            yy = ccteta  * c_i + steta  * s_i * np.cos(phi_rad  - fi_j)
            zz = ccteta1 * c_i + steta1 * s_i * np.cos(phi1_rad - fi_j)

            # 乘积
            zz *= yy

            # 根据正负号选择反射或透射
            if zz <= 0:
                xx += rol * abs(zz)
            else:
                xx += taul * zz

            fi_j += h_fi

        xx *= h_fi
        yy = 1.0 + a * np.cos(b * theta_i)     # 注意：这里是简化形式
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
    计算积分 G(Omega) = (1/(2*pi)) ∫∫ g_L(θ_L) |cos(Ω, Ω_L)| dθ_L dφ_L
    参数:
        iorien: 叶片法向分布类型 (1~6)，同 get_gFun
        theta: 太阳/观测方向的天顶角 (度)
        fi: 太阳/观测方向的方位角 (度)
    返回:
        G_Fun: 投影函数 G(Omega)
    注意：使用 Simpson 规则或梯形积分的简单数值积分（原代码使用矩形积分）
    """
    # 转换为弧度
    theta = np.deg2rad(theta)
    fi = np.deg2rad(fi)

    # 积分网格设置（与原 MATLAB 一致）
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
        g_val = get_gFun(iorien, theta_i)  # 调用前面定义的 get_gFun
        G_Fun += g_val * xx
        theta_i += h_theta

    G_Fun *= h_theta
    G_Fun /= (2 * np.pi)

    return G_Fun
def area_scatter_phase_function(rho_l, tau_l, theta_0, phi_0, theta, phi, leaf_class, ALA=57.5):
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
        b = 0
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
    计算 Kuusk 模型中的热点因子 (Hot Spot Factor) Chs

    参数:
        par   : 叶片平均尺寸 (bl, leaf size parameter)
        SZA   : 太阳天顶角 (度)
        SAA   : 太阳方位角 (度)
        VZA   : 观测天顶角 (度)
        VAA   : 观测方位角 (度)
        Gs    : 太阳方向的 G 函数值 (投影函数)
        Gv    : 观测方向的 G 函数值
        uL    : 平均叶面积密度 (LAI / 冠层高度) 或 FAVD
        z     : 相对高度（从冠层顶部向下，单位与冠层总高度一致）
        ls    : 太阳方向单位冠层深度相对光程 (mu0)
        lv    : 观测方向单位冠层深度相对光程 (muv)

    返回:
        Chs   : 在高度 z 处的热点因子值
    """
    # 角度转弧度
    SZA = np.deg2rad(SZA)
    SAA = np.deg2rad(SAA)
    VZA = np.deg2rad(VZA)
    VAA = np.deg2rad(VAA)

    mu0 = 1/ls
    muv = 1/lv

    # f1 部分
    f1 = np.sqrt(Gs * Gv / (mu0 * muv))

    # 计算相对方位角 phi 并规范化到 [0, 360)
    phi = VAA - SAA
    if phi < 0:
        phi += 2 * np.pi
    if phi > 2 * np.pi:
        phi -= 2 * np.pi

    # 计算 cos(gamma) - 两个方向的夹角余弦
    cosgamma = (np.cos(SZA) * np.cos(VZA) +
                np.sin(SZA) * np.sin(VZA) * np.cos(VAA - SAA))
    cosgamma = np.clip(cosgamma, -1.0, 1.0)  # guard against fp > 1 at exact hotspot

    # 计算 delta（Kuusk 热点模型中的关键参数）
    delta = np.sqrt(np.maximum(0.0, 1/np.cos(SZA)**2 + 1/np.cos(VZA)**2 -
                    2 * cosgamma / (np.cos(SZA) * np.cos(VZA))))

    # 防止除以零或极小值
    delta = np.maximum(delta, 1e-5)   # use np.maximum (not max) to handle NaN-safety

    # 计算 f2
    f2 = (uL * par / delta) * (1 - np.exp(-z * delta / par))

    # 最终热点因子
    Chs = np.exp(f1 * f2)

    return Chs
def _as_float_array(value: float | np.ndarray) -> np.ndarray:
    return np.asarray(value, dtype=float)
def geo_transform(
    s_zenith: float | np.ndarray,
    a_azimuth: float | np.ndarray,
    alpha: float | np.ndarray,
    beta: float | np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    s_zenith = _as_float_array(s_zenith)
    a_azimuth = _as_float_array(a_azimuth)
    alpha = _as_float_array(alpha)
    beta = _as_float_array(beta)

    x_value = np.sin(s_zenith) * np.cos(alpha) * np.cos(a_azimuth - beta) - np.cos(s_zenith) * np.sin(alpha)
    y_value = np.sin(a_azimuth - beta) * np.sin(s_zenith)
    z_value = np.sin(s_zenith) * np.sin(alpha) * np.cos(a_azimuth - beta) + np.cos(s_zenith) * np.cos(alpha)
    z_value = np.clip(z_value, -1.0, 1.0)
    t_zenith = np.arccos(z_value)
    ta_azimuth = np.mod(np.arctan2(y_value, x_value), 2.0 * np.pi)
    return t_zenith, ta_azimuth

def calc_gap_prob(path_array, favd, G):
    if path_array.size == 0:
        return 1.0
    return np.mean(np.exp(-G * favd * path_array))

def calc_sun_position(lat, lon, dt, timezone=8, pressure=1010, temp=10):

    # --- 1. 时间处理 (统一转为 UTC) ---
    if dt.tzinfo is None:
        # 如果是 Naive Time，减去时区得到 UTC
        dt_utc = dt - timedelta(hours=timezone)
    else:
        # 如果是 Aware Time，直接转 UTC
        dt_utc = dt.astimezone(timezone.utc)

    # 计算儒略日 (Julian Day)
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
    t = (jd - 2451545.0) / 36525.0  # 儒略世纪数 (J2000起)

    # 太阳平黄经 (L0)
    L0 = 280.46646 + 36000.76983 * t + 0.0003032 * t**2
    L0 = L0 % 360

    # 太阳平近点角 (M)
    M = 357.52911 + 35999.05029 * t - 0.0001537 * t**2
    M_rad = np.radians(M)

    # 轨道离心率 (e)
    e = 0.016708634 - 0.000042037 * t - 0.0000001267 * t**2

    # 太阳方程中心差 (C)
    C = (1.914602 - 0.004817 * t - 0.000014 * t**2) * np.sin(M_rad) + \
        (0.019993 - 0.000101 * t) * np.sin(2 * M_rad) + \
        0.000289 * np.sin(3 * M_rad)

    # 太阳真黄经 (True Longitude, theta)
    theta = L0 + C

    omega = 125.04 - 1934.136 * t
    lambda_sun = theta - 0.00569 - 0.00478 * np.sin(np.radians(omega))
    lambda_rad = np.radians(lambda_sun)

    # 平均黄赤交角 (Mean Obliquity, epsilon0)
    epsilon0 = 23 + 26/60.0 + 21.448/3600.0 - 46.8150/3600.0 * t - \
               0.00059/3600.0 * t**2 + 0.001813/3600.0 * t**3

    epsilon = epsilon0 + 0.00256 * np.cos(np.radians(omega))
    eps_rad = np.radians(epsilon)


    alpha_rad = np.arctan2(np.cos(eps_rad) * np.sin(lambda_rad), np.cos(lambda_rad))
    # 赤纬 (delta)
    delta_rad = np.arcsin(np.sin(eps_rad) * np.sin(lambda_rad))


    GMST0 = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + \
            0.000387933 * t**2 - (t**3) / 38710000
    GMST0 = GMST0 % 360
    # 视恒星时需要加上章动修正 (此处省略微小章动项，误差<0.001度，如需极高精度可补)
    GAST = GMST0

    # --- 5. 转换到地平坐标系 (Az, El) ---
    # 地方时角 (H)
    H = (GAST + lon) - np.degrees(alpha_rad)
    H = H % 360
    if H > 180: H -= 360
    H_rad = np.radians(H)

    lat_rad = np.radians(lat)

    # 高度角 (Elevation)
    sin_el = np.sin(lat_rad) * np.sin(delta_rad) + np.cos(lat_rad) * np.cos(delta_rad) * np.cos(H_rad)
    el_rad = np.arcsin(np.clip(sin_el, -1, 1))
    el_deg = np.degrees(el_rad)

    # 几何天顶角 (Geometric SZA)
    sza_geo = 90 - el_deg

    # --- 6. 太阳方位角 (Azimuth) ---
    # 北=0, 东=90 标准
    # 使用 atan2 确保象限正确
    y = -np.cos(delta_rad) * np.sin(H_rad)
    x = np.sin(delta_rad) * np.cos(lat_rad) - np.cos(delta_rad) * np.cos(H_rad) * np.sin(lat_rad)

    az_rad = np.arctan2(y, x)
    saa = np.degrees(az_rad)
    if saa < 0: saa += 360

    if el_deg > -1:
        refraction = 1.02 / np.tan(np.radians(el_deg + 10.3 / (el_deg + 5.11)))
        # 温度气压修正
        refraction = refraction * (pressure / 1010.0) * (283.0 / (273.0 + temp))
        # 这里的 refraction 是角分，需要转为度
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
                 SZA=30.0, SAA=180.0, FAVD=0.4, base_ratio =0.3,
                 custom_chm=None, custom_trunk=None, tif_path=None,
                 homogeneous_h=15.0, homogeneous_base=5.0,
                 veg_cover=None, brratio=None, crown_base_frac=None):

        self.res = res
        self.slope = slope
        self.aspect = aspect
        self.rows = int(scale / res)
        self.cols = int(scale / res)
        self.margin = margin
        self.mode = mode
        self.base_ratio = base_ratio


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
        # 枝下高计算：
        #   homogeneous模式：严格按输入的 homogeneous_h / homogeneous_base 计算，不使用额外高度阈值
        #   custom/forest模式：沿用经验比例
        if mode == 'homogeneous':
            base_ratio = homogeneous_base / homogeneous_h if homogeneous_h > 0 else 0.0
            self.chm_bot = np.clip(self.chm_top * base_ratio, 0.0, self.chm_top)
        else:
            self.chm_bot = np.where(self.chm_top < 1.0, 0.0, self.chm_top * self.base_ratio)

        valid_chm = self.chm_top[self.roi_mask]
        valid_chm_bot = self.chm_bot[self.roi_mask]
        self.max_h = np.max(valid_chm) if valid_chm.size > 0 else 1.0
        self.mean_h = np.mean(valid_chm - valid_chm_bot) if valid_chm.size > 0 else 0.0

        chm_cover = np.sum(valid_chm > 0.5) / valid_chm.size if valid_chm.size > 0 else 0.
        self.veg_cover = veg_cover if veg_cover is not None else chm_cover
        self.crownLength = self.mean_h / self.veg_cover if self.veg_cover > 0 else 0
        self.brratio = max(1, min(3, brratio if brratio is not None else self._estimate_brratio(crown_base_frac=crown_base_frac)))

        self.std_h = np.std(valid_chm)
        base_c = (self.mean_h - self.std_h) / self.max_h
        # self.c_index = np.mean(valid_chm) / self.max_h if self.max_h > 0 else 0
        self.c_index = base_c ** 2
        self.svf = (1 + cosd(slope)) / 2.0
        self.path_normal, self.normal_vector = self.get_slope_normal_params()
        self.h_normal, self.normal = self.path_normal, self.normal_vector
        self.path_s = self.path_angles(SZA, SAA)
        self.FAVD = FAVD
        self.LAI = self.mean_h * FAVD
        # print(f"Terrain LAI: {self.LAI:.2f}, FAVD: {FAVD:.2f}")
        self.iD_val = self.precompute_diffuse_iD(FAVD)

    def _estimate_brratio(self, crown_base_frac=None, fallback=2.0):
        """
        Estimate crown b/r from CHM labeled-patch statistics.

        Primary estimator — variance method (Li & Strahler 1992 spheroid geometry):
            For a prolate spheroid with vertical semi-axis b, the area-weighted
            variance of CHM heights over the projected crown circle satisfies
                sigma_h^2 = b^2 / 18  =>  b = 3*sqrt(2) * sigma_h
            r is the equivalent-circle radius: r = sqrt(N * res^2 / pi).

        Connection to crown-base fraction f (枝下高比例, H_base/H_tree):
            b = h_max*(1-f)/2  and  sigma_h/h_max = (1-f) / (6*sqrt(2))
            So  f  can be read from CHM stats, or used as a prior to constrain b.

        Fallback order when CHM has no height variation (flat-top synthetic):
            1. crown_base_frac (e.g. 1/3): b = h_max*(1-f)/2
            2. constant `fallback`

        # --- Original range estimator (kept for reference) ---
        # b_i = 3.0 * (float(h.max()) - float(h.mean()))   # b = 3*(h_max - h_mean)
        """
        _lbl = label(self.chm_top > 0.5)
        assert isinstance(_lbl, tuple)
        labeled, n = _lbl[0], int(_lbl[1])
        if n == 0:
            return fallback
        r_list, b_list = [], []
        _sqrt18 = float(np.sqrt(18.0))
        for i in range(1, n + 1):
            px = labeled == i
            n_px = int(np.sum(px))
            if n_px < 4:
                continue
            r_i = np.sqrt(n_px * self.res ** 2 / np.pi)
            h = self.chm_top[px]
            h_std_i  = float(h.std())
            h_max_i  = float(h.max())

            b_var   = _sqrt18 * h_std_i                          # variance estimator
            b_ratio = h_max_i * (1.0 - crown_base_frac) / 2.0 \
                      if crown_base_frac is not None else None    # crown-base prior

            # prefer variance estimator; fall back to crown-base prior for flat CHMs
            if b_var > 0.2 * r_i:
                b_i = b_var
            elif b_ratio is not None:
                b_i = b_ratio
            else:
                continue
            r_list.append(r_i)
            b_list.append(b_i)
        if not r_list:
            return fallback
        med_r = float(np.median(r_list))
        med_b = float(np.median(b_list))
        if med_r < 0.1 or med_b < 0.1 * med_r:
            return fallback
        return med_b / med_r

    def _generate_forest(self, count, max_h=15, min_h=5, crown_r=3):
        # ... (保持不变) ...
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


    def precompute_diffuse_iD(self, FAVD, leaf_class=6, ALA=57.5):
        za        = np.arange(0, 86, 5)
        za        = np.append(za, 89)
        azimuths  = np.arange(0, 360, 45)
        iv        = np.ones_like(za, dtype=float)
        for i in range(len(iv)):
            G       = Ross_G_function(za[i], 0, leaf_class, ALA)
            # P_inter = float(np.mean([self.gap_fraction_dir(float(az), float(za[i]))
            #                          for az in azimuths]))
            P_intra = calc_gap_prob(self.path_angles_normal(za[i]), FAVD, G)
            iv[i]   = (1 - P_intra) * sind(2 * za[i])
        iD = np.trapz(iv, np.deg2rad(za))
        return iD

    def get_fast_geometry(self, sza, saa, vza, vaa):

        if self.mode == 'homogeneous':
            return (1, 0, 0, 0, 0)

        # --- 树冠结构参数 ---
        # crownradius = 2    # RAMI HET09 horizontal semi-axis (m)
        # crown_b    = 4.0   # RAMI HET09 vertical semi-axis: crown depth 8 m → b = 4 m
        # crownLength = mean_h / veg_cover gives ≈ tree height, not crown depth, because
        # mean_h is averaged over crown pixels only (not area-averaged). Use correct b/r:
        brratio = self.brratio  # = 2.0
        # 树干高度（冠层底部高度均值），与 Li-Strahler 模型中的 htdif 对应
        veg_bot  = self.chm_bot[self.chm_top > 0.5]
        htdif    = float(np.mean(veg_bot)) if veg_bot.size > 0 else 5.0
        tsp = np.deg2rad(sza)
        phsp = np.deg2rad(saa)
        tvp = np.deg2rad(vza)
        phvp = np.deg2rad(vaa)
        slope = np.deg2rad(self.slope)
        aspect = np.deg2rad(self.aspect)
        tvp0, phvp0 = geo_transform(0, 0, slope, aspect)
        lamda = -np.log(1 - self.veg_cover) / (np.pi * np.cos(tvp0))

        slope = np.arctan(np.tan(slope) / brratio)
        tsp = np.arctan(np.maximum(np.tan(tsp) * brratio, 0.0))
        tsp, phsp = geo_transform(tsp, phsp, slope, aspect)
        tvp = np.arctan(np.maximum(np.tan(tvp) * brratio, 0.0))
        tvp, phvp = geo_transform(tvp, phvp, slope, aspect)

        gammai = np.pi / np.cos(tsp)
        gammav = np.pi / np.cos(tvp)
        kgv = np.exp(-gammav * lamda)
        kgs = np.exp(-gammai * lamda)

        f1=np.sqrt(kgs*kgv*(1-kgs)*(1-kgv))
        cosgamma= np.cos(tsp)*np.cos(tvp) + np.sin(tsp)*np.sin(tvp)*np.cos(phvp-phsp)
        cosgamma = np.clip(cosgamma, -1.0, 1.0)  # guard against fp > 1 at exact hotspot
        delta=np.sqrt(np.maximum(0.0,
                      1/(np.cos(tsp)**2)+1/(np.cos(tvp)**2)
                      -2*cosgamma/(np.cos(tsp)*np.cos(tvp))))
        # f2= np.exp(-delta/crownradius/2*self.crownLength)
        f2= np.exp(-delta*brratio)
        kg = kgv * kgs + f1*f2
        kz = kgv - kg
        fcv = 1  - kgv
        c = self.c_index
        phi = acosd(cosgamma)
        delta_mod = cosd(phi * (1 - np.sin(np.pi * c / 2.0)))

        kc = 0.5 * (1 + delta_mod) * fcv

        kt = fcv - kc
        return (kc, kt, kg, kz, kgs)

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
                    sza, saa, vza, vaa, branchFactor = 1,
                    leaf_class= 6, ALA=57.5, Hotspot=0.02):

    LAI = terrain.LAI
    FAVD = terrain.FAVD
    iD = terrain.iD_val
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


    # inv_ls_loc = 1.0 / cos_s_local
    # inv_lv_loc = 1.0 / cos_v_local
    inv_ls_loc = paths_s_withincrown / h_eff_slope
    inv_lv_loc = paths_v_withincrown / h_eff_slope
    # 3. Diffuse Light

    # 4. Veg Single Scattering
    gamma_val = area_scatter_phase_function(
        rho_l, tau_l, sza, saa, vza, vaa, leaf_class, ALA)

    z_steps = 100

    i_arr = np.arange(1, z_steps+1)
    z_arr = (h_eff_slope / z_steps) * (i_arr)

    pz_s_arr = np.exp(-Gs * FAVD * z_arr * inv_ls_loc)
    pz_v_arr = np.exp(-Gv * FAVD * z_arr * inv_lv_loc)

    ps_dir = pz_s_arr[-1]
    pv_dir = pz_v_arr[-1]

    chs_arr = calc_hotspot_factor(Hotspot, sza, saa, vza, vaa, Gs, Gv,
                                  FAVD, z_arr, inv_ls_loc, inv_lv_loc)

    p_bi_arr = pz_s_arr * pz_v_arr * chs_arr
    p_bi_kt_arr = pz_s_arr * pz_v_arr

    int_sun = np.trapz(p_bi_arr, z_arr)
    # int_sun = np.sum(p_bi_arr) * (H_eff / z_steps)
    int_shade = np.sum(p_bi_kt_arr) * (h_eff_slope / z_steps)
    # brf_veg_sun = kc * (gamma_val * int_sun * FAVD / cosd(sza)) / cosd(vza)
    # brf_veg_shade = kt * (gamma_val * int_shade * FAVD * 0.2 / cosd(sza))/ cosd(vza)
    brf_veg_sun = kc * (gamma_val * int_sun * FAVD / cosd(sza)) * inv_lv_loc
    brf_veg_shade = kt * (gamma_val * int_shade * FAVD *
                          np.sqrt(ps_dir) / cosd(sza))*inv_lv_loc
    # 5. Soil BRF

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
