"""
BRDF 矩阵分析脚本
======================
绘制多维度BRDF折线图
- 固定条件：坡度30°, SZA=30°, SAA=0°, 波段=NIR
- 4个坡向：0°, 90°, 180°, 270°
- 3行分析：
  行1: LAI变化 (0.5, 1, 2, 4, 6)
  行2: 叶倾角/LeafClass (6种)
  行3: 天空散射比 (0, 0.2, 0.4, 0.6, 0.8, 1.0)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

# 导入模型
import PATH_RT_T_v6 as MyModel

# 设置全局字体
plt.rcParams['font.family'] = 'arial'
plt.rcParams['font.size'] = 9

# ============================================================================
# 配置参数
# ============================================================================

# 固定参数
SLOPE = 30          # 坡度 (度)
SZA = 30            # 太阳天顶角 (度)
SAA = 0             # 太阳方位角 (度)
BAND = 'NIR'        # 波段

# NIR参数
RHO_L_NIR = 0.55    # 叶片反射率
TAU_L_NIR = 0.405    # 叶片透射率
SOIL_R_NIR = 0.211   # 土壤反射率

# 变化参数
ASPECTS = [0, 90, 180, 270]      # 坡向 (度)

# 行1：LAI变化
LAI_VALUES = [0.5, 1, 2, 4, 6]

# 行2：叶倾角 - 6种叶倾角分布类型
LEAF_CLASSES = {
    1: {'name': 'Planophile', 'params': (1, 2)},
    2: {'name': 'Erectophile', 'params': (-1, -2)},
    3: {'name': 'Plagiophile', 'params': (-1, 4)},
    4: {'name': 'Extremophile', 'params': (1, 4)},
    5: {'name': 'Uniform', 'params': (0, 0)},
    6: {'name': 'Spherical', 'params': (0, 0)},  # Special case
}

# 行3：天空散射比 (sky_ratio = 直接光占比, 1-sky_ratio = 散射光占比)
SKY_RATIOS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]  # 从完全直射到完全散射

# ============================================================================
# 观测方向采样
# ============================================================================

def get_observation_angles(saa, n_samples=5):
    """
    生成观测方向
    从后向散射(-60°)到前向散射(+60°)
    """
    vzas = np.linspace(-60, 60, n_samples)
    vaas = []
    for vza in vzas:
        if vza < 0:
            vaa = saa          # 后向散射
        else:
            vaa = saa + 180    # 前向散射
        vaas.append(vaa)
    
    return vzas, np.array(vaas)


# ============================================================================
# 计算BRF
# ============================================================================

def compute_brf_curve(terrain, leaf_class, sky_ratio, sza, saa, 
                      vzas, vaas, rho_l, tau_l, soil_r):
    """
    计算一条BRF曲线
    """
    brf_values = []
    
    for vza, vaa in zip(vzas, vaas):
            # 获取几何参数
        geo_comps = terrain.get_fast_geometry(terrain.sun_mask, sza, saa, vza, vaa)
            
            # 计算BRF
        brf = MyModel.PATH_RT_Terrain(
                terrain, tau_l, rho_l, soil_r,
                geo_comps, sky_ratio,
                sza, saa, vza, vaa,
                leaf_class=leaf_class, ALA=57.5, Hotspot=0.02
            )
            
        brf_values.append(brf)
    return np.array(brf_values)


# ============================================================================
# 绘图函数
# ============================================================================

def plot_lai_sensitivity(slope, sza, saa, aspect, vzas, vaas,
                         rho_l, tau_l, soil_r, ax, col_idx=0):
    """
    行1：LAI敏感性分析
    """
    vzasabs = np.abs(vzas)
    colors_lai = plt.cm.viridis(np.linspace(0, 1, len(LAI_VALUES)))
    
    for lai, color in zip(LAI_VALUES, colors_lai):
        # 修改terrain的LAI
        terrain = MyModel.TerrainGeometry(
            slope=slope, aspect=aspect,
            homogeneous_h=5,      # 树高
            homogeneous_base=1,   # 枝下高
            mode='homogeneous',
            FAVD=lai/4,
            SZA=SZA, SAA=SAA,
            margin=0
        )
        
        brf = compute_brf_curve(terrain, leaf_class=6, sky_ratio=0,
                              sza=sza, saa=saa, vzas=vzasabs, vaas=vaas,
                              rho_l=rho_l, tau_l=tau_l, soil_r=soil_r)
        
        ax.plot(vzas, brf, marker='o', markersize=1, linewidth=0.5,
               label=f'LAI={lai}', color=color)
    
    ax.axvline(x=-sza, color='red', linestyle='--', alpha=0.5, linewidth=0.5)
    ax.set_xlabel('VZA [°]', fontsize=9)
    if col_idx == 0:
        ax.set_ylabel('Reflectance', fontsize=9)
    else:
        ax.set_yticklabels([])
    ax.set_title('LAI Sensitivity', fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=7)
    ax.set_xlim([-60, 60])
    ax.set_ylim([0.0, 0.8])


def plot_leaf_class_sensitivity(slope, sza, saa, aspect, vzas, vaas,
                               rho_l, tau_l, soil_r, ax, col_idx=0):
    """
    行2：叶倾角敏感性分析
    """
    vzasabs = np.abs(vzas)
    colors_leaf = plt.cm.tab10(np.linspace(0, 1, len(LEAF_CLASSES)))
    
    terrain = MyModel.TerrainGeometry(
            slope=slope, aspect=aspect,
            homogeneous_h=5,      # 树高
            homogeneous_base=1,   # 枝下高
            mode='homogeneous',
            FAVD=0.25,
            SZA=SZA, SAA=SAA,
            margin=0
        )
    
    for (leaf_class, info), color in zip(LEAF_CLASSES.items(), colors_leaf):
        brf = compute_brf_curve(terrain, leaf_class=leaf_class,
                              sky_ratio=0.0, sza=sza, saa=saa,
                              vzas=vzasabs, vaas=vaas,
                              rho_l=rho_l, tau_l=tau_l, soil_r=soil_r)
        
        ax.plot(vzas, brf, marker='s', markersize=1, linewidth=0.5,
               label=info['name'], color=color)

    ax.axvline(x=-sza, color='red', linestyle='--', alpha=0.5, linewidth=0.5)
    ax.set_xlabel('VZA [°]', fontsize=9)
    if col_idx == 0:
        ax.set_ylabel('Reflectance', fontsize=9)
    else:
        ax.set_yticklabels([])
    ax.set_title('Leaf Angle Class Sensitivity', fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=6, ncol=2)
    ax.set_xlim([-60, 60])
    ax.set_ylim([0.0, 0.8])


def plot_sky_ratio_sensitivity(slope, sza, saa, aspect, vzas, vaas,
                              rho_l, tau_l, soil_r, ax, col_idx=0):
    """
    行3：天空散射比敏感性分析
    """
    # sky_ratio: 1.0 = 完全直射, 0.0 = 完全散射
    colors_sky = plt.cm.coolwarm(np.linspace(0, 1, len(SKY_RATIOS)))
    vzasabs = np.abs(vzas)
    terrain = MyModel.TerrainGeometry(
            slope=slope, aspect=aspect,
            homogeneous_h=5,      # 树高
            homogeneous_base=1,   # 枝下高
            mode='homogeneous',
            FAVD=0.25,
            SZA=SZA, SAA=SAA,
            margin=0
        )
    
    for sky_ratio, color in zip(SKY_RATIOS, colors_sky):
        brf = compute_brf_curve(terrain, leaf_class=6,
                              sky_ratio=sky_ratio, sza=sza, saa=saa,
                              vzas=vzasabs, vaas=vaas,
                              rho_l=rho_l, tau_l=tau_l, soil_r=soil_r)
        
        ax.plot(vzas, brf, marker='^', markersize=1, linewidth=0.5,
               label=f'D$_{{\mathrm{{dif}}}}$={sky_ratio:.1f}', color=color)
    
    ax.axvline(x=-sza, color='red', linestyle='--', alpha=0.5, linewidth=0.5)
    ax.set_xlabel('VZA [°]', fontsize=9)
    if col_idx == 0:
        ax.set_ylabel('Reflectance', fontsize=9)
    else:
        ax.set_yticklabels([])
    ax.set_title('Sky Diffuse Sensitivity', fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=7)
    ax.set_xlim([-60, 60])
    ax.set_ylim([0, 0.8])


# ============================================================================
# 主程序
# ============================================================================

def main():
    print("=" * 70)
    print("BRDF 矩阵分析 - NIR波段")
    print("=" * 70)
    print(f"固定参数: 坡度={SLOPE}°, SZA={SZA}°, SAA={SAA}°")
    print(f"分析参数: LAI={LAI_VALUES}, 叶倾角=6类, 天空散射={SKY_RATIOS}")
    print(f"坡向分析: {ASPECTS}°")
    print("=" * 70)
    
    # 获取观测方向
    vzas, vaas = get_observation_angles(SAA, n_samples=121)
    
    # 为每个坡向创建图表
    # 调整尺寸以适应论文排版，4列子图需要更宽的画布
    fig = plt.figure(figsize=(10, 10))
    # fig.suptitle(f'BRDF Analysis: Slope={SLOPE}°, SZA={SZA}°, SAA={SAA}°, Band=NIR\n' + 
    #              'Rows: LAI Sensitivity | Leaf Angle Class | Sky Diffuse',
    #              fontsize=10, fontweight='bold', y=0.995)
    
    for col_idx, aspect in enumerate(ASPECTS):
        print(f"\n处理坡向 {aspect}°...")
        
        # 创建地
        # 行1：LAI敏感性
        ax1 = plt.subplot(3, 4, col_idx + 1)
        print(f"  计算LAI敏感性...")
        plot_lai_sensitivity(SLOPE, SZA, SAA, aspect, vzas, vaas,
                           RHO_L_NIR, TAU_L_NIR, SOIL_R_NIR, ax1, col_idx)
        ax1.set_title(f'LAI\n(Slope=30°, Aspect={aspect}°)', fontsize=10)
        
        # 行2：叶倾角敏感性
        ax2 = plt.subplot(3, 4, col_idx + 5)
        print(f"  计算叶倾角敏感性...")
        plot_leaf_class_sensitivity(SLOPE, SZA, SAA, aspect, vzas, vaas,
                                   RHO_L_NIR, TAU_L_NIR, SOIL_R_NIR, ax2, col_idx)
        ax2.set_title(f'LAD\n(Slope=30°, Aspect={aspect}°)', fontsize=10)
        
        # 行3：天空散射敏感性
        ax3 = plt.subplot(3, 4, col_idx + 9)
        print(f"  计算天空散射敏感性...")
        plot_sky_ratio_sensitivity(SLOPE, SZA, SAA, aspect, vzas, vaas,
                                  RHO_L_NIR, TAU_L_NIR, SOIL_R_NIR, ax3, col_idx)
        ax3.set_title(f'D$_{{\mathrm{{dif}}}}$\n(Slope=30°, Aspect={aspect}°)', fontsize=10)
    
    plt.tight_layout(rect=[0, 0, 1, 0.98], h_pad=0.5, w_pad=0.5)
    
    # 保存图表
    output_path = 'BRDF_Matrix_Analysis.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ 图表已保存: {output_path}")
    
    plt.show()
    print("\n完成!")


if __name__ == "__main__":
    main()
