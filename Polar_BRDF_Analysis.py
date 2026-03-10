import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import matplotlib.patches as patches
import PATH_RT_T as MyModel

# Set global font to Arial
plt.rcParams['font.family'] = 'arial'
plt.rcParams['font.size'] = 10

# ==============================================================================
# Polar BRDF Analysis - Slope and Aspect Sensitivity Analysis
# ==============================================================================

def calculate_polar_brdf(slope, aspect, sza=30, saa=0, FAVD=0.25, 
                        rho_l=0.55, tau_l=0.405, soil_r=0.211):
    """
    Calculate BRF values for different viewing directions
    
    Parameters:
        slope:   Slope angle (degrees)
        aspect:  Aspect angle (degrees, 0=North)
        sza:     Solar zenith angle (degrees)
        saa:     Solar azimuth angle (degrees, 0=North)
        FAVD:    Foliage Area Volume Density
        rho_l:   Leaf reflectance
        tau_l:   Leaf transmittance
        soil_r:  Soil reflectance
    
    Returns:
        vza_array:   Viewing zenith angle array (degrees)
        vaa_array:   Viewing azimuth angle array (degrees)
        brf_matrix:  BRF value matrix (shape: [n_vza, n_vaa])
    """
    
    # Create terrain
    terrain = MyModel.TerrainGeometry(
        mode='homogeneous',
        homogeneous_h=5,      # Tree height
        homogeneous_base=1,   # Branch height
        FAVD=FAVD,            # Foliage Area Volume Density
        scale=100, res=0.5,
        slope=slope, aspect=aspect, margin=0
    )
    
    # Viewing direction sampling (polar coordinates)
    n_vza = 60  # Zenith angle sampling (0-60 degrees)
    n_vaa = 240  # Azimuth angle sampling (0-360 degrees)
    
    vza_array = np.linspace(0, 90, n_vza)  # Zenith angle: 0-60 degrees
    vaa_array = np.linspace(0, 360, n_vaa, endpoint=True)  # Azimuth angle: 0-360 degrees (including endpoints for循环)
    
    # Create BRF matrix
    brf_matrix = np.zeros((n_vza, n_vaa))
    
    print(f"  Calculating BRDF for slope {slope}°, aspect {aspect}°...")
    
    # Traverse all viewing directions
    for i, vza in enumerate(vza_array):
        for j, vaa in enumerate(vaa_array):
            try:
                # Get geometric parameters
                geo_comps = terrain.get_fast_geometry(
                    terrain.sun_mask, sza, saa, vza, vaa
                )
                
                # Calculate BRF (NIR band)
                brf = MyModel.PATH_RT_Terrain(
                    terrain, tau_l, rho_l, soil_r,
                    geo_comps, 0, sza, saa, vza, vaa, 
                    leaf_class=6, Hotspot=0.02
                )
                brf_matrix[i, j] = brf
            except Exception as e:
                print(f"    Warning: Slope {slope}, aspect {aspect}, VZA {vza}, VAA {vaa} calculation failed: {e}")
                brf_matrix[i, j] = np.nan
    
    return vza_array, vaa_array, brf_matrix


def plot_polar_brdf_grid(slopes, aspects, sza = 30, saa = 0, FAVD = 0.25,
                         rho_l = 0.55, tau_l = 0.405, soil_r = 0.211):
    """
    Plot polar BRDF grid (slope x aspect)
    
    Parameters:
        slopes:   Slope angle array [degrees]
        aspects:  Aspect angle array [degrees], 0=North
        sza, saa, FAVD, rho_l, tau_l, soil_r: Fixed parameters
    """
    from matplotlib.gridspec import GridSpec
    
    n_slopes = len(slopes)
    n_aspects = len(aspects)
    
    print("\n" + "="*60)
    print(f"开始计算BRDF (SZA={sza}°, SAA={saa}°, LAI={0.25})")
    print(f"坡度: {slopes}")
    print(f"坡向: {aspects}")
    print("="*60)
    
    # 先计算所有数据以确定值范围
    all_data = {}
    for slope in slopes:
        for aspect in aspects:
            vza, vaa, brf = calculate_polar_brdf(
                slope, aspect, sza, saa, FAVD, rho_l, tau_l, soil_r
            )
            all_data[(slope, aspect)] = (vza, vaa, brf)
    
    # 确定全局值范围
    all_brf_values = np.concatenate([
        data[2].flatten() for data in all_data.values()
    ])
    vmin = np.nanmin(all_brf_values)
    vmax = np.nanmax(all_brf_values)
    print(f"\nBRF值范围: [{vmin:.4f}, {vmax:.4f}]")
    
    # Create figure and GridSpec
    fig = plt.figure(figsize=(10, 7))
    
    # GridSpec: rows=slope count+2 (top title+row label), columns=aspect count+1 (left label+columns)
    gs = GridSpec(n_slopes, n_aspects, figure=fig,
                  left=0.12, right=0.88, top=0.92, bottom=0.08,
                  wspace=0.25, hspace=0.3)
    
    # Define colormap
    cmap = cm.get_cmap('RdYlGn_r')
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    
    # 绘制子图
    for i, slope in enumerate(slopes):
        for j, aspect in enumerate(aspects):
            ax = fig.add_subplot(gs[i, j], projection='polar')
            
            vza, vaa, brf = all_data[(slope, aspect)]
            
            # 转换为极坐标网格
            vaa_rad = np.deg2rad(vaa)
            r = vza
            
            # 创建网格
            R, Theta = np.meshgrid(r, vaa_rad, indexing='ij')
            
            # 绘制contourf
            levels = np.linspace(vmin, vmax, 60)
            contour = ax.contourf(Theta, R, brf, levels=levels, cmap=cmap, norm=norm, linewidths=0.5)
            
            # 标注太阳位置
            sun_r = sza
            sun_theta = np.deg2rad(saa)
            ax.plot(sun_theta, sun_r, marker='*', markersize=6, 
                   color='yellow', markeredgecolor='black', markeredgewidth=0.5, markerfacecolor = None,
                   zorder=100)
            ax.set_theta_offset(np.pi/2)
            ax.xaxis.set_tick_params(pad=-4)
            # 设置径向标签
            ax.set_ylim(0, 90)
            ax.set_yticks([0, 30, 60, 90])
            ax.set_yticklabels(['0°', '30°', '60°', '90°'], fontsize=8, family='arial')
            ax.tick_params(axis='y', labelsize=8)
            
            # 设置角度标签
            ax.set_theta_zero_location('N')
            ax.set_theta_direction(-1)
            ax.set_xticks(np.deg2rad([0, 90, 180, 270]))
            ax.set_xticklabels(['N', 'E', 'S', 'W'], fontsize=8, family='arial')
            ax.tick_params(axis='x', labelsize=8)
            
            # 调整极坐标外轮廓线宽度
            ax.spines['polar'].set_linewidth(0.5)
            
            # ax.set_xticklabels(ax.get_xticklabels(), frac=1.05)
            # 在子图内部显示参数 (右上角)
            # param_text = f'SZA={sza}°\nSAA={saa}°'
            # ax.text(0.98, 0.98, param_text, transform=ax.transAxes,
            #        fontsize=6, family='Times New Roman',
            #        ha='right', va='top', 
            #        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', 
            #                 alpha=0.7, edgecolor='none'))
            
            # 添加网格
            ax.grid(True, linestyle='--', alpha=0.3, linewidth=0.5)
            
            # 添加列标题 (坡向)
            if i == 0:
                aspect_name = {0: 'N', 90: 'E', 180: 'S', 270: 'W'}[aspect]
                ax.text(0.5, 1.1, f'Aspect {aspect}°', transform=ax.transAxes,
                       fontsize=10, family='arial', fontweight='bold',
                       ha='center', va='bottom')
            
            # 添加行标题 (坡度) - 在左侧
            if j == 0:
                ax.text(-0.1, 0.5, f'Slope {slope}°', transform=ax.transAxes,
                       fontsize=10, family='arial', fontweight='bold',
                       ha='right', va='center', rotation=90)
    
    # 添加colorbar - 在右侧
    cbar_ax = fig.add_axes([0.90, 0.12, 0.012, 0.78])
    cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), 
                        cax=cbar_ax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Reflectance', fontsize=10, family='arial')
    
    # 总标题
    # fig.text(0.5, 0.98, 
    #         f'BRDF Polar Plot Analysis (NIR, LAI = 1)',
    #         ha='center', fontsize=10, family='Times New Roman', fontweight='bold')
    
    return fig, all_data


# ==============================================================================
# Main Program
# ==============================================================================

if __name__ == "__main__":
    # Parameter settings
    SLOPES = [15, 30, 45]              # Slope (degrees)
    ASPECTS = [0, 90, 180, 270]        # Aspect (degrees, 0=North, 90=East, 180=South, 270=West)
    SZA = 30                           # Solar zenith angle (degrees)
    SAA = 0                            # Solar azimuth angle (degrees)
    FAVD = 0.25                        # Foliage Area Volume Density
    
    # NIR band optical parameters
    RHO_L_NIR = 0.55
    TAU_L_NIR = 0.405
    SOIL_R_NIR = 0.211
    
    # Plot polar diagram
    fig, data = plot_polar_brdf_grid(
        slopes=SLOPES,
        aspects=ASPECTS,
        sza=SZA,
        saa=SAA,
        FAVD=FAVD,
        rho_l=RHO_L_NIR,
        tau_l=TAU_L_NIR,
        soil_r=SOIL_R_NIR
    )
    
    plt.show()
    
    # Save image
    output_path = r'D:\PhD\matlab\TerrainBRDF\output\polar_brdf_analysis.png'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nImage saved: {output_path}")
