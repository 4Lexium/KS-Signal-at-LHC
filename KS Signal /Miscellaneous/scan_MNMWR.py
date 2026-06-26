import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the function
def f(x, y):
    # return (x**2-y**2)/(2*x) -y
    return 0.5*(x-2*y) -0.5*y**2/x

# Create a grid of points
x = np.linspace(3000, 7000, 100)  # x range: 3000 to 6000
y = np.linspace(100, 2000, 100)   # y range: 100 to 1500
X, Y = np.meshgrid(x, y)
Z = f(X, Y)

# Create a single figure with contour plot
fig, ax = plt.subplots(1, 1, figsize=(10, 8))

# Filled contour plot
contourf = ax.contourf(X, Y, Z, levels=20, cmap='plasma')
cbar = plt.colorbar(contourf, ax=ax, label=r'$p_T(\ell_1) - p_T(\ell_2)$')
cbar.set_label(r'$p_T(\ell_1) - p_T(\ell_2)$', size=14)  

# Add contour lines
contour_lines = ax.contour(X, Y, Z, levels=20, colors='white', linewidths=0.5, alpha=0.7)
ax.clabel(contour_lines, inline=True, fontsize=8, colors='white')

# Labels and title
ax.set_xlabel(r'$m_{W_R}$ / GeV',  fontsize=14)
ax.set_ylabel(r'$m_N$ / GeV',  fontsize=14)
# ax.set_title('$\Delta p^l_T(m_{W_R},m_N) = |m_{W_R}-2m_N| / 2-m_N^2/2m_{W_R}$')

# Grid and axes
ax.grid(True, alpha=0.3)
# ax.axhline(y=500, color='red', linestyle='--', linewidth=1, alpha=0.5, label='y=500')
# ax.axvline(x=4500, color='red', linestyle='--', linewidth=1, alpha=0.5, label='x=4500')
# ax.legend()

# Set aspect ratio to auto since x and y ranges are very different
ax.set_aspect('auto')

plt.tight_layout()
plt.savefig("/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/Delta_lepPT_contour.png", dpi=300, bbox_inches='tight')
plt.show()

def delta_log_pt(x, y):
    pt1 = (x**2 - y**2) / (2 * x)
    pt2 = y

    # avoid log of non-positive values
    mask = (pt1 > 0) & (pt2 > 0)
    out = np.full_like(pt1, np.nan, dtype=float)
    out[mask] = np.log(pt1[mask]) - np.log(pt2[mask])
    return out

Z_log = delta_log_pt(X, Y)

fig2, ax2 = plt.subplots(figsize=(10, 8))

contourf2 = ax2.contourf(X, Y, Z_log, levels=20, cmap='viridis')
cbar = plt.colorbar(contourf2, ax=ax2, label=r'$\log p_T(\ell_1) - \log p_T(\ell_2)$')
cbar.set_label(r'$\log p_T(\ell_1) - \log p_T(\ell_2)$', size=14)  
contour_lines2 = ax2.contour(
    X, Y, Z_log,
    levels=20,
    colors='white',
    linewidths=0.5,
    alpha=0.7
)
ax2.clabel(contour_lines2, inline=True, fontsize=8, colors='white')

ax2.set_xlabel(r'$m_{W_R}$ / GeV',  fontsize=14)
ax2.set_ylabel(r'$m_N$ / GeV',  fontsize=14)
# ax2.set_title(r'$\log p_T(\ell_1) - \log p_T(\ell_2)$')

ax2.grid(True, alpha=0.3)
# ax2.axhline(y=500, color='red', linestyle='--', linewidth=1, alpha=0.5)
# ax2.axvline(x=4500, color='red', linestyle='--', linewidth=1, alpha=0.5)

plt.tight_layout()
plt.savefig("/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/Delta_lepPT_log_contour.png",dpi=300, bbox_inches='tight')
plt.show()