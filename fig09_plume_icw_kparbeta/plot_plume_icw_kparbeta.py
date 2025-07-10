import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import scienceplots
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.ticker import LogFormatterMathtext
plt.style.use('science')
# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_columns', None)
# pd.set_option('display.precision', 8)
plt.rcParams.update({'font.size': 25})
plt.rcParams.update({'lines.linewidth': 3.0})

col_names_first_57 = ['kperp', 'kpar', 'betap', 'vtp', 'w', 'g', 
                      'bxr', 'bxi', 'byr', 'byi', 'bzr', 'bzi', 
                      'exr', 'exi', 'eyr', 'eyi', 'ezr', 'ezi',
                      'ux1r', 'ux1i', 'uy1r', 'uy1i', 'uz1r', 'uz1i', 
                      'ux2r', 'ux2i', 'uy2r', 'uy2i', 'uz2r', 'uz2i',
                      'n1r', 'n1i', 'n2r', 'n2i', 
                      'ps1', 'ps2', 
                      'p1yy', 'p1yz', 'p1zy', 'p1zz', 'p1n0', 'p1cd',
                      'p2yy', 'p2yz', 'p2zy', 'p2zz', 'p2n0', 'p2cd',
                      'tau1', 'mu1', 'alph1', 'col_52', 'col_53', 'col_54',  
                      'tau2', 'mu2', 'alph2',
                     ]

# Combine column names and fill in the gap with placeholders
col_names = col_names_first_57 + [f'col_{i+58}' for i in range(60-57)]

def find_closest_row(data, varname, value):
    # Calculate the absolute difference between each value in 'kperp' and the given value
    abs_diff = abs(data[varname] - value)
    
    # Find the index of the smallest absolute difference
    closest_index = abs_diff.idxmin()
    print(closest_index)
    
    # Return the row corresponding to this index
    return data.loc[closest_index]

# Let me try to plot omega to see if they are ICWs
def get_quantities_ICW(data_ICW):
    data_ICW['kpardi'] = data_ICW['kpar'] / np.sqrt(data_ICW['betap'])
    data_ICW['analy_w'] = (np.abs(data_ICW['kpardi'])/2.) * (np.sqrt(data_ICW['kpardi']**2. + 4.) - np.abs(data_ICW['kpardi']))
    data_ICW['waveT'] = 2. * np.pi / data_ICW['w']
    data_ICW['wbar'] = data_ICW['w'] / data_ICW['kpardi']
    data_ICW['vres_n0'] = data_ICW['wbar']/np.sqrt(data_ICW['betap'])
    data_ICW['1kpar'] = 1. / data_ICW['kpar']
    data_ICW['ratio'] = data_ICW['vres_n0'] / data_ICW['1kpar']
    data_ICW['g_over_w'] = data_ICW['g'] / data_ICW['w']
    # n = 1 mode
    data_ICW['vres_n1'] = (data_ICW['w'] - 1.) / data_ICW['kpar']
    # n = -1 mode
    data_ICW['vres_n-1'] = (data_ICW['w'] + 1.) / data_ICW['kpar']
    
    data_ICW['Eright'] = np.sqrt( ( (data_ICW['exr'] + data_ICW['eyi'])**2. + (data_ICW['exi'] - data_ICW['eyr'])**2. ) / 2.)
    data_ICW['Eleft'] = np.sqrt( ( (data_ICW['exr'] - data_ICW['eyi'])**2. + (data_ICW['exi'] + data_ICW['eyr'])**2. ) / 2.)
    data_ICW['V_I'] = (data_ICW['Eright']**2. - data_ICW['Eleft']**2.) / (data_ICW['Eright']**2. + data_ICW['Eleft']**2.)
    data_ICW['PolarE'] = (np.abs(data_ICW['Eright']) - np.abs(data_ICW['Eleft'])) / (np.abs(data_ICW['Eright']) + np.abs(data_ICW['Eleft']))

def log_fmt(value):
    """
    If value is close to 10^n (for some integer n), 
    label as 10^{n}.
    Otherwise, use a normal decimal format.
    """
    # Compute the log10 of the value
    power_float = np.log10(value)
    # Round to the nearest integer
    power_rounded = np.round(power_float)

    # Define a small tolerance to decide if a number
    # is "close enough" to an exact power of ten
    tolerance = 1e-8

    if abs(power_float - power_rounded) < tolerance:
        # It's effectively an exact power of ten
        return f"$10^{{{int(power_rounded)}}}$"
    else:
        # Otherwise, just return the decimal with a suitable format
        # Use 'g' or '.2f' or whichever format you prefer
        return f"{value:g}"



# Load the data and assign column names
filename = "./lmfpc_icw_BetaKpar_betap_kpar.mode1"
data_bkpar = pd.read_table(filename, delim_whitespace=True, header=None, names=col_names)

get_quantities_ICW(data_bkpar)

chosen_wave_modes = np.array(
    [
        [0.52480750, 0.1],
        [0.52480750, 0.3],
        [0.52480750, 1],
        [0.52480750, 3],
        [0.52480750, 10],
    ]
)

# Set up plotting grids
betap_unique = np.sort(data_bkpar['betap'].unique())
kpar_unique = np.sort(data_bkpar['kpar'].unique())

# Recreate the mesh grid for betap and kpar
kpar_grid, betap_grid = np.meshgrid(kpar_unique, betap_unique)

# Compute kpardi grid
kpardi_grid = kpar_grid / np.sqrt(betap_grid)


# =======================================================================
# Plot $- \gamma / \omega$ on 
# $(k_\parallel \rho_i, \beta_i)$ and $(k_\parallel d_i, \beta_i)$ Planes

# Reshape - g_over_w
g_vres_grid = -data_bkpar['g_over_w'].values.reshape(len(kpar_unique), len(betap_unique))

# Plot - g_over_w on (kparrhoi, betap)
# Plot the 2D color plot with logarithmic axes
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111)

c = ax.pcolormesh(kpar_grid, betap_grid, g_vres_grid, cmap=plt.get_cmap('Wistia'), norm=LogNorm(vmin=1e-4, vmax=1e2))

# Add contours
contour_levels = [1e-4, 1e-3, 1e-2, 1e-1, 0.3, 0.7, 1e1, 1e2, 1e3]

cs = ax.contour(kpar_grid, betap_grid, g_vres_grid, levels=contour_levels, colors='black', linestyles='--', linewidths=1)
ax.clabel(cs, inline=False, fontsize=25, fmt=log_fmt)  # Optional: label the contour lines

g_cs = ax.contour(kpar_grid, betap_grid, g_vres_grid, levels=[1], colors='#0000a7', linestyles='-', linewidths=2)
ax.clabel(g_cs, inline=False, fontsize=40, fmt=r'$\boldsymbol{-\gamma/\omega = 1}$')

# Plot betap = 0.1, 0.3, 1, 3, 10 lines
for i in [0.1, 0.3, 1, 3, 10]:
    ax.axhline(i, color='k', linestyle='-', linewidth=1, alpha=0.3)


# Plot kparrhoi = 0.1, 0.2, 0.3, ..., 1 lines
for j in np.arange(0.1, 1.1, 0.1):
    ax.axvline(j, color='k', linestyle='-', linewidth=1, alpha=0.3)

for i in np.arange(chosen_wave_modes.shape[0]):
    ax.plot(chosen_wave_modes[i][0], chosen_wave_modes[i][1], marker='o', color='k', markersize=7)

plt.text(0.04, 15, r"$(a)$", fontsize=30, ha='center', va='center', fontweight='bold')
    
# Set logarithmic scales
ax.set_xscale('log')
ax.set_yscale('log')

# ax.set_xlim(0.1, 10)
# Add color bar and labels
cbar = plt.colorbar(c, ax=ax)
ax.set_xlabel(r"$k_\parallel \rho_i$")
ax.set_ylabel(r"$\beta_i$")
# ax.set_xlim(0.1, 10)
ax.set_title(r"$-\gamma/\omega$")
plt.tight_layout()
plt.savefig("./G_BetaKparRhoi_WithContoursAndWaveModes.pdf", dpi=200)
# plt.show()


# =======================================================================
# Plot $v_{res, n = 1}$ on 
# $(k_\parallel \rho_i, \beta_i)$ and $(k_\parallel d_i, \beta_i)$ Planes

# Reshape Vresn1
vres_grid = data_bkpar['vres_n1'].values.reshape(len(kpar_unique), len(betap_unique))
# Plot vresn1 on (kparrhoi, betap)
# Plot the 2D color plot with logarithmic axes
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111)

c = ax.pcolormesh(kpar_grid, betap_grid, vres_grid, cmap=plt.get_cmap('Wistia'), vmin=-3.1, vmax=0.1)

# Add contours
contour_levels = [-4.0, -3.0, -2.0, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.1]

cs = ax.contour(kpar_grid, betap_grid, vres_grid, levels=contour_levels, colors='black', linestyles='--', linewidths=1)
ax.clabel(cs, inline=False, fontsize=25, fmt='%3.2f')  # Optional: label the contour lines

g_cs = ax.contour(kpar_grid, betap_grid, g_vres_grid, levels=[1], colors='#0000a7', linestyles='-', linewidths=2)
ax.clabel(g_cs, inline=False, fontsize=40, fmt=r'$\boldsymbol{-\gamma/\omega = 1}$')

# Plot betap = 0.1, 0.3, 1, 3, 10 lines
for i in [0.1, 0.3, 1, 3, 10]:
    ax.axhline(i, color='k', linestyle='-', linewidth=1, alpha=0.3)


# Plot kparrhoi = 0.1, 0.2, 0.3, ..., 1 lines
for j in np.arange(0.1, 1.1, 0.1):
    ax.axvline(j, color='k', linestyle='-', linewidth=1, alpha=0.3)

for i in np.arange(chosen_wave_modes.shape[0]):
    ax.plot(chosen_wave_modes[i][0], chosen_wave_modes[i][1], marker='o', color='k', markersize=7)

plt.text(0.04, 15, r"$(b)$", fontsize=30, ha='center', va='center', fontweight='bold')
    
# Set logarithmic scales
ax.set_xscale('log')
ax.set_yscale('log')

# ax.set_xlim(0.1, 10)
# Add color bar and labels
cbar = plt.colorbar(c, ax=ax)
ax.set_xlabel(r"$k_\parallel \rho_i$")
ax.set_ylabel(r"$\beta_i$")
ax.set_title(r"$n = 1, v_{res}/v_{ti}$")
plt.tight_layout()
plt.savefig("Vresn1_BetaKparRhoi_WithContoursAndWaveModes.pdf", dpi=200)
# plt.show()

# =======================================================================
# Plot $\omega/\Omega_i$ on 
# $(k_\parallel \rho_i, \beta_i)$ and $(k_\parallel d_i, \beta_i)$ Planes

# I know that the direct output w is omega/Omega_i
# Do I need to reshape w?
wO = data_bkpar['w'].values.reshape(len(kpar_unique), len(betap_unique))

# Plot - g_over_w on (kparrhoi, betap)
# Plot the 2D color plot with logarithmic axes
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111)

c = ax.pcolormesh(kpar_grid, betap_grid, wO, cmap=plt.get_cmap('Wistia'), vmax=1)

# Add contours
contour_levels = [1e-2, 0.1, 0.2, 0.4, 0.6, 0.8, 1]

cs = ax.contour(kpar_grid, betap_grid, wO, levels=contour_levels, colors='black', linestyles='--', linewidths=1)
ax.clabel(cs, inline=False, fontsize=25, fmt=log_fmt)  # Optional: label the contour lines

g_cs = ax.contour(kpar_grid, betap_grid, g_vres_grid, levels=[1], colors='#0000a7', linestyles='-', linewidths=2)
ax.clabel(g_cs, inline=False, fontsize=40, fmt=r'$\boldsymbol{-\gamma/\omega = 1}$')

# Plot betap = 0.1, 0.3, 1, 3, 10 lines
for i in [0.1, 0.3, 1, 3, 10]:
    ax.axhline(i, color='k', linestyle='-', linewidth=1, alpha=0.3)


# Plot kparrhoi = 0.1, 0.2, 0.3, ..., 1 lines
for j in np.arange(0.1, 1.1, 0.1):
    ax.axvline(j, color='k', linestyle='-', linewidth=1, alpha=0.3)

for i in np.arange(chosen_wave_modes.shape[0]):
    ax.plot(chosen_wave_modes[i][0], chosen_wave_modes[i][1], marker='o', color='k', markersize=7)


plt.text(0.04, 15, r"$(c)$", fontsize=30, ha='center', va='center', fontweight='bold')
# Set logarithmic scales
ax.set_xscale('log')
ax.set_yscale('log')

# ax.set_xlim(0.1, 10)
# Add color bar and labels
cbar = plt.colorbar(c, ax=ax)
ax.set_xlabel(r"$k_\parallel \rho_i$")
ax.set_ylabel(r"$\beta_i$")
# ax.set_xlim(0.1, 10)
ax.set_title(r"$\omega/\Omega_i$")
plt.tight_layout()
plt.savefig("wO_BetaKparRhoi_WithContoursAndWaveModes.pdf", dpi=200)
# plt.show()

# =======================================================================
# Plot $\mathcal{P}_E$ on 
# $(k_\parallel \rho_i, \beta_i)$ and $(k_\parallel d_i, \beta_i)$ Planes

# Reshape `Ep` into the same shape as the grids
vres_grid = data_bkpar['PolarE'].values.reshape(len(kpar_unique), len(betap_unique))

# Plot the 2D color plot with logarithmic axes
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111)
c = ax.pcolormesh(kpar_grid, betap_grid, vres_grid, cmap='Wistia', vmax = 1)

g_cs = ax.contour(kpar_grid, betap_grid, g_vres_grid, levels=[1], colors='#0000a7', linestyles='-', linewidths=2)
ax.clabel(g_cs, inline=False, fontsize=40, fmt=r'$\boldsymbol{-\gamma/\omega = 1}$')


# Plot betap = 0.1, 0.3, 1, 3, 10 lines
for i in [0.1, 0.3, 1, 3, 10]:
    ax.axhline(i, color='k', linestyle='-', linewidth=1, alpha=0.3)


# Plot kparrhoi = 0.1, 0.2, 0.3, ..., 1 lines
for j in np.arange(0.1, 1.1, 0.1):
    ax.axvline(j, color='k', linestyle='-', linewidth=1, alpha=0.3)


for i in np.arange(chosen_wave_modes.shape[0]):
    ax.plot(chosen_wave_modes[i][0], chosen_wave_modes[i][1], marker='o', color='k', markersize=7)

plt.text(0.04, 15, r"$(d)$", fontsize=30, ha='center', va='center', fontweight='bold')

# Set logarithmic scales
ax.set_xscale('log')
ax.set_yscale('log')
# ax.set_xlim(0.1, 10)
# Add color bar and labels
cbar = plt.colorbar(c, ax=ax)
ax.set_xlabel(r"$k_\parallel \rho_i$")
ax.set_ylabel(r"$\beta_i$")
ax.set_title(r"$\mathcal{P}_E$", fontsize=25)
plt.tight_layout()
plt.savefig("PolarE_BetaKparRhoi_WithContoursAndWaveModes.pdf", dpi=200)
# plt.show()



