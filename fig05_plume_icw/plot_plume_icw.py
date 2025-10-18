import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import scienceplots
plt.style.use('science')
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.precision', 8)
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

# Load the data and assign column names
filename = "./lmfpc_b1_kperp001_mime1836_kpar_10_10000.mode1"
data = pd.read_table(filename, delim_whitespace=True, header=None, names=col_names)

# Calculate other quantities
data['kpardi'] = data['kpar'] / np.sqrt(data['betap'])
data['analy_w'] = (np.abs(data['kpardi'])/2.) * (np.sqrt(data['kpardi']**2. + 4.) - np.abs(data['kpardi']))
data['wbar'] = data['w'] / data['kpardi']
data['g_over_w'] = data['g'] / data['w']
# n = 1 mode
data['vres_n1'] = (data['w'] - 1.) / data['kpar']
# n = -1 mode
data['vres_n-1'] = (data['w'] + 1.) / data['kpar']

data['Eright'] = np.sqrt( ( (data['exr'] + data['eyi'])**2. + (data['exi'] - data['eyr'])**2. ) / 2.)
data['Eleft'] = np.sqrt( ( (data['exr'] - data['eyi'])**2. + (data['exi'] + data['eyr'])**2. ) / 2.)
data['Ep_VoverI'] = (data['Eright']**2. - data['Eleft']**2.) / (data['Eright']**2. + data['Eleft']**2.)
data['Ep'] = (data['Eright'] - data['Eleft']) / (data['Eright'] + data['Eleft'])


# Plot wave frequency versus perpendicular wavevector
plt.figure(figsize=(10, 3))
plt.loglog(data['kpar'], data['w'], color='k', label="PLUME")
plt.loglog(data['kpar'], data['analy_w'], color="r", linestyle="dashed", label="Empirical")
# plt.loglog(data['kpar'], data['analy_w'], color="r", linestyle="dashed", label=r"$\frac{k_\parallel \rho_i}{2 \sqrt{\beta_i}} \left[\sqrt{\left(\frac{k_\parallel \rho_i}{\sqrt{\beta_i}}\right)^2 + 4} - \frac{k_\parallel \rho_i}{\sqrt{\beta_i}} \right]$")
plt.axvline(0.525, color='k', linestyle="dashed")
plt.text(0.01, 0.7, r"$(a)$", fontsize=30, ha='center', va='center', fontweight='bold')
plt.legend(fontsize=24)
plt.xlabel(r"$k_\parallel \rho_i$")
plt.ylabel(r"$\omega/\Omega_i$")
plt.grid()
plt.savefig("wbar_wAnalyticalGKwbar_ICW_narrow.pdf", dpi=200)
# plt.show()

# Plot eletric field polarization
plt.figure(figsize=(10, 3))
plt.plot(data['kpar'], data['Ep'], color='k')
plt.axvline(0.525, color='k', linestyle="dashed")
plt.text(0.01, -0.15, r"$(b)$", fontsize=30, ha='center', va='center', fontweight='bold')
# plt.legend()
plt.xscale('log')
plt.xlabel(r"$k_\parallel \rho_i$")
plt.ylabel(r"$\mathcal{P}_E$")
# plt.xlim(1e-1, 1e1)
plt.grid()
plt.savefig("PE_ICW.pdf", dpi=200)
# plt.show()


# Plot total and separated damping rates versus perpendicular wavevector
plt.figure(figsize=(10, 7))
plt.axhline(1, color="gray", linestyle="dotted", linewidth=2)
plt.axhline(0.1, color="gray", linestyle="dashed", linewidth=2)
plt.loglog(data['kpar'], -data['g_over_w'], color='k', linestyle="dashed", label=r"$\gamma$")
plt.loglog(data['kpar'], data['ps1'], color="r", label=r"$\gamma_i$", linewidth=2)
plt.loglog(data['kpar'], data['p1zy'] + data['p1zz'], label=r"$\gamma_{i, LD}$", color="r", linestyle="dashed")
plt.loglog(data['kpar'], data['p1yy'] + data['p1yz'], label=r"$\gamma_{i, TTD}$", color="r", linestyle="dotted")
# plt.loglog(data['kpar'], data['p1yy'] + data['p1yz'], label=r"$\gamma_{i, TTD}$", color="c", linestyle="dash")
plt.loglog(data['kpar'], data['ps2'], color="b", label=r"$\gamma_e$", linewidth=2)
plt.loglog(data['kpar'], data['p1cd'], label=r"$\gamma_{i, CD}$", color="g", linestyle="dashed")
plt.axvline(0.525, color="gray", linestyle="--", linewidth=2)
plt.text(0.01, 0.7, r"$(c)$", fontsize=30, ha='center', va='center', fontweight='bold')

plt.legend(fontsize=24)
plt.xlabel(r"$k_\parallel \rho_i$")
plt.ylabel(r"$|\gamma|/\omega$")
plt.grid()
plt.ylim(1e-6, 2)
plt.savefig("GoverW_ICW.pdf", dpi=200)
# plt.show()



