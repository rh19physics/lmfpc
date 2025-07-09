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
filename = "./lmfpc_b1_kpar005_mime1836_kperp_1_100000.mode1"
data = pd.read_table(filename, delim_whitespace=True, header=None, names=col_names)



# Plot wave frequency versus perpendicular wavevector
plt.figure(figsize=(10, 3))
plt.loglog(data['kperp'], np.sqrt(data['betap'])*data['w']/data['kpar'], color='k', label="data from PLUME")
# plt.loglog(data['kperp'], np.sqrt(1 + data['kperp']**2 / data['betap']), color="r", label=r"$\sqrt{1 + \frac{(k_\perp \rho_i)^2}{\beta_i + \frac{2}{1 + T_e / T_i}}}$")
plt.axvline(1, color='k', linestyle="dashed")
plt.text(0.12, 12, r"$(a)$", fontsize=30, ha='center', va='center', fontweight='bold')
plt.xlabel(r"$k_\perp \rho_i$")
plt.ylabel(r"$\omega/(k_\parallel v_A)$")
plt.grid()
plt.xlim(1e-1, 1e1)
plt.savefig("wbar_kpar005.pdf", dpi=200)
# plt.show()


# Plot total and separated damping rates versus perpendicular wavevector
plt.figure(figsize=(10, 7))
plt.axhline(1, color="gray", linestyle="dotted", linewidth=2)
plt.axhline(0.1, color="gray", linestyle="dashed", linewidth=2)
plt.loglog(data['kperp'], -data['g']/data['w'], color='k', linestyle="-", label=r"$\gamma$")
plt.loglog(data['kperp'], data['ps1'], color="#c1272d", label=r"$\gamma_i$", linewidth=2)
plt.loglog(data['kperp'], data['p1zy'] + data['p1zz'], label=r"$\gamma_{i, LD}$", color="#c1272d", linestyle="dashed")
plt.loglog(data['kperp'], data['ps2'], color="#065EDA", label=r"$\gamma_e$", linewidth=2)
# plt.loglog(data['kperp'], data['p1yy'] + data['p1yz'], label=r"$\gamma_{i, TTD}$", color="r", linestyle="dashdot")
# plt.loglog(data['kperp'], data['p1cd'], label=r"$\gamma_{i, CD}$", color="g", linestyle="dashed")
plt.axvline(1, color="gray", linestyle="--", linewidth=2)
plt.text(0.12, 5, r"$(b)$", fontsize=30, ha='center', va='center', fontweight='bold')
plt.legend(fontsize=21)
plt.xlabel(r"$k_\perp \rho_i$")
plt.ylabel(r"$|\gamma|/\omega$")
plt.grid()
plt.xlim(1e-1, 1e1)
plt.ylim(1e-6, 1e1)
plt.savefig("GoverW_KAW.pdf", dpi=200)
# plt.show()



