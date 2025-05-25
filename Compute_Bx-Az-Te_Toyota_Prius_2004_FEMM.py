# ============================================================
# Magnetic Field Mapping - B and A from FEMM Simulation
#
# This script runs a magnetostatic simulation in FEMM to 
# extract the magnetic flux density (B) and magnetic vector 
# potential (A) for a PMSM across different current magnitudes 
# and rotor positions. Results are stored in a .mat file.
#
# Author: Oseias de Paula Ferreira
# Institution: Universidade Federal de Minas Gerais (UFMG)
# Date: May 2025
# ============================================================

import femm
from femm import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import string
import scipy.io as sio
import math

# -----------------------------
# Input Data - Currents
# -----------------------------
iss = [10, 12.96, 21.85, 36.67, 57.41, 84.07, 116.67, 155.19, 199.63, 250]  # [A]

# Open FEMM and load the model
femm.openfemm(1)
femm.opendocument('Toyota-Prius-2004.FEM')
femm.mi_smartmesh(0)
femm.mi_saveas('temp.fem')
femm.mi_probdef(0, 'millimeters', 'planar', 1.e-8, 83.6, 15, 1)

# -----------------------------
# Simulation Parameters
# -----------------------------
Npolos = 4
RotorMagnets = 8
n = 360
dk = n / 30
N = len(iss)
Nsteps = round(n / dk)

# Initialize arrays
k_vec = np.zeros(Nsteps)
Is_vec = np.zeros((N, Nsteps))
tq = np.zeros(Nsteps)
Tq = np.zeros((N, Nsteps))
inic = 0

# -----------------------------
# FEMM Simulation Loop
# -----------------------------
for i in range(N):
    Beta = torque_max_angles[iss[i]]  # Ensure this dict is available beforehand

    ids = -iss[i] * np.sin(np.radians(Beta))
    iqs =  iss[i] * np.cos(np.radians(Beta))

    for j in range(Nsteps):
        starttime = time.time()
        theta_r = j * dk
        k_vec[j] = theta_r
        Is_vec[i, j] = np.sqrt(ids**2 + iqs**2)

        # Update rotor angle
        femm.mi_modifyboundprop('SlidingBand', 10, theta_r)

        # Park transform dq → abc
        theta_e = (RotorMagnets / 2) * np.radians(theta_r)
        Id = np.array([np.sin(theta_e - 2 * k * np.pi / 3) for k in range(3)])
        Iq = np.array([np.cos(theta_e - 2 * k * np.pi / 3) for k in range(3)])
        Itot = ids * Id + iqs * Iq

        femm.mi_setcurrent('A', Itot[0])
        femm.mi_setcurrent('B', Itot[1])
        femm.mi_setcurrent('C', Itot[2])

        femm.mi_analyze(1)
        femm.mi_loadsolution()

        # Collect geometry and allocate on first run
        if inic == 0:
            nn = femm.mo_numelements()
            b = np.zeros((N, Nsteps, nn), dtype=complex)
            A = np.zeros((N, Nsteps, nn))
            z = np.zeros(nn, dtype=complex)
            a = np.zeros(nn)
            g = np.zeros(nn)

            for m in range(1, nn + 1):
                elm = femm.mo_getelement(m)
                z[m - 1] = elm[3] + 1j * elm[4]
                a[m - 1] = elm[5]
                g[m - 1] = elm[6]
        inic = 1

        # Store flux data
        for m in range(nn):
            if g[m] > 10:  # Rotor magnet
                A[i, j, m] = femm.mo_geta(z[m].real, z[m].imag)
            elif g[m] > 0:  # Iron (stator or rotor)
                Bxy = femm.mo_getb(z[m].real, z[m].imag)
                b[i, j, m] = Bxy[0] + 1j * Bxy[1]

        tq[j] = femm.mo_gapintegral('SlidingBand', 0)
        Tq[i, j] = tq[j]

        elapsed_time = time.time() - starttime
        print(f'Theta_r={theta_r:.2f}° :: Corrente={iss[i]:.2f} A :: Beta={Beta:.2f}° '
              f':: Tempo={elapsed_time:.2f}s :: Torque={tq[j]:.2f} Nm :: '
              f'id={ids:.2f}, iq={iqs:.2f} :: ia={Itot[0]:.2f}, ib={Itot[1]:.2f}, ic={Itot[2]:.2f}')
    print()

# -----------------------------
# Mesh Volume Estimation
# -----------------------------
probinfo = femm.mo_getprobleminfo()
h = probinfo[2]            # Depth (mm)
unit_scale = probinfo[3]   # Conversion to meters
v = np.zeros(nn)
for m in range(nn):
    v[m] = a[m] * h * unit_scale**2  # [m³]

# -----------------------------
# Save results to .mat file
# -----------------------------
sio.savemat('simulacao_fem.mat', {
    'iss': np.array(iss),
    'k_vec': k_vec,
    'Is_vec': Is_vec,
    'Tq': Tq,
    'b': b,
    'A': A,
    'z': z,
    'a': a,
    'g': g,
    'v': v,
    'probinfo': probinfo
})
