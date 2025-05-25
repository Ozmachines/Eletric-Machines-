# ============================================================
# Magnetic Flux Density (B) and Vector Potential (A) Mapping
# using FEMM for a SPM Motor (Surface Permanent Magnet)
#
# This script runs a magnetic simulation in FEMM for different
# stator current magnitudes and rotor angles, applying dq-abc 
# transformations. It records:
#   - Electromagnetic torque
#   - Magnetic flux density vector B
#   - Magnetic vector potential A
#
# Results are stored in a .mat file for post-processing (e.g., FFT).
#
# Author: Oseias de Paula Ferreira
# Institution: Universidade Federal de Minas Gerais (UFMG)
# Date: May 2025
# Requirements: FEMM, Python 3.x, scipy.io, numpy, matplotlib
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

# Motor currents (MotorCAD reference)
Loc = list(string.ascii_uppercase[0:10])
iss = [0.1, 0.12, 0.19, 0.31, 0.48, 0.69, 0.94, 1.25, 1.6, 2]

# Open FEMM and prepare geometry
femm.openfemm(1)
femm.opendocument('my_SPM_motor_AGB.FEM')
femm.mi_smartmesh(0)
femm.mi_saveas('temp.fem')
femm.mi_probdef(0, 'inches', 'planar', 1.e-8, 0.3, 10, 1)
femm.mi_saveas('temp.fem')

# Simulation settings
RotorMagnets = 16
n = 361
dk = n / 30
Nsteps = round(n / dk)

# Initialize storage
N = len(iss)
k_vec = np.zeros(Nsteps)
Is_vec = np.zeros((N, Nsteps))
tq = np.zeros(Nsteps)
Tq = np.zeros((N, Nsteps))
inic = 0

# FEMM group data placeholder (to be filled during simulation)
b = A = z = a = g = v = None

# Simulation loop
for i in range(N):
    Beta = torque_max_angles[iss[i]]  # Optimal beta angle for current (from previous step)
    ids = -iss[i] * np.sin(np.radians(90 + Beta))
    iqs =  iss[i] * np.cos(np.radians(90 + Beta))

    for j in range(Nsteps):
        theta_r = j * dk
        k_vec[j] = theta_r
        Is_vec[i, j] = np.sqrt(ids**2 + iqs**2)

        femm.mi_modifyboundprop('AGE', 10, theta_r)

        # dq → abc transformation
        theta_e = (RotorMagnets / 2) * np.radians(theta_r)
        Id = np.array([np.sin(theta_e - k * 2 * np.pi / 3) for k in range(3)])
        Iq = np.array([np.cos(theta_e - k * 2 * np.pi / 3) for k in range(3)])
        Itot = ids * Id + iqs * Iq

        femm.mi_setcurrent('A', Itot[0])
        femm.mi_setcurrent('B', Itot[1])
        femm.mi_setcurrent('C', Itot[2])

        starttime = time.time()
        femm.mi_analyze(1)
        femm.mo_loadsolution()

        tq[j] = femm.mo_gapintegral('AGE', 0)
        Tq[i, j] = tq[j]

        # Collect mesh data at first step
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

        # Store flux density and vector potential
        for m in range(nn):
            if g[m] > 10:  # Rotor magnet
                A[i, j, m] = femm.mo_geta(z[m].real, z[m].imag)
            elif g[m] > 0:  # Stator/rotor iron
                Bxy = femm.mo_getb(z[m].real, z[m].imag)
                b[i, j, m] = Bxy[0] + 1j * Bxy[1]

        elapsed_time = time.time() - starttime
        print(f'Theta_r={theta_r:.2f}°  ::  Current={iss[i]:.2f} A  ::  Beta={Beta:.2f}°  '
              f'::  {elapsed_time:.2f} s  ::  Torque={tq[j]:.2f} Nm  '
              f'::  id={ids:.2f}  ::  iq={iqs:.2f}  '
              f'::  ia={Itot[0]:.2f}, ib={Itot[1]:.2f}, ic={Itot[2]:.2f}')
    print()

# Retrieve geometry info and compute element volumes
probinfo = femm.mo_getprobleminfo()
nn = femm.mo_numelements()
h = probinfo[2]
unit_scale = probinfo[3]
v = np.zeros(nn)
for m in range(nn):
    v[m] = a[m] * h * unit_scale**2

# Save all data to .mat file
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
