# ============================================================
# Power Loss Computation in PMSM - FEMM-Based Simulation
#
# This script estimates various power losses in a Surface 
# Permanent Magnet Synchronous Machine (SPM), including:
#   - Hysteresis losses
#   - Eddy current losses (iron)
#   - Proximity effect losses
#   - Ohmic (Joule) losses in windings
#   - Eddy current losses in permanent magnets
#
# Author: Oseias de Paula Ferreira
# Institution: Universidade Federal de Minas Gerais (UFMG)
# Date: May 2025
# ============================================================

import numpy as np
import math

# ------------------------------
# Operating point setup
# ------------------------------
SpeedMin = 400      # [RPM]
SpeedMax = 4000     # [RPM]
SpeedStep = 400
NSpeed = 9
RotorMagnets = 16
omag = 0.556e6      # Magnet conductivity in [S/m]

# ------------------------------
# Core loss coefficients (in steel)
# ------------------------------
p = 7700            # Density [kg/m³]
Const = 0.45        # lb to kg conversion
ch = 0.0084         # Hysteresis loss coefficient [W/(lb*T²*Hz)]
ce = 3.12e-5        # Eddy current loss coefficient [W/(lb*T²*Hz²)]
cs = 0.97           # Empirical scaling factor

# Convert coefficients to SI units
ch *= p / Const     # Now in [W/m³/T²/Hz]
ce *= p / Const

# ------------------------------
# Proximity effect parameters (AWG-based)
# ------------------------------
AWG = 25
Winding_Fill = 0.3882
PhaseResistance = 0.7742  # @ 20°C
TemperatureRise = 90

dwire = 0.324861 * 0.0254 * math.exp(-0.115942 * AWG)
owire = 58e6 / (1 + TemperatureRise * 0.004)
cePhase = (math.pi**2 / 8) * dwire**2 * Winding_Fill * owire

# ------------------------------
# Ohmic losses in phase windings
# ------------------------------
Iphase = np.array(Is_vec / np.sqrt(2))
PhaseOhmicLoss = np.zeros((N, round(n / dk)))
for i in range(N):
    for kk in range(int(n / dk)):
        PhaseOhmicLoss[i, kk] = (Iphase[i, kk])**2 * (3 * PhaseResistance)

# ------------------------------
# Remove average current density from magnets
# ------------------------------
Jm_vec = Jm.copy()
for i in range(N):
    for k in range(1, RotorMagnets + 1):
        g3 = (g == (10 + k))
        vmag = np.dot(v.T, g3)
        vg3 = v.T * g3
        vg3 = vg3[:, np.newaxis]
        Jo = (Jm_vec[i] @ vg3) / vmag
        Jm_vec[i] = Jm_vec[i] - Jo * g3

# ------------------------------
# Initialize result matrices
# ------------------------------
total_loss_matrix = np.zeros((N, NSpeed))
stator_loss_matrix = np.zeros((N, NSpeed))
rotor_loss_matrix = np.zeros((N, NSpeed))
magnet_loss_matrix = np.zeros((N, NSpeed))
phase_loss_matrix = np.zeros((N, NSpeed))
prox_loss_matrix = np.zeros((N, NSpeed))
cobre_loss_matrix = np.zeros((N, NSpeed))  # Redundant here, but placeholder

Speed = np.zeros(NSpeed)
Pin = np.zeros((N, NSpeed))
Tin = np.zeros((N, NSpeed))
Pout = np.zeros((N, NSpeed))
Tout = np.zeros((N, NSpeed))
Tloss = np.zeros((N, NSpeed))
Eff = np.zeros((N, NSpeed))

# ------------------------------
# Loop for power loss estimation
# ------------------------------
for i in range(N):
    Bsq = bsq[i]
    Bsq_ex = bsq_ex[i]
    Bsq_btt = bsq_btt[i]
    jm = Jm_vec[i]
    OhmicLoss = PhaseOhmicLoss[i]
    Torque_value = np.mean(Tq[i])

    for j, thisSpeed in enumerate(np.linspace(SpeedMin, SpeedMax, NSpeed)):
        thisFreq = thisSpeed / 60
        w1 = np.arange(ns)
        w = MyLowestHarmonic * thisFreq * w1 * (w1 < (ns / 2))

        # Region masks
        g1 = (g == 10).astype(float)  # Rotor
        g2 = (g == 1).astype(float)   # Stator
        g4 = (g == 2).astype(float)   # Windings

        # Loss
