# ============================================================
# PMSM Torque Optimization via FEMM (Toyota Prius 2004 Model)
# 
# This script automates FEMM to evaluate the electromagnetic 
# torque of a permanent magnet synchronous machine (PMSM) 
# as a function of the current vector angle (Beta). 
# For each current magnitude, the Beta angle that yields 
# the maximum torque is identified.
#
# Author: Oseias de Paula Ferreira
# Institution: Universidade Federal de Minas Gerais (UFMG)
# ============================================================


import femm
import numpy as np
import time

# ================================================
# PMSM Simulation using FEMM - Torque vs Beta Angle
# ================================================

# Open FEMM and load the project
femm.openfemm(0)
femm.opendocument('Toyota-Prius-2004.FEM')
femm.mi_smartmesh(0)  # Coarse mesh for faster simulation
femm.mi_probdef(0, 'millimeters', 'planar', 1.e-8, 83.6, 15, 1)

# ------------------------
# Simulation Parameters
# ------------------------
N_pole_pairs = 4              # Number of pole pairs
N_rotor_poles = 8             # Number of rotor poles
N_theta_steps = 1             # Rotor position steps
N_beta_steps = 15             # Current vector angle steps

# Rotor positions (theta_r) in degrees
theta_r_values = np.linspace(0, 0, N_theta_steps)  # Fixed rotor position

# Stator current magnitudes (A)
iss = [10, 12.96, 21.85, 36.67, 57.41, 84.07, 116.67, 155.19, 199.63, 250]

# Current vector angles Beta (degrees)
BETA = np.linspace(0, 90, N_beta_steps)

# Dictionary to store Beta angle of max torque per current
torque_max_angles = {}

# Vector to store all simulated torque values
torque_vec = []

# ------------------------
# Main Simulation Loop
# ------------------------
for current in iss:
    max_torque = -np.inf
    max_angle = None

    for theta_r in theta_r_values:
        for beta in BETA:
            start_time = time.time()

            # d-q currents
            id_val = -current * np.sin(np.radians(beta))
            iq_val =  current * np.cos(np.radians(beta))

            # Update rotor mechanical angle
            femm.mi_modifyboundprop('SlidingBand', 10, theta_r)

            # Park transformation: dq to abc
            theta_e = (N_rotor_poles / 2) * np.radians(theta_r)  # Electrical angle
            Id = np.array([
                np.sin(theta_e),
                np.sin(theta_e - 2 * np.pi / 3),
                np.sin(theta_e - 4 * np.pi / 3)
            ])
            Iq = np.array([
                np.cos(theta_e),
                np.cos(theta_e - 2 * np.pi / 3),
                np.cos(theta_e - 4 * np.pi / 3)
            ])
            I_abc = id_val * Id + iq_val * Iq

            # Apply phase currents in FEMM
            femm.mi_setcurrent('A', I_abc[0])
            femm.mi_setcurrent('B', I_abc[1])
            femm.mi_setcurrent('C', I_abc[2])

            # Run FEMM analysis
            femm.mi_analyze(0)
            femm.mi_loadsolution()

            # Compute electromagnetic torque
            torque = femm.mo_gapintegral('SlidingBand', 0)
            torque_vec.append(torque)

            # Track maximum torque and associated Beta
            if torque > max_torque:
                max_torque = torque
                max_angle = beta

            elapsed_time = time.time() - start_time
            print(
                f'Theta_r = {theta_r:.2f}°  ::  Current = {current:.2f} A  ::  Beta = {beta:.2f}°  '
                f'::  {elapsed_time:.2f} s  ::  Torque = {torque:.2f} Nm  '
                f'::  id = {id_val:.2f}  ::  iq = {iq_val:.2f}  '
                f'::  ia = {I_abc[0]:.2f}, ib = {I_abc[1]:.2f}, ic = {I_abc[2]:.2f}'
            )
        print()

    # Store angle of max torque for current
    torque_max_angles[current] = max_angle

# ------------------------
# Final Results
# ------------------------
print('Maximum torque angles for each current value:')
for current, beta_max in torque_max_angles.items():
    print(f'Current {current:.2f} A -> Beta = {beta_max:.2f}°')
