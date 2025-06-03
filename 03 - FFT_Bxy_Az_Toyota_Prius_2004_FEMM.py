# ============================================================
# Harmonic Analysis via FFT for Magnetic Flux Density (B)
# and Magnetic Vector Potential (A) from FEMM Simulation
#
# This script processes the simulated magnetic field data (B)
# and vector potential (A) by applying the Fast Fourier 
# Transform (FFT) across rotor positions to extract harmonic 
# content. It estimates various derived magnetic quantities 
# used in iron loss models and eddy current loss estimation.
#
# Author: Oseias de Paula Ferreira
# Institution: Universidade Federal de Minas Gerais (UFMG)
# Date: May 2025
# ============================================================

import numpy as np

# ---------------------------------------------
# Prepare Input Data
# ---------------------------------------------
b_flux = np.array(b)      # Complex magnetic flux density (Bx + j*By)
A_flux = np.array(A)      # Magnetic vector potential A (real)
ns = int(n // dk)         # Number of rotor angle steps (FFT points)

# Decompose B field into real and imaginary parts
b_real = np.real(b_flux)
b_imag = np.imag(b_flux)

# Magnitude of B field
b_mod = np.sqrt(b_real**2 + b_imag**2)

# Initialize result arrays
bsq = np.zeros((N, ns, nn))
bsq_mod = np.zeros((N, ns, nn))
bsq_ex = np.zeros((N, ns, nn))
bsq_btt = np.zeros((N, ns, nn))
Jm = np.zeros((N, ns, nn))  # FFT of vector potential A

# Angular frequency (used later for loss estimation)
omag = 0.556e6  # [rad/s], example value

# ---------------------------------------------
# FFT and Harmonic Feature Extraction Loop
# ---------------------------------------------
for i in range(N):
    # FFT along rotor angle steps (axis=0)
    bxfft = np.abs(np.fft.fft(b_real[i], ns, axis=0)) * 2 / ns
    byfft = np.abs(np.fft.fft(b_imag[i], ns, axis=0)) * 2 / ns

    # Compute harmonic loss-related quantities
    bsq[i] = bxfft**2.045 + byfft**2.045       # General iron loss model (e.g., Bertotti)
    bsq_mod[i] = np.sqrt(bsq[i])               # Resultant flux magnitude
    bsq_ex[i] = bxfft**1.5 + byfft**1.5        # Excess loss estimation
    bsq_btt[i] = bxfft**4 + byfft**4           # Eddy current loss approximation

    # FFT of vector potential A (for eddy current modeling)
    Jm[i] = np.fft.fft(A_flux[i], ns, axis=0) * 2 / ns
