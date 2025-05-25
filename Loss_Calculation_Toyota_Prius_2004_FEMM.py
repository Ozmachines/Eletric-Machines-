# ============================================================
# Power Loss Estimation in PMSM - FEMM Post-Processing
#
# This script calculates losses in a PMSM using FEMM field 
# simulation outputs (B, A, Jm). The losses include:
#   - Stator and rotor iron losses (hysteresis + eddy current)
#   - Proximity effect losses in windings
#   - Ohmic (Joule) losses in windings
#   - Eddy current losses in permanent magnets
#
# Author: Oseias de Paula Ferreira
# Institution: Universidade Federal de Minas Gerais (UFMG)
# Date: May 2025
# ============================================================

# -----------------------------
# Operating Conditions
# -----------------------------
SpeedMin = 500        # RPM
SpeedMax = 1200       # RPM
SpeedStep = 120
NSpeed = 9
omag = 0.556e6        # Magnet conductivity (S/m)

# -----------------------------
# Core Loss Coefficients (Steel)
# -----------------------------
p = 7700              # [kg/m³]
Const = 0.45          # lb to kg conversion
ch = 0.001105         # [W/(lb*T²*Hz)]
ce = 2.872e-6         # [W/(lb*T²*Hz²)]
cex = 0.002725        # [W/(lb*T^1.5*Hz^1.5)]
cs = 0.95             # Correction factor

# Convert to SI units
ch *= p / Const
ce *= p / Const
cex *= p / Const

# -----------------------------
# Wire Geometry and Loss Setup
# -----------------------------
AWG = 19
Winding_Fill = 0.53
PhaseResistance = 0.05146
TemperatureRise = 90
MyLowestHarmonic = 1

dwire = 0.324861 * 0.0254 * math.exp(-0.115942 * AWG)
owire = 58e6 / (1 + TemperatureRise * 0.004)
cePhase = (math.pi**2 / 8) * dwire**2 * Winding_Fill * owire

# -----------------------------
# DC Ohmic Losses
# -----------------------------
Iphase = np.array(Is_vec / np.sqrt(2))
PhaseOhmicLoss = np.zeros((N, round(n / dk)))
for i in range(N):
    for kk in range(int(n / dk)):
        PhaseOhmicLoss[i, kk] = (Iphase[i, kk])**2 * (3 * PhaseResistance)

# -----------------------------
# Magnet Current Correction
# -----------------------------
Jm_vec = Jm.copy()
for i in range(N):
    g3 = (g == 11).astype(float)
    vmag = v @ g3
    Jo = np.dot(Jm_vec[i], v * g3) / vmag
    Jm_vec[i] -= np.outer(Jo, g3)

# -----------------------------
# Allocate Result Matrices
# -----------------------------
total_loss_matrix = np.zeros((N, NSpeed))
stator_loss_matrix = np.zeros((N, NSpeed))
rotor_loss_matrix = np.zeros((N, NSpeed))
magnet_loss_matrix = np.zeros((N, NSpeed))
phase_loss_matrix = np.zeros((N, NSpeed))
prox_loss_matrix = np.zeros((N, NSpeed))
Speed = np.zeros(NSpeed)

Pin = np.zeros((N, NSpeed))
Tin = np.zeros((N, NSpeed))
Pout = np.zeros((N, NSpeed))
Tout = np.zeros((N, NSpeed))
Eff = np.zeros((N, NSpeed))

# -----------------------------
# Loss Calculation Loop
# -----------------------------
for i in range(N):
    Bsq = bsq[i]
    Bsq_ex = bsq_ex[i]
    Bsq_btt = bsq_btt[i]
    jm = Jm_vec[i]
    OhmicLoss = PhaseOhmicLoss[i]
    Torque_value = np.mean(Tq[i][:])

    for j, thisSpeed in enumerate(np.linspace(SpeedMin, SpeedMax, NSpeed)):
        thisFreq = thisSpeed / 60
        w1 = np.arange(ns)
        w = MyLowestHarmonic * thisFreq * w1 * (w1 < ns / 2)

        # Region masks
        g1 = (g == 4).astype(float)  # Rotor
        g2 = (g == 5).astype(float)  # Stator
        g4 = (g == 2).astype(float)  # Coils

        # Loss calculations
        rotor_loss = 8 * ((ch * w + ce * w**2) @ Bsq @ (v * g1)) / cs
        stator_loss = 8 * ((ch * w + ce * w**2) @ Bsq @ (v * g2)) / cs
        prox_loss = 8 * (cePhase * w**2) @ Bsq @ (v * g4)
        magnet_loss = 8 * 0.5 * (omag * (2 * np.pi * w)**2) @ (np.abs(jm)**2) @ v
        phase_loss = OhmicLoss

        # Total loss
        total_loss = rotor_loss + stator_loss + prox_loss + phase_loss + magnet_loss

        # Save losses
        total_loss_matrix[i, j] = total_loss
        stator_loss_matrix[i, j] = stator_loss
        rotor_loss_matrix[i, j] = rotor_loss
        magnet_loss_matrix[i, j] = magnet_loss
        phase_loss_matrix[i, j] = phase_loss
        prox_loss_matrix[i, j] = prox_loss
        Speed[j] = thisSpeed

        # Input/output power and torque
        omega = np.pi * thisSpeed / 30
        Pout[i, j] = Torque_value * omega
        Pin[i, j] = Pout[i, j] + total_loss
        Tout[i, j] = Pin[i, j] / omega
        Tin[i, j] = (Pin[i, j] - total_loss) / omega
        Eff[i, j] = 100 * Pout[i, j] / Pin[i, j]

# -----------------------------
# Print Results for i = 9
# -----------------------------
index = 9
print(f"Total Loss at index {index}, last speed: {total_loss_matrix[index][-1]:.2f} W")
print(f"Stator Loss: {stator_loss_matrix[index][-1]:.2f} W")
print(f"Rotor Loss: {rotor_loss_matrix[index][-1]:.2f} W")
print(f"Magnet Loss: {magnet_loss_matrix[index][-1]:.2f} W")
print(f"Phase Loss: {phase_loss_matrix[index][-1]:.2f} W")
print(f"Proximity Loss: {prox_loss_matrix[index][-1]:.2f} W")
