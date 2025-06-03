# ============================================================
# Efficiency Map Plotting: Torque vs Speed
#
# This script visualizes the efficiency map of an electric 
# machine by plotting efficiency contours as a function of 
# torque and speed. It applies a Gaussian filter to smooth 
# out the results and highlights contour levels.
#
# Author: Oseias de Paula Ferreira
# Institution: Universidade Federal de Minas Gerais (UFMG)
# Date: May 2025
# Requirements: numpy, matplotlib, scipy
# ============================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

# -------------------------------
# Configure Plot Appearance
# -------------------------------
def configure_plot():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": "Latin Modern Roman",
        "axes.labelsize": 18,
        "axes.titlesize": 18,
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "legend.fontsize": 20,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "font.weight": "bold",
        "axes.labelweight": "bold"
    })

# -------------------------------
# Create Filled Contour Plot
# -------------------------------
def create_contour_plot(X, Y, Z):
    fig, ax = plt.subplots(figsize=[6.5, 3.5])

    # Clamp minimum efficiency to 40% for visibility
    Z = np.where(Z < 40, 40, Z)

    # Apply Gaussian filter for smoothness
    Z_smoothed = gaussian_filter(Z, sigma=0.8)

    # Filled contour plot
    n_levels = 17
    contour_filled = ax.contourf(X, Y, Z_smoothed, n_levels, cmap='jet', alpha=0.9)

    # Add contour lines
    contour_lines = ax.contour(X, Y, Z_smoothed, n_levels, colors='black', linewidths=0.5)
    plt.clabel(contour_lines, inline=True, fontsize=10)

    # Color bar
    cbar = plt.colorbar(contour_filled, ax=ax)
    cbar.set_label('Eficiência (%)', fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    # Axis formatting
    ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='x')
    ax.set_yticks(np.arange(0.1, 1.0, 0.1))
    ax.set_xlabel('Velocidade (rpm)', fontsize=12, color='black', weight='bold')
    ax.set_ylabel('Torque (Nm)', fontsize=12, color='black', weight='bold')

    # Apply font weight to ticks
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontsize(12)
        tick.set_weight('bold')

    plt.tight_layout()
    plt.show()

# -------------------------------
# Main Function
# -------------------------------
def main():
    # Define input data arrays
    X = Speed        # Speed [rpm] -> shape (N, M)
    Y = Tq[:, 0]     # Torque [Nm] -> shape (N,)
    Z = Eff          # Efficiency [%] -> shape (N, M)

    configure_plot()
    create_contour_plot(X, Y, Z)

# -------------------------------
# Run Main
# -------------------------------
if __name__ == "__main__":
    main()
