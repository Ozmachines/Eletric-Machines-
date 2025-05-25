
**Summary: Efficiency Map Estimation of Electric Machines**

This work aims to accurately determine the efficiency map of a Permanent Magnet Synchronous Machine (PMSM) using finite element simulations (FEMM) combined with post-processing in Python. The efficiency map represents the machine's performance as a function of electromagnetic torque and rotational speed, enabling evaluation of its operation under various load conditions.

The methodology is structured as follows:

**1. Current Definition and MTPA Control**

   * For each current magnitude (e.g., obtained from MotorCAD), the optimal current angle (beta) that maximizes torque is computed using the MTPA strategy.
   * The corresponding id and iq values are calculated and transformed into phase currents (i\_abc) using the inverse Park transform.

**2. Magnetostatic FEMM Simulation**

   * The machine model is excited with the computed i\_abc currents across a full mechanical revolution (0 to 360 degrees).
   * At each rotor angle:

     * The electromagnetic torque is computed.
     * Magnetic flux density (B) and magnetic vector potential (A) are extracted from each finite element.

**3. Harmonic Analysis via FFT**

   * A Fast Fourier Transform (FFT) is applied to the spatial distribution of B and A to obtain harmonic content.
   * Loss components are derived based on harmonic magnitudes:

     * Iron losses (hysteresis and eddy currents) in stator and rotor cores
     * Proximity effect losses in stator windings
     * Eddy current losses in permanent magnets

**4. Electromagnetic Loss Calculation**

   * Local losses are computed by integrating field quantities over the mesh element volumes.
   * DC ohmic losses in the windings are also calculated using RMS current values and winding resistance.

**5. Efficiency Map Construction**

   * For each combination of current and speed:

     * Input and output powers are determined.
     * Mechanical torque and overall efficiency are computed.
     * Results are stored in a 2D matrix forming the efficiency map.

**Final Outcome**

The resulting efficiency map allows:

* Identification of optimal operating regions
* Quantification of dominant loss mechanisms
* Informed comparison of machine topologies and control strategies
