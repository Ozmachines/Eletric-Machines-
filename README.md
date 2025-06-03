**Summary: Efficiency Map Estimation of Electric Machines**

This work aims to determine the efficiency map of a Permanent Magnet Synchronous Machine (PMSM) using finite element simulations (FEMM) combined with post-processing in Python. The efficiency map represents the machine's performance as a function of electromagnetic torque and rotational speed, enabling evaluation of its operation under various load conditions.

The methodology is structured as follows:

**1. Current Definition and Torque Maximization**

* For each current magnitude the current angle (beta) that maximizes the electromagnetic torque is determined.
* The corresponding *i<sub>d</sub>* and *i<sub>q</sub>* values are computed and transformed into phase currents (*i<sub>abc</sub>*) using the inverse Park transform.

**2. Magnetostatic FEMM Simulation**

* The machine model is excited with the computed *i<sub>abc</sub>* currents across a full mechanical revolution (0 to 360 degrees).
* At each rotor angle:

  * The electromagnetic torque is computed.
  * Magnetic flux density (*B*) and magnetic vector potential (*A*) are extracted from each finite element.

**3. Harmonic Analysis via FFT**

* A Fast Fourier Transform (FFT) is applied to the spatial distribution of *B* and *A* to obtain harmonic content.
* Loss components are derived based on harmonic magnitudes:

  * Iron losses (hysteresis and eddy currents) in stator and rotor cores;
  * Proximity effect losses in stator windings;
  * Eddy current losses in permanent magnets.

**4. Electromagnetic Loss Calculation**

* Local losses are computed by integrating field quantities over the mesh element volumes.
* DC ohmic losses in the windings are also calculated using RMS current values and winding resistance.

**5. Efficiency Map Construction**

* For each combination of torque and speed:

  * Input and output powers are determined.
  * Eletromagnetic losses and overall efficiency are computed.
  * Results are stored in a 2D matrix forming the efficiency map.


