# Parameters of the simulation

## Loop function:

    This function can be used as value of a scalar parameter (two or more parameter and non-scalar values do not accept this function) to perform a loop over a range of values. Basically, each loop value will be simulated considering all other values fixed. The result will be saved in different directories. There are two types of loops:

        - loop[start_value end_value step]: 
          the value start at "start_value" and goes until "end_value" adding "step" in each iteration;

        - loop{val_1, ..., val_n}: 
          the value start at "val_1" and goes until "val_n";

# Atom (atom.csv) (List)

--- 
- symbol      ->	(char[2])	   Atom symbol
- Z           ->	(int) 		   Atomic number
- mass		    ->	(float)      Mass [Da or u]
--- 

# Transition (transition.csv) (Table)

*Columns*

--- 
- idx         ->  (int)       Identification
- gamma       ->  (float)     Natural linewidth [2pi * kHz]
- lambda      ->  (float)     Resonant wavelength [nm]
- J_gnd       ->  (int)       Total angular momentum of the ground state
- J_exc       ->  (int)       Total angular momentum of the excited state
- g_gnd       ->  (float)     g-factor of the ground state
- g_exc       ->  (float)     g-factor of excited state
--- 

*Rows*

Each row is a transition which can be used.

# Magnetic Field (magnetic_field.csv) (list)

We consider three magnetic fields. The first one is the quadrupole field close to the centre whose axial gradient (scalar) and axial direction (vector) are B_0 and B_axial, respectively. The second one is a constant magnetic field B_bias (vector). The last one is a linear magnetic field B_lin_grad (vector).

--- 
- B_0         ->  (float)     Axial gradient [G / cm]
- B_axial     ->  (float[3])  Axial direction [(x, y, z)]
- B_bias      ->  (float[3])  Constant field (biased) [(x, y, z) G]
- B_lin_grad  ->  (float[3])  Gradient [(x, y, z) G / cm]
--- 

# Beams  (beams/)

## Setup (beams/setup.csv) (table)

*Columns*

--- 
- idx         ->  (int)       Identification
- k_dir       ->  (float[3])  Wave vector direction [(x, y, z)]
- pol_amp     ->  (float[3])  Polarization amplitude [(sigma+, sigma-, pi)]
- s_0         ->  (float)     Saturation parameter (I / I_{sat})
- delta       ->  (float)     Laser detuning [Natural Linewidth]
- trans_id    ->  (int)       Transition identification
- sidebands   ->  (float[2])
--- 


	We are considering two bases to handle with the beams. The first one is 
    the basis C = {c1, c2, c3}, which c1, c2, and c3 are complex vectors 
    related to the polarizations sigma+, sigma-, and pi on the beams frame.
    The second one is the basis D = {d1, d2, d3}, which d1, d2, and d3 are
    also complex vectors related to the polarizations sigma+, sigma-, and 
    pi on the magnetic field frame (d3 is parallel to the magnetic field 
    direction). We are not interested in both real and imaginary components
    of the vectors b and c, we only need  its module, therefore we consider
    the non-unit vector eps = (e1, e2, e3) to define the polarization of 
    each beam in the simulation. This vector is defined on the basis C and
    its components only can be 0 or 1.

- k_dic		-> (float[3])   non-unit vector on the basis A parallel to the wave vector
- eps		-> (int[3])     Polarizations on the basis C

# Initial Conditions (initial_conditions.csv)

- T_0       -> (float)      Initial temperature [uK]
- v_0       -> (float)      Initial velocity [cm/s]
- g_bool    -> (int)        Activate gravity [1: On, 0: Off]

# Environment (conditions.csv)

- B_0       -> (float)      Axial magnetic field gradient [G / cm]
- delta     -> (float)(l)   Laser detuning in units of the transition rate (gamma)
- s_0       -> (float)      Peak of the saturation parameter (I_peak / I_sat)
- w         -> (float)      Waist Radius [cm]
- g_bool    -> (int)        1- use gravity, 0 - do not use gravity