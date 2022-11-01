/**
 * Monte Carlos simulation of a single atom in a Magneto-Optical Trap (MOT)
 * Bruno N. Santos <nicolau.bruno@gmail.com>
 * Version 2.0
 */

//  Header
#include "mot_sim.h"

results_t simulate_atom(char *params_path, int opt, long seed_time){
    // Variables
    int i, iter = 0;                            // Auxiliary
    double r, v, t = 0, dt = 0;                 // Dynamics
    double s_0_ini, B_0_ini, delta_ini;         // Compression
    double progress;                            // Progress
    int last_progress = 0;                      // Progress
    results_t res;

    // Seed random variable
    srand(time(0) + seed_time);

    // Parameters of the simulation
    initial_conditions_t ini_conds = get_initial_conditions(params_path);
    performance_t perform = get_performance(params_path);
    magnetic_field_t B_params = get_magnetic_field(params_path, perform);
    beams_setup_t beams_setup = get_beams(params_path);
    zeeman_slower_t ZS = get_zeeman_slower(params_path);
    compression_t compression = get_compression(params_path);
    atom_t atom = get_atom(ini_conds, perform, beams_setup, B_params, ZS, opt, params_path);

    // Compression
    B_0_ini = B_params.B_0;
    s_0_ini = beams_setup.beams[0].s_0;
    delta_ini = beams_setup.beams[0].delta;

    // Print parameters
    //printf("opt = %d\n", opt);
    //print_initial_conditions(ini_conds);
    //print_performance(perform);
    //print_magnetic_field(B_params);
    //print_beams(beams_setup);
    //print_zeeman_slower(ZS);
    //print_compression(compression);
    //print_atom(atom);
    //exit(0);    

    atom.pos[0] = 0.0;
    atom.pos[1] = 0.0;
    atom.pos[2] = 0.5;

    atom.vel[0] = 10;
    atom.vel[1] = 10;
    atom.vel[2] = 10;

    // Set initial values
    //--
    res.time = 0;
    res.trapped_atom = 0;
    res.ini_vel = (double*) calloc(3, sizeof(double));
    //for(i = 0; i < 3; i++) 
    //    res.ini_vel[i] = atom.vel[i];

    // Histograms
    //set_hist(opt, &res, perform);
    //--

    // Distance from origin of the magnetic field
    //r = r3_mod(r3_diff(atom.pos, B_params.B_origin));

    // Speed of the atom
    //v = r3_mod(atom.vel);

    //print_results(res, atom, opt, 0);
    //exit(0);

    //
    // Iterations
    //--
    while(res.time < (perform.recording_time + perform.wait_time)){
        // Execute 30% of the simulation time
        if(res.time > (perform.wait_time/3)){
            // Check whether the atom was trapped
            if(res.time > (perform.wait_time/3)){
                if(!is_inside_threshold(atom, B_params, beams_setup, perform) || (v > perform.max_v)){
                    res.escape_square_vel = v*v;
                    break;
                }
            }
        }

        // Move atom
        dt = move(&atom, beams_setup, perform, ini_conds, B_params, ZS);

        // Iterations numbers
        res.time += dt; // 1 / gamma

        // Distance from origin
        r = r3_mod(r3_diff(atom.pos, B_params.B_origin));

        // Speed of the atom
        v = r3_mod(atom.vel);

        // Waiting the equilibrium
        if(res.time > perform.wait_time && res.time < (perform.wait_time + compression.time)){
            t = res.time - perform.wait_time;

            // Compression
            if(t < compression.time){  
                ZS.beam.s_0 = 0;              
                beams_setup.beams[0].s_0 = ((compression.s_c - s_0_ini) * (t/compression.time) + s_0_ini);
                beams_setup.beams[0].delta = ((compression.delta_c - delta_ini) * (t/compression.time) + delta_ini);
                B_params.B_0 = ((compression.B_c - B_0_ini)/compression.time) * t + B_0_ini;
            } 
        } else if(res.time > (perform.wait_time + compression.time))
            res.trapped_atom = 1;

        //
        // Update results
        if(res.trapped_atom == 1 && r < perform.max_r){
            //
            // Update position and velocity
            //--
            // Marginal histogram
            if(opt == 1 || opt == 2){
                for(i = 0; i < 3; i++) {
                    update_hist(&res.pos_hist[i], atom.pos[i] - B_params.B_origin[i]);
                    if(v < perform.max_v) update_hist(&res.vel_hist[i], atom.vel[i]);
                }
            }

            // 3D-Histograms
            else if(opt == 3){
                update_hist_3d(&res.pos_3Dhist, r3_diff(atom.pos, B_params.B_origin));
                if(v < perform.max_v) update_hist_3d(&res.vel_3Dhist, atom.vel);
            }
            //--
        }

        progress = (100*res.time / (perform.recording_time + perform.wait_time + compression.time));
        if((((int) progress) % 2) == 0 && last_progress < ((int) progress)){
            print_status(atom, res, B_params, beams_setup, progress);
            last_progress = (int) progress;
        }

        //print_status(atom, res, B_params, beams_setup, progress);

        // Number of iterations
        iter++;
    }
    //--

    print_status(atom, res, B_params, beams_setup, progress);
    //print_results(res, atom, opt, 1);

    return res;
}

performance_t get_performance(char *params_path){
    //
    // Variables
    //

    int i = 0, j;
    char *token, *saveptr, *path, **rows;
    performance_t perform;

    // Read lines from the CSV file
    path = str_concatenate(params_path, "performance.csv");
    rows = read_lines(path);

    //
    // Maximum time of simulation
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.recording_time = atof(str_replace(token, ".", ","));

    // C Program
    else perform.recording_time = atof(token);
    //--

    //
    // Maximum distance
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.max_r = atof(str_replace(token, ".", ","));

    // C Program
    else perform.max_r = atof(token);
    //--

    //
    // Maximum speed
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.max_v = atof(str_replace(token, ".", ","));

    // C Program
    else perform.max_v = atof(token);
    //--

    // Skip parameter (number of simulations)
    i += 1;

    //
    // Number of bins
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.num_bins = (int) atof(str_replace(token, ".", ","));

    // C Program
    else perform.num_bins = (int) atof(token);
    //--

    //
    // Waiting time (time to reach the equilibrium)
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.wait_time = atof(str_replace(token, ".", ","));

    // C Program
    else perform.wait_time = atof(token);
    //--

    //
    // Time interval
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) perform.dt = atof(str_replace(token, ".", ","));

    // C Program
    else perform.dt = atof(token);
    //--

    // Release memory
    for(j = 0; j < i; j++) free(rows[j]);
    free(rows);
    free(path);

    return perform;
}

initial_conditions_t get_initial_conditions(char *params_path){
    //
    // Variables
    //

    int i = 0, j;
    char *token, *saveptr, *path, **rows;
    initial_conditions_t ini_conds;

    // Read lines from the CSV file
    path = str_concatenate(params_path, "initial_conditions.csv");
    rows = read_lines(path);

    // Initial temperature 
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) ini_conds.T_0 = atof(str_replace(token, ".", ","));

    // C program
    else ini_conds.T_0 = atof(token);
    //--

    // Module of the initial velocity of the atoms
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) ini_conds.v_0 = atof(str_replace(token, ".", ","));

    // C program
    else ini_conds.v_0 = atof(token);
    //--

    //
    // Gravity
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) ini_conds.g_bool = atoi(str_replace(token, ".", ","));

    // C program
    else ini_conds.g_bool = atoi(token);
    //--

    // Release memory
    for(j = 0; j < i; j++) free(rows[j]);
    free(rows);
    free(path);

    return ini_conds;
}

magnetic_field_t get_magnetic_field(char *params_path, performance_t perform){
    // Variables
    int i = 0, n, j;
    double *angles;
    double *B, *r;
    char *token, *saveptr, *path, **rows;
    magnetic_field_t B_params;

    // Read lines from the CSV file
    path = str_concatenate(params_path, "magnetic_field.csv");
    rows = read_lines(path);

    // Magnetic field gradient
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value
    
    // Python module
    if(Py_MODULE) B_params.B_0 = atof(str_replace(token, ".", ","));

    // C program
    else B_params.B_0 = atof(token);
    //--

    // Magnetic field axial direction
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value
    angles = get_double_array(token, &n);

    // Rotating about x
    B_params.B_basis = rotating_matrix(angles[0], 1);

    // Rotation about y
    B_params.B_basis = r3_operator_product(B_params.B_basis, rotating_matrix(angles[1], 2));

    // Rotation about z
    B_params.B_basis = r3_operator_product(B_params.B_basis, rotating_matrix(angles[2], 3));
    //--

    // Bias of magnetic field
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) B_params.B_bias = get_double_array(token, &n);

    // C program
    else B_params.B_bias = get_double_array(token, &n);
    //--

    // Linear magnetic field gradient
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) B_params.B_lin_grad = get_double_array(token, &n);

    // C program
    else B_params.B_lin_grad = get_double_array(token, &n);
    //--

    // Origin of the magnetic field
    //--
    r = (double*) calloc(3, sizeof(double));
    B_params.B_origin = (double*) calloc(3, sizeof(double));

    for(j = 0; j < 3; j++) r[j] = 1;
    B = magnetic_field(B_params, r);

    for(j = 0; j < 3; j++)
        B_params.B_origin[j] = - B_params.B_bias[j] / (B[j] - B_params.B_bias[j]);

    //--

    // Release memory
    for(j = 0; j < i; j++) free(rows[j]);
    free(r);
    free(rows);
    free(path);
    free(angles);

    return B_params;
}

beams_setup_t get_beams(char *params_path){
    // Variables
    int i = 0, num_beams = 0, n, j;
    char *token, *saveptr, *path, **rows;
    beams_setup_t beams_setup;
    double s_0, delta, w;
    beam_t *beams = (beam_t*) calloc(MAX_BEAMS, sizeof(beam_t));
    beam_t *c_beams;

    //
    // Get main parameters of the beams
    //--
    // Read lines from the CSV file
    path = str_concatenate(params_path, "beams/main.csv");
    rows = read_lines(path);
    i = 0;

    //
    // Laser detuning
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE){
        if(token[0] == '-') delta= -atof(str_replace(token+1, ".", ",")); 
        else delta = atof(str_replace(token, ".", ","));
    } 

    // C program
    else {
        if(token[0] == '-') delta= -atof(token+1); 
        else delta = atof(token);
    }
    //--

    //
    // Peak of the saturation parameter
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) s_0 = atof(str_replace(token, ".", ","));

    // C program
    else s_0 = atof(token);
    //--

    //
    // Waist Radius
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) w = atof(str_replace(token, ".", ","));

    // C program
    else w = atof(token);
    //--
    //--

    //
    // Get all beams in the setup
    //--
    // Read lines from the CSV file
    path = str_concatenate(params_path, "beams/setup.csv");
    rows = read_lines(path);

    //
    // Get beams
    //--
    for(i = 1; !(rows[i] == NULL); i++){
        //
        // Wave vector direction
        //--
        token = strtok_r(rows[i], DELIM, &saveptr);
        beams[num_beams].k_dir = r3_normalize(get_double_array(token, &n));
        //--

        //
        // Polarization vector
        //--
        token = strtok_r(NULL, DELIM, &saveptr); // Value
        beams[num_beams].pol_amp = r3_normalize(get_double_array(token, &n));
        //--

        // Add main parameters
        beams[num_beams].s_0 = s_0;
        beams[num_beams].delta = delta;
        beams[num_beams].w = w;

        num_beams++;
    }
    //--

    c_beams = (beam_t *) calloc(num_beams, sizeof(beam_t));
    for(n = 0; n < num_beams; n++) c_beams[n] = beams[n];

    beams_setup.num = num_beams;
    beams_setup.beams = c_beams;
    //--

    //
    // Get sidebands
    //--
    // Read lines from the CSV file
    path = str_concatenate(params_path, "beams/sidebands.csv");
    rows = read_lines(path);
    i = 0;

    //
    // Number of sidebands
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE)
        beams_setup.sidebands.num = atoi(str_replace(token, ".", ","));

    // C program
    else beams_setup.sidebands.num = atoi(token);
    //--

    //
    // Resonant frequency of the sidebands
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) beams_setup.sidebands.freq = atof(str_replace(token, ".", ","));

    // C program
    else beams_setup.sidebands.freq = atof(token);
    //--
    //--

    // Release memory
    for(j = 0; j < i; j++) free(rows[j]);
    free(beams);
    free(path);
    free(rows);

    // Return
    return beams_setup;
}

transition_t get_transition(char *params_path){
    // Variables
    int i = 0, j;
    char *token, *saveptr, *path, **rows;
    transition_t transition;

    // Read lines from the CSV file
    path = str_concatenate(params_path, "transition.csv");
    rows = read_lines(path);

    //
    // Transition rate
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.gamma = atof(str_replace(token, ".", ","));

    // C program
    else transition.gamma = atof(token);
    //--

    //
    // Resonant wave length
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    if(Py_MODULE) transition.lambda = atof(str_replace(token, ".", ","));

    // C program
    else transition.lambda = atof(token);
    //--

    //
    // Total angular momentum of the ground state
    //--    
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.J_gnd = (int) atof(str_replace(token, ".", ","));

    // C Program
    else transition.J_gnd = (int) atof(token);  
    //--

    //
    // Total angular momentum of the excited state
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.J_exc = (int) atof(str_replace(token, ".", ","));

    // C Program
    else transition.J_exc = (int) atof(token);  
    //--

    //
    // Landè factor of the ground state
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.g_gnd = atof(str_replace(token, ".", ","));

    // C Program
    else transition.g_gnd = atof(token);  
    //--

    //
    // Landè factor of the excited state
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) transition.g_exc = atof(str_replace(token, ".", ","));

    // C Program
    else transition.g_exc = atof(token); 
    //--

    //
    // Check values
    //--
    if((transition.J_exc - transition.J_gnd) < 0){
        printf("J_exc must be grater than J_gnd.\n");
        exit(0);
    } else if(transition.J_exc < 0 || transition.J_gnd < 0){
        printf("J_exc and J_gnd must be positive values.\n");
        exit(0);
    }
    //--

    // Release memory
    for(j = 0; j < i; j++) free(rows[j]);
    free(rows);
    free(path);

    // Return
    return transition;
}

atom_t get_atom(initial_conditions_t ini_conds, performance_t perform, beams_setup_t beams_setup, magnetic_field_t B_params, zeeman_slower_t ZS, int opt, char *params_path){
    //
    // Variables
    // --
    int i = 0, j;
    char *path, **rows;
    char *token, *saveptr;
    atom_t atom;
    //--

    // Read lines from the CSV file
    path = str_concatenate(params_path, "atom.csv");
    rows = read_lines(path);

    //
    // Symbol
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    atom.symbol = (char *) malloc(strlen(token) * sizeof(char));
    strcpy(atom.symbol, token);
    //--

    //
    // Atomic number
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) atom.Z = (int) atof(str_replace(token, ".", ","));

    // C program
    else atom.Z = (int) atof(token);
    //--

    //
    // Mass
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) atom.mass = atof(str_replace(token, ".", ","));

    // C program
    else atom.mass = atof(token);
    //--

    // Optical transition
    atom.transition = get_transition(params_path);

    // Initial state
    atom.J = atom.transition.J_gnd;
    atom.mJ = -atom.transition.J_gnd;

    // Initial atom state
    set_ini_atom_state(&atom, beams_setup, B_params, perform, ini_conds, ZS, opt);

    // Release memory
    for(j = 0; j < i; j++) free(rows[j]);
    free(rows);
    free(path);

    return atom;
}

zeeman_slower_t get_zeeman_slower(char *params_path){
    // Variables
    int i = 0, n, j;
    char *path, **rows;
    char *token, *saveptr;
    zeeman_slower_t ZS;

    // Transition
    path = str_concatenate(params_path, "zeeman_slower/");
    ZS.transition = get_transition(path);

    // Beam
    //--
    path = str_concatenate(params_path, "zeeman_slower/beam.csv");
    rows = read_lines(path);

    // Laser detuning
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE){
        if(token[0] == '-') ZS.beam.delta = -atof(str_replace(token+1, ".", ",")); 
        else ZS.beam.delta = atof(str_replace(token, ".", ","));
    } 

    // C program
    else {
        if(token[0] == '-') ZS.beam.delta = -atof(token+1); 
        else ZS.beam.delta = atof(token);
    }
    //--

    // Resonant saturation parameter
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) ZS.beam.s_0 = atof(str_replace(token, ".", ","));

    // C program
    else ZS.beam.s_0 = atof(token);
    //--

    // Waist
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) ZS.beam.w = atof(str_replace(token, ".", ","));

    // C program
    else ZS.beam.w = atof(token);
    //--

    // Wave vector direction
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr);
    token = strtok_r(NULL, DELIM, &saveptr); // Value
    ZS.beam.k_dir = r3_normalize(get_double_array(token, &n));
    //--

    // Polarization
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr);
    token = strtok_r(NULL, DELIM, &saveptr); // Value
    ZS.beam.pol_amp = r3_normalize(get_double_array(token, &n));
    //--
    //--

    // Release memory
    for(j = 0; j < i; j++) free(rows[j]);
    free(rows);
    free(path);

    return ZS;
}

compression_t get_compression(char *params_path){
    // Variables
    int i = 0, j;
    char *token, *saveptr, *path, **rows;
    compression_t compression;

    // Read lines from the CSV file
    path = str_concatenate(params_path, "compression.csv");
    rows = read_lines(path);

    //
    // Compression time
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) compression.time = atof(str_replace(token, ".", ","));

    // C program
    else compression.time = atof(token);
    //--

    //
    // Resonant saturation parameter after compression
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) compression.s_c = atof(str_replace(token, ".", ","));

    // C program
    else compression.s_c = atof(token);
    //--

    //
    // Axial magnetic field gradient after compression
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE) compression.B_c = atof(str_replace(token, ".", ","));

    // C program
    else compression.B_c = atof(token);
    //--

    //
    // Laser detuning after compression
    //--
    i += 1;
    token = strtok_r(rows[i], DELIM, &saveptr); // Variable name
    token = strtok_r(NULL, DELIM, &saveptr); // Value

    // Python module
    if(Py_MODULE){
        if(token[0] == '-') compression.delta_c = -atof(str_replace(token+1, ".", ",")); 
        else compression.delta_c = atof(str_replace(token, ".", ","));
    } 

    // C program
    else {
        if(token[0] == '-') compression.delta_c = -atof(token+1); 
        else compression.delta_c = atof(token);
    }
    //--

    // Release memory
    for(j = 0; j < i; j++) free(rows[j]);
    free(rows);
    free(path);

    // Return
    return compression;    
}

int set_ini_atom_state(atom_t *atom, beams_setup_t beams_setup, magnetic_field_t B_params, performance_t perform, initial_conditions_t ini_conds, zeeman_slower_t ZS, int opt){
    // Varibles
    int i;
    double *rd_v;

    // Atoms start from the origin
    if(opt == 1) {
        // Set position at origin
        atom->pos = (double*) calloc(3, sizeof(double));
        
        // Initial position
        for(i = 0; i < 3; i++)
            atom->pos[i] = B_params.B_origin[i];

        // Initial velocity
        //--        
        atom->vel = (double *) calloc(3, sizeof(double));
        for(i = 0; i < 3; i++)
            atom->vel[i] = random_norm(0, sqrt(10*k_B * ini_conds.T_0 / (atom->mass * u))) * 100; // cm / s
        //--
    } 

    // Atomic Jet
    else if(opt == 2){
        // Initial velocity
        atom->vel = r3_scalar_product(ini_conds.v_0, r3_scalar_product(-1, ZS.beam.k_dir));

        // Initial position
        //--
        atom->pos = r3_scalar_product(beams_setup.beams[0].w, ZS.beam.k_dir);

        // Random vector
        rd_v = r3_cross_product(random_r3_vector(), atom->pos);
        rd_v = r3_scalar_product(random_norm(0, ZS.beam.w), r3_normalize(rd_v));

        atom->pos = r3_sum(atom->pos, rd_v);
        //--
    }

    return 1;
}

/*
int set_ini_atom_state(atom_t *atom, beams_setup_t beams_setup, magnetic_field_t B_params, performance_t perform, initial_conditions_t ini_conds, zeeman_slower_t ZS, int opt){
    // Variables
    int i;
    double R, z_0 = 0, beta, chi;
    double *r_0, *B, std_dev;
    double delta_zeeman, weight;

    // Theoretical equilibrium position
    if(beams_setup.beams[0].delta < 0){
        // Gravity relevance
        R = 1e7 * (h * 2 * PI * atom->transition.gamma);
        R = R / (2 * atom->mass * u * atom->transition.lambda * g);

        // Magnetic field
        r_0 = (double*) calloc(3, sizeof(double));
        r_0[2] = -1;
        B = magnetic_field(B_params, r_0);
        beta = 2 * PI * (mu_B / h) * B[2] * 1e-4 * 1e10;

        // Transition
        chi = atom->transition.g_exc * (- atom->transition.J_gnd - 1);
        chi += atom->transition.g_gnd * atom->transition.J_gnd;

        // Equilibrium position
        z_0 = sqrt((R - 1)*beams_setup.beams[0].s_0 - 1) / 2;
        z_0 = z_0 + beams_setup.beams[0].delta;
        z_0 = ((2 * PI * atom->transition.gamma * 1e3) / (chi * beta)) * z_0; // cm

        // Release memory
        free(r_0);
        free(B);
    }
    
    // Marginal and 3D histograms option
    if(opt < 2){ 
        // Initial position
        //--
        // Check regime and threshold
        if((fabs(beams_setup.beams[0].delta) < sqrt(1 + beams_setup.beams[0].s_0)) || (fabs(z_0) > perform.max_r) || z_0 > 0){
            atom->pos = r3_scalar_product(random_norm(0, perform.max_r / 2), random_r3_vector()); // cm
        } else {
            atom->pos = (double*) calloc(3, sizeof(double));
            atom->pos[2] = z_0;
        } 
        //--

        // Initial velocity
        //--        
        atom->vel = (double *) calloc(3, sizeof(double));
        for(i = 0; i < 3; i++){
            std_dev = sqrt(k_B * ini_conds.T_0 / (atom->mass * u)) * 10; // cm / s
            atom->vel[i] = random_norm(0, std_dev); // cm / s
        }  
        //--

    // Trapped atoms option
    } else if(opt == 2){
        // Initial velocity
        atom->vel = r3_scalar_product(ini_conds.v_0, r3_scalar_product(-1, ZS.beam.k_dir));

        // Initial position
        atom->pos = r3_scalar_product(perform.max_r, ZS.beam.k_dir);

        // Zeeman shift
        delta_zeeman = zeeman_shift(*atom, B_params);
        weight = 1.0*exp(-2 * r3_mod(atom->pos) / perform.max_r);
    
        // Check Zeeman shift 
        while(delta_zeeman > weight*fabs(beams_setup.beams[0].delta)){
            if((r3_mod(atom->pos) - 0.01) < 0){
                atom->pos = (double*) calloc(3, sizeof(double));
                break;
        
            } else {
                atom->pos = r3_scalar_product(r3_mod(atom->pos) - 0.01, ZS.beam.k_dir);
                delta_zeeman = zeeman_shift(*atom, B_params);
                weight = 1.0*exp(-2 * r3_mod(atom->pos) / perform.max_r);
            }
        }
    }
    //-- 

    return 1;
}
*/

int set_hist(int opt, results_t *res, performance_t perform){
    // Variables
    int i, j, k;

    // Only marginals
    //--
    if(opt == 1 || opt == 2){
        res->pos_hist = (histogram_t*) calloc(3, sizeof(histogram_t));
        res->vel_hist = (histogram_t*) calloc(3, sizeof(histogram_t));

        for(i = 0; i < 3; i++){
            // Position
            res->pos_hist[i].num_bins = perform.num_bins;
            res->pos_hist[i].bin_size = 2 * perform.max_r / res->pos_hist[i].num_bins;
            res->pos_hist[i].coord0 = - perform.max_r;
            res->pos_hist[i].freqs = (int*) calloc(res->pos_hist[i].num_bins, sizeof(int));

            // Velocity
            res->vel_hist[i].num_bins = perform.num_bins;
            res->vel_hist[i].bin_size = 2 * perform.max_v / res->vel_hist[i].num_bins;
            res->vel_hist[i].coord0 = - perform.max_v;
            res->vel_hist[i].freqs = (int*) calloc(res->vel_hist[i].num_bins, sizeof(int));
        }
    }
    //--

    // Complete histograms
    //--
    else if(opt == 3){
        // Position
        res->pos_3Dhist.num_bins = (int*) calloc(3, sizeof(int));
        res->pos_3Dhist.bins_size = (double*) calloc(3, sizeof(double));
        res->pos_3Dhist.coord0 = (double*) calloc(3, sizeof(double));

        // Velocity
        res->vel_3Dhist.num_bins = (int*) calloc(3, sizeof(int));
        res->vel_3Dhist.bins_size = (double*) calloc(3, sizeof(double));
        res->vel_3Dhist.coord0 = (double*) calloc(3, sizeof(double));

        for(i = 0; i < 3; i++){
            // Position
            res->pos_3Dhist.num_bins[i] = perform.num_bins;
            res->pos_3Dhist.bins_size[i] = 2 * perform.max_r / res->pos_3Dhist.num_bins[i];
            res->pos_3Dhist.coord0[i] = - perform.max_r;

            // Velocity
            res->vel_3Dhist.num_bins[i] = perform.num_bins;
            res->vel_3Dhist.bins_size[i] = 2 * perform.max_v / res->vel_3Dhist.num_bins[i];
            res->vel_3Dhist.coord0[i] = - perform.max_v;
        }

        //
        // Position
        //--
        res->pos_3Dhist.freqs = (int***) calloc(res->pos_3Dhist.num_bins[0], sizeof(int**));

        for(i = 0; i < res->pos_3Dhist.num_bins[0]; i++){
            res->pos_3Dhist.freqs[i] = (int**) calloc(res->pos_3Dhist.num_bins[1], sizeof(int*));
            for(j = 0; j < res->pos_3Dhist.num_bins[1]; j++){
                res->pos_3Dhist.freqs[i][j] = (int*) calloc(res->pos_3Dhist.num_bins[2], sizeof(int));
                for(k = 0; k < res->pos_3Dhist.num_bins[2]; k++){
                   res->pos_3Dhist.freqs[i][j][k] = 0; 
                }
            }
        }
        //--

        //
        // Velocity
        //--
        res->vel_3Dhist.freqs = (int***) calloc(res->vel_3Dhist.num_bins[0], sizeof(int**));

        for(i = 0; i < res->vel_3Dhist.num_bins[0]; i++){
            res->vel_3Dhist.freqs[i] = (int**) calloc(res->vel_3Dhist.num_bins[1], sizeof(int*));
            for(j = 0; j < res->vel_3Dhist.num_bins[1]; j++){
                res->vel_3Dhist.freqs[i][j] = (int*) calloc(res->vel_3Dhist.num_bins[2], sizeof(int));
                for(k = 0; k < res->vel_3Dhist.num_bins[2]; k++){
                   res->vel_3Dhist.freqs[i][j][k] = 0; 
                }
            }
        }
        //--
    }
    //--

    return 0;
}

double move(atom_t *atom, beams_setup_t beams_setup, performance_t perform, initial_conditions_t ini_conds, magnetic_field_t B_params, zeeman_slower_t ZS){
    //
    // Variables
    int i, j;                           // Auxiliary variables
    double dt = 0;                      // Time interval
    double R_max = -1, *R, *R_aux;      // Scattering rates
    double *probs;                      // Probability to absorb a beam
    int chosen_beam = 0;                // Absorption variables
    double *a_B;                        // Acceleration
    double *rd_v, vel_mod;              // Photonic recoil
    int num_beams = 0;
    beam_t *valid_beams;
    int check_transition = 0;
    double *valid_R;

    // Allocate variables
    probs = (double*) calloc(beams_setup.num + 2, sizeof(double));
    R = (double*) calloc(beams_setup.num + 1, sizeof(double));
    valid_R = (double*) calloc(beams_setup.num + 1, sizeof(double));
    valid_beams = (beam_t*) calloc(beams_setup.num + 1, sizeof(beam_t));

    // Get scattering rates
    //--
    R_aux = get_scatt_rates(beams_setup, B_params, *atom);

    for(i = 0; i < beams_setup.num; i++){
        if(beams_setup.beams[i].s_0 > 0){
            R[i] = R_aux[i];
            if(R_aux[i] > R_max || R_max < 0) R_max = R_aux[i];
            valid_beams[num_beams] = beams_setup.beams[i];
            valid_R[num_beams] = R[i];
            num_beams += 1;
        }
    }

    if(ZS.beam.s_0 > 0){
        R[i] = get_ZS_scatt_rates(ZS, B_params, *atom);
        if(R[i + 1] > R_max || R_max < 0) R_max = R[i];
        valid_beams[num_beams] = ZS.beam;
        valid_R[num_beams] = R[i];
        num_beams += 1;
        check_transition = num_beams;
    }
    //--

    // Time interval
    //--

    // Defining dt based on the parameter perform.dt
    //aux_dt = perform.dt / (2 * PI * atom->transition.gamma*1e3);
    //if(aux_dt > (1 / ((beams_setup.num + 1)*R_max))) dt = 1 / ((beams_setup.num + 1)*R_max);
    //else dt = aux_dt;

    // Defining dt based on the number of beams
    if(num_beams > 6) dt = 1 / ((num_beams)*R_max);
    else dt = 1 / (6*R_max);

    // Check maximum value
    if((dt * 2 * PI * atom->transition.gamma*1e3) > perform.dt)
        dt = perform.dt / (2 * PI * atom->transition.gamma*1e3);
    //--

    //printf("dt [tau] = %.10f\n", dt * (2 * PI * atom->transition.gamma * 1e3));

    // Get probabilities to absorb a beam
    for(i = 0; i < num_beams; i++){
        probs[i + 1] = valid_R[i] * dt;
        probs[0] += probs[i + 1];

        //r3_print(valid_beams[i].k_dir, "k");
        //printf("R[%d] = %.10f\n", i+1, valid_R[i]);
        //printf("probs[%d] = %.10f\n\n", i+1, probs[i + 1]);
    }

    // Probability of the atom does not absorb a beam
    probs[0] = 1 - probs[0];
    //printf("probs[%d] = %.10f\n\n", 0, probs[0]);
    
    // Pick a beam
    chosen_beam = random_pick(probs, num_beams + 1);
    //printf("chosen_beam = %d\n", chosen_beam);

    // Movement
    //--
    // Magnetic acceleration
    a_B = magnetic_acceleration(*atom, B_params);

    // Update position
    // --
    for(i = 0; i < 3; i++) {
        // Previous velocity
        atom->pos[i] += atom->vel[i] * dt;

        // Magnetic acceleration
        atom->pos[i] += (a_B[i] * dt*dt) / 2;
    }

    // Gravity
    if(ini_conds.g_bool)
        atom->pos[2] += -(g * dt*dt) / 2;
    // --

    // Update velocity
    //--
    // Gravitational acceleration
    if(ini_conds.g_bool)
        atom->vel[2] += - g * dt;

    // Magnetic acceleration
    for(i = 0; i < 3; i++)
        atom->vel[i] += a_B[i] * dt;

    //r3_print(atom->vel, "vel");

    // Photonic recoil
    //--
    if(chosen_beam > 0){
        // Random vector
        rd_v = random_r3_vector();

        // Add velocity
        //--
        if(chosen_beam == check_transition)
            vel_mod = 1e4 * h / (ZS.transition.lambda * atom->mass * u); // cm / s
            
        else
            vel_mod = 1e4 * h / (atom->transition.lambda * atom->mass * u); // cm / s

        atom->vel = r3_sum(atom->vel, r3_scalar_product(vel_mod, valid_beams[chosen_beam-1].k_dir)); // cm / s
        atom->vel = r3_sum(atom->vel, r3_scalar_product(vel_mod, rd_v)); // cm / s
        //--

        //r3_print(r3_scalar_product(vel_mod, valid_beams[chosen_beam-1].k_dir), "Absorption [cm/s]");
        //r3_print(r3_scalar_product(vel_mod, rd_v), "Emission [cm/s]");
        //r3_print(atom->vel, "vel");

        // Release memory
        free(rd_v);
    }
    //--
    //--
    
    // Convert time unit
    dt = dt * (2 * PI * atom->transition.gamma*1e3);
    
    // Release memory
    free(probs);
    free(valid_R);
    free(R);
    free(valid_beams);
    
    //exit(0);

    return dt;
}

double *magnetic_field(magnetic_field_t B_params, double *r){
    //
    // Variables
    //
    int i;
    double *B;          // Magnetic field vector

    // Allocate memory for the magnetic field vector
    B = (double*) calloc(3, sizeof(double));

    // Change basis of the position vector
    B = r3_apply_operator(B_params.B_basis, r);

    // Get magnetic field from B_basis
    B[0] *= 0.5;
    B[1] *= 0.5;
    B[2] *= -1.0;
    B = r3_scalar_product(B_params.B_0, B);

    // Get back from the lab basis
    B = r3_apply_operator(r3_transposed_operator(B_params.B_basis), B);

    // Linear magnetic field gradient
    for(i = 0; i < 3; i++)
        B[i] += B_params.B_lin_grad[i]*r[i];

    // Bias magnetic field
    B = r3_sum(B_params.B_bias, B);

    return B;
}

double *magnetic_acceleration(atom_t atom, magnetic_field_t B_params){
    // Variables
    int i;
    double *a_B, *del_B, norm;  // Magnetic field
    double g_lande;             // Transition
    double *r_prime;

    // Magnetic field gradient
    if(r3_mod(atom.pos) > 0){
        del_B = (double*) calloc(3, sizeof(double)); // G / cm
        r_prime = r3_apply_operator(B_params.B_basis, atom.pos); // cm

        // Anti-helmholtz coils (Magnetic field frame)
        norm = sqrt(pow(r_prime[0], 2) + pow(r_prime[1], 2) + 4 * pow(r_prime[2], 2));
        del_B[0] = r_prime[0] / (2 * norm);
        del_B[1] = r_prime[1] / (2 * norm);
        del_B[2] = - 2 * r_prime[2] / norm;
        del_B = r3_scalar_product(B_params.B_0, del_B); 

        // Change basis (Lab frame)
        del_B = r3_apply_operator(r3_transposed_operator(B_params.B_basis), del_B);

        // Add linear magnetic field gradient
        del_B = r3_sum(B_params.B_lin_grad, del_B);

        // Magnetic acceleration
        a_B = (double*) calloc(3, sizeof(double)); // cm / s^2

        // Atom state
        if(atom.J == atom.transition.J_gnd)
            g_lande = atom.transition.g_gnd;
        else g_lande = atom.transition.g_exc;


        for(i = 0; i < 3; i++)
            a_B[i] = - (mu_B * g_lande * atom.mJ * del_B[i] / (atom.mass * u)) * 1e3; // cm / s^2

        // Release memory
        free(del_B);
        free(r_prime);
    } else {
        a_B = (double*) calloc(3, sizeof(double)); // cm / s^2
        for(i = 0; i < 3; i++) a_B[i] = 0;
    }

    return a_B;
}

double *get_scatt_rates(beams_setup_t beams_setup, magnetic_field_t B_params, atom_t atom){
    // Variables          
    int i, j, l, m;                                             // Auxiliary variables
    double *B, *eB;                                             // Magnetic field     
    double *R;                                                  // Scattering rates
    double s, s_0, r;                                           // Saturation parameter
    double *zee_shift, doppler_shift, laser_detuning, delta;    // Detuning
    double lambda, gamma, g_gnd, g_exc;                         // Transition
    int mj_gnd, mj_exc;                                         // Transition
    double **C;                                                 // Basis of the beam frame
    int pol[] = {+1, -1, 0};                                    // All polarizations
    beam_t beam;                                                // Beam     

    // Allocate memory
    R = (double*) calloc(beams_setup.num, sizeof(double));
    zee_shift = (double*) calloc(3, sizeof(double));

    // Transition
    gamma = atom.transition.gamma; // 2pi kHz
    lambda = atom.transition.lambda; // nm

    //
    // Magnetic field
    //--
    B = magnetic_field(B_params, atom.pos);
    if(r3_mod(B) == 0){
        eB = B_params.B_basis[2];
    } else eB = r3_normalize(B);
    //r3_print(B, "B [G]");
    //--

    // Zeeman shift
    for(i = 0; i < 3; i++){
        zee_shift[i] = zeeman_shift(atom, B_params, pol[i]); // 2pi kHz
        //printf("zeeman_shift(%d) [Gamma] = %.10f\n", pol[i], zee_shift[i] / gamma);
    }
    //printf("\n");

    //
    // Check each beam
    //--
    for(i = 0; i < beams_setup.num; i++){
        // Get beam
        beam = beams_setup.beams[i];
        //r3_print(beam.k_dir, "k");

        // Polarizations
        set_polarizations_amplitudes(&beam, eB);
        //r3_print(beam.pol_amp, "pol");
        //printf("\n");

        // Initial Saturation parameter
        //--
        // Basis of the beam frame
        C = orthonormal_basis(r3_normalize(beam.k_dir));

        // Distance from the propagation axis
        r = pow(r3_inner_product(C[0], atom.pos), 2);
        r += pow(r3_inner_product(C[1], atom.pos), 2);
        r = sqrt(r);

        s_0 = beam.s_0;
        s_0 = s_0 * exp(-2 * pow((r / beam.w), 2));
        //--
        
        // Laser detuning    
        laser_detuning = beam.delta*gamma; // 2pi kHz
        //printf("laser_detuning = %.10f\n", beam.delta*gamma);

        // Doppler shift
        doppler_shift = - 1e4 * r3_inner_product(atom.vel, C[2]) / lambda; // 2pi kHz
        //printf("doppler_shift [Gamma] = %.10f\n", doppler_shift);

        // Check all possible transitions
        //--
        // Polarizations
        for(j = 0; j < 3; j++){
            // Saturation parameter considering sidebands
            s = beam.pol_amp[j] * s_0 / (2*beams_setup.sidebands.num + 1);

            // Main beam
            delta = laser_detuning + zee_shift[j] + doppler_shift;
            //printf("zeeman_shift(%d) [Gamma] = %.10f\n", pol[j], zee_shift[j]);
            //printf("delta(%d) [Gamma] = %.10f\n", pol[j], delta / gamma);
            R[i] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

            // Get all scattering rates due to the sidebands
            for(m = 0; m < beams_setup.sidebands.num; m++){
                // Right sideband
                laser_detuning = (beam.delta + (m+1)*beams_setup.sidebands.freq)*gamma; // 2pi kHz
                delta = laser_detuning + zee_shift[j] + doppler_shift;
                R[i] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz

                // Left sideband
                laser_detuning = (beam.delta - (m+1)*beams_setup.sidebands.freq)*gamma; // 2pi kHz
                delta = laser_detuning + zee_shift[j] + doppler_shift;
                R[i] += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz
            }
        }

        // Release memory  
        for(l = 0; l < 3; l++) free(C[l]); 
        free(C);

        //printf("R[%d] = %.10f\n", i, R[i]);
        //printf("\n");
    }
    //--

    // Release memory
    free(B);
    
    return R;    
}

double get_ZS_scatt_rates(zeeman_slower_t ZS, magnetic_field_t B_params, atom_t atom){
    // Variables          
    int j;                                                      // Auxiliary variables
    double *B, *eB;                                             // Magnetic field     
    double R = 0;                                               // Scattering rates
    double s, s_0, r;                                           // Saturation parameter
    double zee_shift, doppler_shift, laser_detuning, delta;     // Detuning
    double lambda, gamma;                                       // Transition
    double **C;                                                 // Basis of the beam frame
    int pol[] = {+1, -1, 0};                                    // All polarizations
    beam_t beam;

    //
    // Magnetic field
    //--
    B = magnetic_field(B_params, atom.pos);
    if(r3_mod(B) == 0) eB = B_params.B_basis[2];
    else eB = r3_normalize(B);
    //--

    // Polarizations
    beam = ZS.beam;
    set_polarizations_amplitudes(&beam, eB);

    // Initial Saturation parameter
    //--
    // Basis of the beam frame
    C = orthonormal_basis(r3_normalize(beam.k_dir));

    // Distance from the propagation axis
    r = pow(r3_inner_product(C[0], atom.pos), 2);
    r += pow(r3_inner_product(C[1], atom.pos), 2);
    r = sqrt(r);

    s_0 = beam.s_0;
    s_0 = s_0 * exp(-2 * pow((r / beam.w), 2));
    //--

    // Transition
    gamma = ZS.transition.gamma; // 2pi kHz
    lambda = ZS.transition.lambda; // nm

    // Doppler shift
    doppler_shift = - 1e4 * r3_inner_product(atom.vel, C[2]) / lambda; // 2pi kHz

    // Laser detuning
    laser_detuning = beam.delta*gamma; // 2pi kHz

    // Check all possible transitions
    //--
    // Polarizations
    for(j = 0; j < 3; j++){
        // Zeeman shift
        zee_shift = zeeman_shift(atom, B_params, pol[j]);

        // Saturation parameter considering sidebands
        s = beam.pol_amp[j] * s_0;

        // Total detuning
        delta = laser_detuning + zee_shift + doppler_shift;

        // Scattering rate
        R += ((2*PI*gamma * 1e3)/2) * s / (1 + s + 4*(delta*delta)/(gamma*gamma)); // Hz
    }
    //--
        
    // Release memory  
    for(j = 0; j < 3; j++) free(C[j]); 
    free(C);
    free(B);

    return R;    
}

/*
int set_polarizations_amplitudes_old(beam_t *beam, double *eB){
    //
    // Variables
    int i, j;
    double **R1, **R2;
    double *pol_amp;
    complex_t **C1, **C2, **A_s_c, **A;
    complex_t *C1_eps, *C2_eps, *C2_s_eps;

    // Bases of the beam frame
    R1 = orthonormal_basis(beam->k_dir);
    C1 = r3_oper_to_c3_oper(R1);

    // Bases of the B frame
    R2 = orthonormal_basis(eB);
    C2 = r3_oper_to_c3_oper(R2);

    //
    // Change-of-Basis matrix of the Spherical basis to the Cartesian basis
    //--
    A_s_c = c3_operator_zeros();

    A_s_c[0][0].re = - 1 / sqrt(2);
    A_s_c[0][1].re = 1 / sqrt(2);

    A_s_c[1][0].im = 1 / sqrt(2);
    A_s_c[1][1].im = 1 / sqrt(2);

    A_s_c[2][2].re = 1;
    //--

    //
    // Change-of-Basis matrix of the Cartesian beam frame to the Cartesian B frame
    //--
    A = c3_operator_zeros();

    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            A[i][j] = c3_inner_product(C2[i], C1[j]);
        }
    }
    //--

    //
    // Get polarization amplitudes
    //--
    // Polarization on the Cartesian beam frame
    C1_eps = c3_apply_operator(A_s_c, r3_to_c3(beam->pol_amp));

    // Polarization on the Cartesian B frame
    C2_eps = c3_apply_operator(A, C1_eps);

    // Polarization on the Spherical B frame
    C2_s_eps = c3_apply_operator(c3_operator_dagger(A_s_c), C2_eps);

    // Set Polarization amplitudes
    pol_amp = (double*) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++) pol_amp[i] = pow(c_mod(C2_s_eps[i]), 2);
    beam->pol_amp = pol_amp;
    //--

    //
    // Release memory    
    //--
    for(i = 0; i < 3; i++){
        free(A[i]);
        free(A_s_c[i]);
        free(C1[i]);
        free(C2[i]);
        free(R1[i]);
        free(R2[i]);
    }

    free(A);
    free(A_s_c);
    free(C1);
    free(C1_eps);
    free(C2);
    free(C2_eps);
    free(C2_s_eps);
    free(R1);
    free(R2);
    //--

    return 1;
}
*/

int set_polarizations_amplitudes(beam_t *beam, double *eB){
    // Variables
    int i;
    double **R3_basis;
    complex_t **M_P, **M_P_dagger, **M_L, **M_B;
    complex_t *pol_amp;

    // Polarization amplitudes
    pol_amp = r3_to_c3(beam->pol_amp);

    // Change-of-basis matrix from polarization frame P to cartesian frame C
    M_P = change_of_basis_from_P_to_C();
    M_P_dagger = c3_operator_dagger(M_P);

    // Change-of-basis matrix from cartesian beam frame C_k to cartesian lab frame C_L
    R3_basis = orthonormal_basis(beam->k_dir);
    M_L = c3_operator_dagger(r3_oper_to_c3_oper(R3_basis));

    // Change-of-basis matrix from cartesian beam frame C_L to cartesian lab frame C_B
    R3_basis = orthonormal_basis(eB);
    M_B = r3_oper_to_c3_oper(R3_basis);

    // Change bases
    pol_amp = c3_apply_operator(M_P, pol_amp);
    pol_amp = c3_apply_operator(M_L, pol_amp);
    pol_amp = c3_apply_operator(M_B, pol_amp);
    pol_amp = c3_apply_operator(M_P_dagger, pol_amp);

    // Get polarization amplitudes
    for(i = 0; i < 3; i++)
        beam->pol_amp[i] = pow(c_mod(pol_amp[i]), 2);


    // Free memory
    for(i = 0; i < 3; i++){
        free(M_P[i]);
        free(M_P_dagger[i]);
        free(M_L[i]);
        free(M_B[i]);
    }

    free(M_P);
    free(M_P_dagger);
    free(M_L);
    free(M_B);

    return 1;
}

double zeeman_shift(atom_t atom, magnetic_field_t B_params, int pol){
    // Variables
    double chi, g_gnd, g_exc, mj_gnd, mj_exc;       // Transition
    double *B, zeeman_shift;                        // Zeeman shift

    // Transition
    g_gnd = atom.transition.g_gnd;
    g_exc = atom.transition.g_exc;
    mj_gnd = -atom.transition.J_gnd;
    mj_exc = mj_gnd + pol;
    chi = g_gnd * mj_gnd - g_exc * mj_exc;

    // Magnetic field [Gauss]
    B = magnetic_field(B_params, atom.pos);

    // Zeeman shift (units of Gamma)
    zeeman_shift = 1e3 * (mu_B / h) * r3_mod(B) * chi; // 2pi kHz

    return zeeman_shift;
}

int is_inside_threshold(atom_t atom, magnetic_field_t B_params, beams_setup_t beams_setup, performance_t perform){
    // Variables
    //int i, in = 1;
    //double s_square;
    int i, in = 0;

    //for(i = 0; i < beams_setup.num; i++){
    //    s_square = r3_inner_product(atom.pos, atom.pos) - pow(r3_inner_product(atom.pos, beams_setup.beams[i].k_dir), 2);
    //
    //    if(s_square > pow(beams_setup.beams[0].w, 2)){
    //        in = 0;
    //        break;
    //    }
    //}

    for(i = -1; i < 2; i++)
        if((zeeman_shift(atom, B_params, i) - beams_setup.beams[0].delta) < 0) 
            in = 1;

    //if(in && (fabs(zeeman_shift(atom, B_params)) > fabs()))
    //    in = 0;

    return in;
}

complex_t **change_of_basis_from_P_to_C(){
    // Variables
    int i, j;
    double N = 1 / sqrt(2);
    complex_t **M_P;

    // Get 3x3 complex matrix of zeros
    M_P = c3_operator_zeros();

    // Rows
    //--
    M_P[0][0].re = 1;
    M_P[0][1].re = 1;

    M_P[1][0].im = 1;
    M_P[1][1].im = -1;

    M_P[2][2].re = sqrt(2);
    //--

    // Multiply matrix by N
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++){
            M_P[i][j].re = N * M_P[i][j].re;
            M_P[i][j].im = N * M_P[i][j].im;
        }
    }

    return M_P;
}

//
// Utility functions
//

char *str_concatenate(char *str1, char *str2){
    // Variables    
    int i, j, size;
    char *str;

    //
    // Concatenation
    //

    size = (int) (strlen(str1) + strlen(str2) - 1);
    str = (char*) calloc(STRING_BUFFER_SIZE, sizeof(char));
    
    for(i = 0; i < ((int) strlen(str1)); i++)
        str[i] = str1[i];

    for(j = 0; j < ((int) strlen(str2)); j++)
        str[i + j] = str2[j];

    str[size+1] = '\0';

    return str;
}

char **read_lines(char *path){
    //
    // Variables
    int i = 0, j = 0, k = 0;
    FILE *fp;
    char *row, *aux_row;
    char **rows;

    // Alocate memory
    rows = (char**) calloc(MAX_LINES, sizeof(char*));

    // Open file
    if((fp = fopen(path, "r"))){
        // Read lines
        while(!feof(fp)){
            // Allocate memory
            aux_row = (char*) calloc(STRING_BUFFER_SIZE, sizeof(char));
            row = (char*) calloc(STRING_BUFFER_SIZE, sizeof(char));

            // Try to read line
            if(fgets(row, STRING_BUFFER_SIZE, fp) != NULL){
                // Get row length and remove \n
                k = 0;
                for(j = 0; row[j] != '\0'; j++) {
                    if(row[j] != '\n'){
                        aux_row[k] = row[j];
                        k += 1;
                    }
                }

                // Copy the read row
                rows[i] = (char*) calloc(k+1, sizeof(char));
                strcpy(rows[i], aux_row);

                i += 1;
            } else break;

            // Release memory
            free(row);
            free(aux_row);
        }

        // Indicate the file ending
        rows[i+1] = (char*) calloc(3, sizeof(char));
        rows[i+1] = NULL;

        // Close file
        fclose(fp); 
    } else {
        printf("File \"%s\" does not exist\n", path);
        exit(0);
    }

    return rows;
}

char *str_replace(char *orig, char *rep, char *with){
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep (the string to remove)
    int len_with; // length of with (the string to replace rep with)
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    // sanity checks and initialization
    if (!orig || !rep)
        return NULL;
    len_rep = strlen(rep);
    if (len_rep == 0)
        return NULL; // empty rep causes infinite loop during count
    if (!with)
        with = "";
    len_with = strlen(with);

    // count the number of replacements needed
    ins = orig;
    for (count = 0; (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }

    tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    // first time through the loop, all the variable are set correctly
    // from here on,
    //    tmp points to the end of the result string
    //    ins points to the next occurrence of rep in orig
    //    orig points to the remainder of orig after "end of rep"
    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);

    return result;
}

double *get_double_array(char *str, int *size){
    //
    // Variables
    //
    int i, j, max_size = 124;
    char *token;
    double aux_arr[max_size];
    double *arr;

    str = str + 1;
    str[strlen(str)-1] = '\0';
    token = strtok_r(str, " ", &str);    

    //
    // Parse string
    //

    i = 0;
    while(token && (i < max_size)){
        // Py MODULE
        if(Py_MODULE){
            if(token[0] == '-') aux_arr[i] = -atof(str_replace(token+1, ".", ",")); 
            else aux_arr[i] = atof(str_replace(token, ".", ","));
        }

        // C Program
        else {
            if(token[0] == '-') aux_arr[i] = -atof(token+1); 
            else aux_arr[i] = atof(token);
        }

        token = strtok_r(str, " ", &str);
        i++;
    }

    arr = (double *) malloc(i * sizeof(double));
    for(j = 0; j < i; j++) arr[j] = aux_arr[j];

    // Release memory
    free(token);

    return arr;
}

double **orthonormal_basis(double *v){
    //
    // Variables
    //

    int i;
    double **B;             // Desired basis
    double *v1, *v2, *v3;   // Auxiliary vectors

    // Normalize vector v
    v3 = r3_normalize(v);

    // Generate a random vector  
    v1 = (double*) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++) v1[i] = ((double) rand()) / ((double) RAND_MAX);

    // Define a orthonormal vector
    v2 = r3_scalar_product(r3_inner_product(v1, v3), v3);
    v1 = r3_normalize(r3_diff(v1, v2));

    // Define the last vector of the basis
    v2 = r3_cross_product(v3, v1);

    // Define basis
    B = (double **) calloc(3, sizeof(double *));

    for(i = 0; i < 3; i++)
        B[i] = (double *) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++){
        B[0][i] = v1[i];
        B[1][i] = v2[i];
        B[2][i] = v3[i];
    }

    // Free
    free(v1);
    free(v2);
    free(v3);

    // Return
    return B;
}

double random_exp(double mean){
    return (- mean * log(1 - ((double) rand()) / ((double) RAND_MAX)));
}

double random_norm(double mean, double std_dev){
    //
    // Variables
    int i;
    double *v, r, theta;     // Variables for Box-Muller method
    double std_norm;         // Normal(0, 1)
    double norm;             // Adjusted normal

    v = (double*) calloc(2, sizeof(double));

    //
    // Box-Muller transform
    //--
    // Generate uniform random numbers
    for(i = 0; i < 2; i++) v[i] = ((double) rand()) / ((double) RAND_MAX);

    // Compute r
    r = sqrt(-2 * log(v[0]));

    // Generate theta
    theta = 0.0;
    while(theta == 0.0) theta = 2.0 * PI * v[1];

    // Generate std_norm value
    std_norm = r * cos(theta);

    // Adjust std_norm
    norm = (std_norm * std_dev) + mean;
    //--

    // Release memory
    free(v);

    return norm;
}

double random_uniform(double a, double b){
    return (a + (b - a) * ((double) rand()) / ((double) RAND_MAX));
}

int random_pick(double *probs, int size){
    //
    // Variables
    int i, picked = 0;
    double module;
    double rd_n;          // Random number
    double *cum_probs;    // Cumulative probabilities
    int *idx;

    // Indexes
    idx = (int*) calloc(size, sizeof(int));
    for(i = 0; i < size; i++) idx[i] = i;

    // Normalize probabilities
    module = 0;
    for(i = 0; i < size; i++) module += probs[i];
    for(i = 0; i < size; i++) probs[i] = 100*probs[i] / module;

    // Generate a random number  
    rd_n = 100*((double) rand()) / ((double) RAND_MAX);

    // Cumulative probabilities
    cum_probs = (double*) calloc(size, sizeof(double));
    cum_probs[0] = probs[0];
    for(i = 1; i < size; i++){
        cum_probs[i] += cum_probs[i-1] + probs[i];
    }

    // Pick
    for(i = 0; i < size; i++){
        if(rd_n < cum_probs[i]){
            picked = idx[i];
            break;
        }
    }

    // Release memory
    free(cum_probs);
    free(idx);

    // Return
    return picked;
}

double *random_r3_vector(){
    // Variables
    double *v;
    int i;

    v = (double*) calloc(3, sizeof(double));

    for(i = 0; i < 3; i++) {
        v[i] = ((double) rand()) / ((double) RAND_MAX);
        if((((double) rand()) / ((double) RAND_MAX)) < 0.5) v[i] = -v[i];
    }

    return r3_normalize(v);
}

int update_hist(histogram_t *hist, double val){
    //
    // Variables
    //

    int bin;
    double lower_lim, upper_lim;

    // Add frequency
    for(bin = 0; bin < hist->num_bins; bin++){
        lower_lim = hist->coord0 + bin * hist->bin_size;
        upper_lim = lower_lim + hist->bin_size;

        if((val >= lower_lim) && (val < upper_lim)){
            hist->freqs[bin] += 1;
            break;
        }
    }

    return 1;
}

int update_hist_3d(histogram_3d_t *hist, double *vals){
    //
    // Variables
    //

    int i, bin, *coord, dim = 3;
    double lower_lim, upper_lim;

    // Coordinates
    coord = (int*) calloc(dim, sizeof(int));

    // Check each coordinate
    for(i = 0; i < dim; i++){
        // Add frequency
        for(bin = 0; bin < (*hist).num_bins[i]; bin++){
            lower_lim = (*hist).coord0[i] + bin * (*hist).bins_size[i];
            upper_lim = lower_lim + (*hist).bins_size[i];

            if((vals[i] >= lower_lim) && (vals[i] < upper_lim)){
                coord[i] = bin;
                break;
            }
        }  
    }

    (*hist).freqs[coord[0]][coord[1]][coord[2]] += 1;

    // Release memory
    free(coord);

    return 1;
}

double **rotating_matrix(double theta, int axis){
    // Variables
    int i;       // Utility
    double **R;

    // Allocate memory
    R = (double**) calloc(3, sizeof(double*));
    for(i = 0; i < 3; i++) R[i] = (double*) calloc(3, sizeof(double));

    // Convert to radian
    theta = theta * PI / 180;

    // Rotating about x
    if(axis == 1){
        R[0][0] = 1;

        R[1][1] = cos(theta);
        R[2][2] = cos(theta);
        R[1][2] = -sin(theta);
        R[2][1] = sin(theta);

    // Rotating about y
    } else if(axis == 2){
        R[1][1] = 1;

        R[0][0] = cos(theta);
        R[2][2] = cos(theta);
        R[0][2] = sin(theta);
        R[2][0] = -sin(theta);

    // Rotating about z
    } else if(axis == 3){
        R[2][2] = 1;

        R[0][0] = cos(theta);
        R[1][1] = cos(theta);
        R[0][1] = -sin(theta);
        R[1][0] = sin(theta);
    }

    return R;
}

//
// Debug
//

int print_performance(performance_t perform){
    printf("Performance\n--\n");
    printf("recording_time [1/Gamma] = %f\n", perform.recording_time);
    printf("max_r [cm] = %f\n", perform.max_r);
    printf("max_v [cm/s] = %f\n", perform.max_v);
    printf("num_bins = %d\n", perform.num_bins);
    printf("wait_time [1/Gamma] = %f\n", perform.wait_time);
    printf("time interval [1/Gamma] = %f\n", perform.dt);
    printf("\n");

    return 1;
}

int print_initial_conditions(initial_conditions_t ini_conds){
    printf("Initial Conditions\n--\n");
    printf("T_0 [uK] = %f\n", ini_conds.T_0);
    printf("v_0 [m/s] = %f\n", ini_conds.v_0);
    printf("g_bool = %d\n", ini_conds.g_bool);
    printf("\n");

    return 1;
}

int print_magnetic_field(magnetic_field_t B_params){
    printf("Magnetic field\n--\n");
    printf("B_0 [G/cm] = %f\n", B_params.B_0);
    r3_operator_print(B_params.B_basis, "B_basis");
    r3_print(B_params.B_bias, "B_bias [G]");
    r3_print(B_params.B_lin_grad, "B_lin_grad [G/cm]");
    r3_print(B_params.B_origin, "B_origin [G]");
    printf("\n");

    return 1;
}

int print_atom(atom_t atom){
    printf("Atom\n--\n");
    printf("symbol = %s\n", atom.symbol);
    printf("Z = %d\n", atom.Z);
    printf("mass [u] = %f\n", atom.mass);
    r3_print(atom.pos, "pos [cm/s]");
    printf("r [cm] = %f\n", r3_mod(atom.pos));
    r3_print(atom.vel, "vel [cm/s]");
    printf("v [cm/s] = %f\n", r3_mod(atom.vel));
    printf("J = %d\n", atom.J);
    printf("mJ = %d\n", atom.mJ);
    printf("\n");

    return 1;
}

int print_beams(beams_setup_t beams_setup){
    // Variables
    int i;

    printf("\nSidebands\n");
    printf("num = %d\n", beams_setup.sidebands.num);
    printf("freq = %f\n", beams_setup.sidebands.freq);
    printf("\n");
    printf("All Beams\n--\n");
    for(i = 0; i < beams_setup.num; i++){
        printf("Beam %d\n", i+1);
        printf("s_0 = %f\n", beams_setup.beams[i].s_0);
        printf("delta [1/Gamma] = %f\n", beams_setup.beams[i].delta);
        printf("w [cm] = %f\n", beams_setup.beams[i].w);
        r3_print(beams_setup.beams[i].k_dir, "k");
        r3_print(beams_setup.beams[i].pol_amp, "pol_amp");
        printf("--\n");
    }
    printf("\n");

    return 1;
}

int print_compression(compression_t compression){
    printf("Compression\n--\n");
    printf("t [tau] = %f\n", compression.time);
    printf("s_c = %f\n", compression.s_c);
    printf("B_c [G/cm] = %f\n", compression.B_c);
    printf("\n");

    return 1;

}

int print_results(results_t res, atom_t atom, int opt, int show_hist){
    // Variables
    int i, j;

    printf("Simulation status\n--\n");
    r3_print(res.ini_vel, "Initial Velocity [cm/s]");
    printf("total time [ms] = %f\n", res.time / (atom.transition.gamma));
    printf("total time [1/Gamma] = %f\n", res.time);
    printf("Atom trapped [ms] = %d\n", res.trapped_atom);
    printf("\n");

    if(show_hist && (opt == 1 || opt == 2)){
        //
        // Position frequencies
        //-
        for(i = 0; i < 3;i++){
            printf("pos[%d] = [\n", i+1);
            for(j = 0; j < res.pos_hist[i].num_bins; j++)
                printf("%d ", res.pos_hist[i].freqs[j]);
            printf("]\n\n");
        }
        printf("\n");
        //--

        //
        // Velocity frequencies
        //--
        for(i = 0; i < 3;i++){
            printf("vel[%d] = [\n", i+1);
            for(j = 0; j < res.vel_hist[i].num_bins; j++)
                printf("%d ", res.vel_hist[i].freqs[j]);
            printf("]\n\n");
        }
        printf("\n");
        //--
    }

    return 1;
}

int print_status(atom_t atom, results_t res, magnetic_field_t B_params, beams_setup_t beams_setup, double progress){
    printf("Simulation status (%f)\n--\n", progress);
    printf("Time = %f ms\n", 1e3*res.time / (2*PI*1e3 * atom.transition.gamma));
    printf("Time = %f tau\n", res.time);
    r3_print(r3_diff(atom.pos, B_params.B_origin), "atom position [cm]");
    r3_print(atom.vel, "atom velocity [cm / s]");
    printf("distance from origin = %f cm\n", r3_mod(r3_diff(atom.pos, B_params.B_origin)));
    printf("atom speed = %f cm/s\n", r3_mod(atom.vel));
    printf("trapped_atom = %d\n", res.trapped_atom);
    printf("B_0 = %f\n", B_params.B_0);
    printf("s_0 = %f\n", beams_setup.beams[0].s_0);
    printf("delta = %f\n", beams_setup.beams[0].delta);
    printf("\n");

    return 1;
}

int print_zeeman_slower(zeeman_slower_t ZS){
    printf("Zeeman Slower\n--\n");

    printf("\nBeam\n");
    printf("s_0 = %f\n", ZS.beam.s_0);
    printf("delta [tau] = %f\n", ZS.beam.delta);
    printf("w [cm] = %f\n", ZS.beam.w);
    r3_print(ZS.beam.k_dir, "k");
    r3_print(ZS.beam.pol_amp, "pol_amp");

    printf("\n");

    printf("Transition\n");
    printf("gamma [2pi kHz] = %f\n", ZS.transition.gamma);
    printf("lambda [nm] = %f\n", ZS.transition.lambda);
    printf("J_gnd = %d\n", ZS.transition.J_gnd);
    printf("J_exc = %d\n", ZS.transition.J_exc);
    printf("g_gnd = %f\n", ZS.transition.g_gnd);
    printf("g_exc = %f\n", ZS.transition.g_exc);

    return 1;
}