//  Header
#include "mot_sim.h"

// Main
//--
results_t simulate_atom(char *params_path, int opt, long seed_time){
    // Variables
    //--
    // Auxiliaries
    int i, iter = 0, max_iter = 1e6;
    int progress, last_progress = 0;
    double dt, v, r;
    scattering_t scatt;
    
    // Parameters
    initial_conditions_t ini_conds;
    performance_t perform;
    magnetic_field_t B_params;
    transitions_set_t trans_set;
    beams_setup_t beams_setup;
    atom_t atom;

    results_t res;
    //--

    // Seed random variable
    srand(time(0) + seed_time);

    // Parameters
    //--
    set_initial_conditions(params_path, &ini_conds);
    set_performance(params_path, &perform);
    set_magnetic_field(params_path, &B_params);
    set_transitions(params_path, &trans_set);
    set_beams_setup(params_path, &beams_setup);
    set_atom(params_path, &atom);

    // (Debug) Print parameters
    //print_initial_conditions(ini_conds);
    //print_performance(perform);
    //print_magnetic_field(B_params);
    //print_transitions(trans_set);
    //print_beams_setup(beams_setup);
    //print_atom(atom);

    // Set initial state of the atom
    set_ini_atom_state(&atom, ini_conds, B_params, beams_setup, trans_set.trans[0]);
    //r3_print(atom.pos, "pos [cm]");
    //r3_print(atom.vel, "vel [cm/s]");
    //exit(0);
    //--

    // Set results
    res.time = 0;
    res.trapped_atom = 0;
    res.escape_vel = 0;
    res.ini_vel = (double*) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++) res.ini_vel[i] = atom.vel[i];

    // Set histogram
    set_hist(&res, perform);
    //print_results(res, atom, 1);
    //exit(0);

    // Time evolution
    while(res.time < (perform.recording_time + perform.wait_time) && (iter < max_iter)){
        // Update position, velocity and time
        scatt = scattering_event(atom, trans_set, beams_setup, B_params, perform.dt);
        dt = move(&atom, trans_set, scatt, B_params, ini_conds.g_bool); // ms
        res.time += dt; // ms

        // Distance from origin
        r = r3_mod(r3_diff(atom.pos, B_params.B_origin));

        // Speed
        v = r3_mod(atom.vel);

        // Check if the atom pass through the threshold after 30% of the waiting time
        if(res.time > (perform.wait_time / 3)){
            // Check whether the atom was trapped
            if(!is_inside_threshold(atom, trans_set, B_params, beams_setup)){
                res.escape_vel = v;
                break;
            }
        }

        // Equilibrium
        if(res.time > perform.wait_time){
            res.trapped_atom = 1;

            if(r < perform.max_r){
                for(i = 0; i < 3; i++) {
                    update_hist(&res.pos_hist[i], atom.pos[i] - B_params.B_origin[i]);
                    if(v < perform.max_v) update_hist(&res.vel_hist[i], atom.vel[i]);
                }
            }
        }

        // Print progress
        progress = (int) (100 * res.time / (perform.recording_time + perform.wait_time));
        if((progress % 1 == 0) && (progress > last_progress)){
            //print_progress(atom, res, B_params, progress);
            last_progress = progress;
        }

        // Iterations
        iter++;
    }

    //print_progress(atom, res, B_params, progress);
    //print_results(res, atom, 1);

    return res;
}

void set_initial_conditions(char *params_path, initial_conditions_t *ini_conds){
    // Variables
    int i = 0, j = 0;;
    char *path, *token, *saveptr, **rows;

    // File path
    path = str_concatenate(params_path, "initial_conditions.csv");

    // Read lines
    rows = read_lines(path);

    // Check lines
    for(i = 1; !(rows[i] == NULL); i++){
        // Variable
        token = strtok_r(rows[i], DELIM, &saveptr);
        
        // Initial temperature
        //--
        if(strcmp(token, "T_0") == 0){
            token = strtok_r(NULL, DELIM, &saveptr); // Value
            ini_conds->T_0 = str_to_double(token);
        }
        //--
        
        // Initial velocity
        //--
        if(strcmp(token, "v_0") == 0){
            token = strtok_r(NULL, DELIM, &saveptr); // Value
            ini_conds->v_0 = str_to_double(token);
        }
        //--
        
        // Gravity
        //--
        if(strcmp(token, "g_bool") == 0){
            token = strtok_r(NULL, DELIM, &saveptr); // Value
            ini_conds->g_bool = str_to_int(token);
        }
        //--
    }

    // Free memory
    for(j = 0; j < i; j++)
        free(rows[j]);

    free(path);
    free(rows);
}

void set_performance(char *params_path, performance_t *perform){
    // Variables
    int i = 0, j = 0;
    char *path, *token, *saveptr, **rows;

    // File path
    path = str_concatenate(params_path, "performance.csv");

    // Read lines
    rows = read_lines(path);

    // Check lines
    for(i = 1; !(rows[i] == NULL); i++){
        // Variable
        token = strtok_r(rows[i], DELIM, &saveptr);
        
        // Time interval
        //--
        if(strcmp(token, "dt") == 0){
            token = strtok_r(NULL, DELIM, &saveptr); // Value
            perform->dt = str_to_double(token);
        }
        //--
        
        // Waiting time
        //--
        if(strcmp(token, "wait_time") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            perform->wait_time = str_to_double(token);
        }
        //--
        
        // Recording time
        //--
        if(strcmp(token, "recording_time") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            perform->recording_time = str_to_double(token);
        }
        //--
        
        // Maximum distance (resolution of the position grid)
        //--
        if(strcmp(token, "max_r") == 0){
            token = strtok_r(NULL, DELIM, &saveptr); // Value
            perform->max_r = str_to_double(token);
        }
        //--
        
        // Maximum velocity (resolution of the velocity grid)
        //--
        if(strcmp(token, "max_v") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            perform->max_v = str_to_double(token);
        }
        //--
        
        // Number of bins
        //--
        if(strcmp(token, "num_bins") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            perform->num_bins = str_to_int(token);
        }
        //--
    }

    // Free memory
    for(j = 0; j < i; j++)
        free(rows[j]);

    free(path);
    free(rows);
}

void set_magnetic_field(char *params_path, magnetic_field_t *B_params){
    // Variables
    int i = 0, j = 0, n;
    char *path, *token, *saveptr, **rows;
    double *r, *B;

    // Allocate memory
    B_params->B_basis = (double**) calloc(3, sizeof(double*));

    // File path
    path = str_concatenate(params_path, "magnetic_field.csv");

    // Read lines
    rows = read_lines(path);

    // Check lines
    for(i = 1; !(rows[i] == NULL); i++){
        // Variable
        token = strtok_r(rows[i], DELIM, &saveptr);
        
        // Magnetic field gradient
        //--
        if(strcmp(token, "B_0") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            B_params->B_0 = str_to_double(token);
        }
        //--
        
        // Magnetic field basis
        //--
        if(strcmp(token, "B_axial") == 0){
            // Value
            token = strtok_r(NULL, DELIM, &saveptr); // Value
            B_params->B_basis = orthonormal_basis(get_double_array(token, &n));
        }
        //--
        
        // Magnetic field basis
        //--
        if(strcmp(token, "B_bias") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            B_params->B_bias = get_double_array(token, &n);
        }
        //--
        
        // Magnetic field linear gradient
        //--
        if(strcmp(token, "B_lin_grad") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            B_params->B_lin_grad = get_double_array(token, &n);
        }
        //--
    }

    // Origin of the magnetic field
    //--
    r = (double*) calloc(3, sizeof(double));
    B_params->B_origin = (double*) calloc(3, sizeof(double));

    for(j = 0; j < 3; j++) r[j] = 1.0;
    B = magnetic_field(*B_params, r);

    for(j = 0; j < 3; j++)
        B_params->B_origin[j] = - B_params->B_bias[j] / (B[j] - B_params->B_bias[j]);
    //--

    // Free memory
    for(j = 0; j < i; j++)
        free(rows[j]);

    free(path);
    free(rows);
    free(r);
    free(B);
}

void set_transitions(char *params_path, transitions_set_t *trans_set){
    // Variables
    int i = 0, j = 0;
    char *path, *token, *saveptr, **rows;
    transition_t *trans;

    // Allocate memory
    trans = (transition_t*) calloc(MAX_SIZE_ARRAY, sizeof(transition_t));

    // File path
    path = str_concatenate(params_path, "transitions.csv");

    // Read lines
    rows = read_lines(path);

    // Check lines
    for(i = 1; !(rows[i] == NULL); i++){
        // Index
        token = strtok_r(rows[i], DELIM, &saveptr);
        trans[i-1].idx = str_to_int(token);

        // Line width
        token = strtok_r(NULL, DELIM, &saveptr);
        trans[i-1].gamma = str_to_double(token);

        // Wavelength
        token = strtok_r(NULL, DELIM, &saveptr);
        trans[i-1].lambda = str_to_double(token);

        // Total angular momentum of the ground state
        token = strtok_r(NULL, DELIM, &saveptr);
        trans[i-1].J_gnd = str_to_int(token);

        // Total angular momentum of the excited state
        token = strtok_r(NULL, DELIM, &saveptr);
        trans[i-1].J_exc = str_to_int(token);

        // Lande factor of the ground state
        token = strtok_r(NULL, DELIM, &saveptr);
        trans[i-1].g_gnd = str_to_double(token);

        // Lande factor of the excited state
        token = strtok_r(NULL, DELIM, &saveptr);
        trans[i-1].g_exc = str_to_double(token);
    }

    // Transitions set
    trans_set->num = i - 1;
    trans_set->trans = (transition_t*) calloc(trans_set->num, sizeof(transition_t));

    for(j = 0; j < trans_set->num; j++)
        trans_set->trans[j] = trans[j];

    // Free memory
    for(j = 0; j < i; j++)
        free(rows[j]);

    free(path);
    free(rows);
    free(trans);
}

void set_beams_setup(char *params_path, beams_setup_t *beams_setup){
    // Variables
    int i = 0, j = 0, n;
    char *path, *token, *saveptr, **rows;
    double delta = 0, s_0 = 0, w = 0;
    int trans_idx = -1;
    double *arr_sidebands;
    beams_modulation_t sidebands;
    beam_t *all_beams;

    // Allocate memory
    all_beams = (beam_t*) calloc(MAX_SIZE_ARRAY, sizeof(beam_t));

    // Read General Values
    //--
    // File path
    path = str_concatenate(params_path, "beams/general.csv");

    // Read lines
    rows = read_lines(path);

    // Check lines
    for(i = 1; !(rows[i] == NULL); i++){
        // Variable
        token = strtok_r(rows[i], DELIM, &saveptr);
        
        // Detuning
        //--
        if(strcmp(token, "delta") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            delta = str_to_double(token);
        }
        //--
        
        // Saturation parameter
        //--
        if(strcmp(token, "s_0") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            s_0 = str_to_double(token);
        }
        //--
        
        // Waist
        //--
        if(strcmp(token, "w") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            w = str_to_double(token);
        }
        //--
        
        // Transition ID
        //--
        if(strcmp(token, "trans_id") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            trans_idx = str_to_int(token);
        }
        //--
        
        // Sidebands
        //--
        if(strcmp(token, "sidebands") == 0){
            token = strtok_r(NULL, DELIM, &saveptr);
            arr_sidebands = get_double_array(token, &n);
            sidebands.num = (int) arr_sidebands[0];
            sidebands.freq = arr_sidebands[1];
        }
        //--
    }
    //--

    // Read all beams
    //--
    // File path
    path = str_concatenate(params_path, "beams/setup.csv");

    // Read lines
    rows = read_lines(path);

    // Check lines
    for(i = 1; !(rows[i] == NULL); i++){
        // Index
        token = strtok_r(rows[i], DELIM, &saveptr);
        all_beams[i-1].idx = str_to_int(token);

        // Wave vector
        token = strtok_r(NULL, DELIM, &saveptr);
        all_beams[i-1].k_dir = r3_normalize(get_double_array(token, &n));

        // Polarization
        token = strtok_r(NULL, DELIM, &saveptr);
        all_beams[i-1].pol_amp = r3_normalize(get_double_array(token, &n));

        // Saturation parameter
        token = strtok_r(NULL, DELIM, &saveptr);
        if(strcmp(token, "--") == 0) all_beams[i-1].s_0 = s_0;
        else all_beams[i-1].s_0 = str_to_double(token);

        // Waist
        token = strtok_r(NULL, DELIM, &saveptr);
        if(strcmp(token, "--") == 0) all_beams[i-1].w = w;
        else all_beams[i-1].w = str_to_double(token);

        // Detuning
        token = strtok_r(NULL, DELIM, &saveptr);
        if(strcmp(token, "--") == 0) all_beams[i-1].delta = delta;
        else all_beams[i-1].delta = str_to_double(token);

        // Sidebands
        token = strtok_r(NULL, DELIM, &saveptr);
        if(strcmp(token, "--") == 0) all_beams[i-1].sidebands = sidebands;
        else{
            arr_sidebands = get_double_array(token, &n);
            all_beams[i-1].sidebands.num = (int) arr_sidebands[0];
            all_beams[i-1].sidebands.freq = arr_sidebands[1];
        }

        // Transition
        token = strtok_r(NULL, DELIM, &saveptr);
        if(strcmp(token, "--") == 0) all_beams[i-1].trans_idx = trans_idx;
        else all_beams[i-1].trans_idx = str_to_double(token);
    }

    // Set beams setup
    beams_setup->num = i - 1;
    beams_setup->beams = (beam_t*) calloc(beams_setup->num, sizeof(beam_t));

    for(j = 0; j < beams_setup->num; j++)
        beams_setup->beams[j] = all_beams[j];
    //--

    // Free memory
    for(j = 0; j < i; j++)
        free(rows[j]);

    free(path);
    free(rows);
    free(all_beams);
}

void set_atom(char *params_path, atom_t *atom){
    // Variables
    int i = 0, j = 0;
    char *path, *token, *saveptr, **rows;

    // File path
    path = str_concatenate(params_path, "atom.csv");

    // Read lines
    rows = read_lines(path);

    // Check lines
    for(i = 1; !(rows[i] == NULL); i++){
        // Variable
        token = strtok_r(rows[i], DELIM, &saveptr);
        
        // Symbol
        //--
        if(strcmp(token, "symbol") == 0){
            token = strtok_r(NULL, DELIM, &saveptr); // Value
            atom->symbol = (char*) calloc(strlen(token), sizeof(char));
            strcpy(atom->symbol, token);
        }
        //--
        
        // Atomic Number
        //--
        if(strcmp(token, "Z") == 0){
            token = strtok_r(NULL, DELIM, &saveptr); // Value
            atom->Z = str_to_int(token);
        }
        //--
        
        // Mass
        //--
        if(strcmp(token, "mass") == 0){
            token = strtok_r(NULL, DELIM, &saveptr); // Value
            atom->mass = str_to_double(token);
        }
        //--
    }

    // Position and Velocity
    atom->pos = (double*) calloc(3, sizeof(double));
    atom->vel = (double*) calloc(3, sizeof(double));

    // Free memory
    for(j = 0; j < i; j++)
        free(rows[j]);

    free(path);
    free(rows);  
}

void set_ini_atom_state(atom_t *atom, initial_conditions_t ini_conds, magnetic_field_t B_params, beams_setup_t beams_setup, transition_t trans){
    // Variables
    int i;
    double beta, B_0, chi, delta, R, s_0, gamma;
    double z0 = 0;
    double **B_basis_t;

    // Verify theoretical position to optimize the simulation
    // ---
    if(are_all_beams_similars(beams_setup)){
        B_basis_t = r3_transposed_operator(B_params.B_basis);
        B_0 = B_params.B_0 *(0.5 * (B_basis_t[2][0] + B_basis_t[2][1]) - B_basis_t[2][2]);
        beta = (mu_B / h) * B_0 * 1e8;

        chi = - trans.g_gnd * trans.J_gnd;
        if(B_0 > 0) chi = trans.g_exc * (trans.J_gnd + 1);
        else chi = trans.g_exc * (trans.J_gnd - 1);

        gamma = 2*PI*trans.gamma*1e3;
        delta = beams_setup.beams[0].delta*gamma;
        s_0 = beams_setup.beams[0].s_0;

        R = (h * PI*trans.gamma) / (trans.lambda * atom->mass * u * g) * 1e5;
        z0 = (delta + gamma*0.5 * sqrt((R - 1)*s_0 - 1)) / (beta * chi);

        //printf("slope = %.15f um/kHz\n", (1 / (beta*chi))*1e9);
        //printf("z0 = %f m\n", z0);
        //exit(0);
    }
    // ---

    // Set position at origin
    for(i = 0; i < 3; i++)
        atom->pos[i] = B_params.B_origin[i]; // cm;

    if(z0 < 0) atom->pos[2] += z0;

    //r3_print(atom->pos, "pos [cm]");
    //exit(0);

    //atom->pos[0] = 0.5;
    //atom->pos[1] = 0.5;
    //atom->pos[2] = 0.5;

    // Set velocity according a Maxwell-Boltzmann distribution    
    for(i = 0; i < 3; i++)
        atom->vel[i] = random_norm(0, sqrt(1e-2 * k_B * ini_conds.T_0 / (atom->mass * u))*100); // cm / s

    //atom->vel[0] = 5.0;
    //atom->vel[1] = 5.0;
    //atom->vel[2] = 5.0;
}

void set_hist(results_t *res, performance_t perform){
    // Variables
    int i;

    // Allocate memory
    res->pos_hist = (histogram_t*) calloc(3, sizeof(histogram_t));
    res->vel_hist = (histogram_t*) calloc(3, sizeof(histogram_t));

    for(i = 0; i < 3; i++){
        // Position
        res->pos_hist[i].num_bins = perform.num_bins;
        res->pos_hist[i].bin_size = 2 * perform.max_r / res->pos_hist[i].num_bins;
        res->pos_hist[i].coord0 = - perform.max_r;
        res->pos_hist[i].freqs = (long int*) calloc(res->pos_hist[i].num_bins, sizeof(long int));

        // Velocity
        res->vel_hist[i].num_bins = perform.num_bins;
        res->vel_hist[i].bin_size = 2 * perform.max_v / res->vel_hist[i].num_bins;
        res->vel_hist[i].coord0 = - perform.max_v;
        res->vel_hist[i].freqs = (long int*) calloc(res->vel_hist[i].num_bins, sizeof(long int));
    }
}

double *magnetic_field(magnetic_field_t B_params, double *r){
    // Variables
    int i;
    double *B;          // Magnetic field vector

    // Allocate memory for the magnetic field vector
    B = (double*) calloc(3, sizeof(double)); // G

    // Change basis of the position vector
    B = r3_apply_operator(B_params.B_basis, r); // mm

    // Get magnetic field from B_basis
    B[0] = B[0] / 2;
    B[1] = B[1] / 2;
    B[2] = -B[2];
    B = r3_scalar_product(B_params.B_0, B);

    // Get back from the lab basis
    B = r3_apply_operator(r3_transposed_operator(B_params.B_basis), B);

    // Linear magnetic field gradient
    for(i = 0; i < 3; i++)
        B[i] += B_params.B_lin_grad[i]*r[i]; // G

    // Bias magnetic field
    B = r3_sum(B_params.B_bias, B); // G

    return B;
}

double move(atom_t *atom, transitions_set_t trans_set, scattering_t scatt, magnetic_field_t B_params, int g_bool){
    // Variables
    int i;
    double dt;
    double *a;
    double photon_speed;            // Speed gained due to photon absorption or emission

    // Allocate memory
    a = (double*) calloc(3, sizeof(double));

    // Time interval
    dt = scatt.dt*1e-3; // s

    // Kinematics
    //--
    // Accelerations
    //--
    // Gravity
    if(g_bool) a[2] = - g * 1e2; // cm / s^2

    // Magnetic acceleration
    for(i = 0; i < trans_set.num; i++)
        a = r3_sum(a, magnetic_acceleration(*atom, trans_set.trans[i], B_params)); // cm / s^2
    //--

    // Update position [cm]
    atom->pos = r3_sum(atom->pos, r3_scalar_product(dt, atom->vel)); // cm
    atom->pos = r3_sum(atom->pos, r3_scalar_product((dt*dt / 2), a)); // cm

    // Update velocity [cm/s]
    //--
    // Acceleration effect
    atom->vel = r3_sum(atom->vel, r3_scalar_product(dt, a)); // cm / s

    // Scattering event
    if(scatt.beam.idx >= 0){
        photon_speed = 1e4 * h  / (atom->mass * u * scatt.trans.lambda); // cm / s

        // Absorption
        atom->vel = r3_sum(atom->vel, r3_scalar_product(photon_speed, r3_normalize(scatt.beam.k_dir))); // cm / s

        // Emission
        atom->vel = r3_sum(atom->vel, r3_scalar_product(photon_speed, random_r3_vector())); // cm / s
    }
    //--
    //--

    // Release memory
    free(a);

    return dt*1e3; // ms
}

scattering_t scattering_event(atom_t atom, transitions_set_t trans_set, beams_setup_t beams_setup, magnetic_field_t B_params, double max_dt){
    // Variables
    int i;                      // Auxiliaries
    scattering_t scatt;         // Return
    double *R, R_max;           // Scattering rates
    double *probs;              // Probabilities to absorb each beam
    int chosen_beam = 0;        // Chosen beam

    // Allocate memory
    probs = (double*) calloc(beams_setup.num + 1, sizeof(double));

    // Get scattering rates
    R = get_scatt_rates(beams_setup, atom, trans_set, B_params); // Hz
    R_max = max_arr(R, beams_setup.num); // Hz

    // Set time interval [ms]
    scatt.dt = 1 / ((beams_setup.num)*R_max);
    scatt.dt *= 1e3;

    // Check maximum value
    if(scatt.dt > max_dt) scatt.dt = max_dt;

    // Get probabilities
    //--
    for(i = 0; i < beams_setup.num; i++){
        probs[i + 1] = R[i] * scatt.dt * 1e-3;
        probs[0] += probs[i + 1];

        //r3_print(beams_setup.beams[i].k_dir, "k");
        //printf("R[%d] = %f\n\n", i + 1, R[i]);
    }

    // Probability of the beam not being absorbed
    probs[0] = 1 - probs[0];
    //--

    // Pick or not a beam
    chosen_beam = random_pick(probs, beams_setup.num + 1);
    if(chosen_beam > 0) {
        scatt.beam = beams_setup.beams[chosen_beam-1];
        scatt.trans = get_transition_by_idx(scatt.beam.trans_idx, trans_set);
    } else scatt.beam.idx = -1;

    // Free memory
    free(R);
    free(probs);

    // Return
    return scatt;
}

double *get_scatt_rates(beams_setup_t beams_setup, atom_t atom, transitions_set_t trans_set, magnetic_field_t B_params){
    // Variables
    int i, j, k, m;                                         // Auxiliaries
    double *R;                                              // (Return) All scattering rates
    double *B, *eB;                                         // Magnetic field
    double **C;                                             // Orthonormal basis
    double r;                                               // Distance
    double s, s_0, gamma, lambda;                           // Parameters
    double doppler_shift, zee_shift, delta;                 // Detunings
    double laser_detuning, laser_delta;                     // Laser detunings
    transition_t trans;                                     // Transition
    beam_t beam;                                            // Beam
    int pol[] = {+1, -1, 0};                                // Polarizations

    // Allocate memory
    R = (double*) calloc(beams_setup.num, sizeof(double));
    eB = (double*) calloc(3, sizeof(double));

    // Magnetic field
    B = magnetic_field(B_params, atom.pos);
    if(r3_mod(B) == 0) 
        for(i = 0; i < 3;i++) 
            eB[i] = B_params.B_basis[2][i];
    else eB = r3_normalize(B);

    //r3_print(B, "B");
    //r3_print(eB, "eB");
    //exit(0);

    // Checking all beams out
    for(i = 0; i < beams_setup.num; i++){
        R[i] = 0;

        // Set beam
        beam = copy_beam(beams_setup.beams[i]);
        //r3_print(beam.k_dir, "k");

        // Set polarization amplitudes
        set_polarizations_amplitudes(&beam, B, eB);
        //r3_print(beam.pol_amp, "pol");

        // Get transition
        trans = get_transition_by_idx(beam.trans_idx, trans_set);

        // Transition parameters
        gamma = trans.gamma; // 2pi kHz
        laser_detuning = beam.delta * gamma; // 2pi kHz
        lambda = trans.lambda; // nm

        //printf("gamma [2pi kHz] = %f\n", gamma);
        //printf("delta [2pi kHz] = %f\n", laser_detuning);
        //printf("lambda [nm] = %f\n", lambda);
        //printf("\n");
        //exit(0);

        // Effective saturation parameter
        s_0 = beam.s_0 / (2*beam.sidebands.num + 1);
        atom.pos[0] = 0;
        atom.pos[1] = 0;
        atom.pos[2] = -1;
        C = orthonormal_basis(r3_normalize(beam.k_dir));
        r = pow(r3_inner_product(C[0], atom.pos), 2);
        r += pow(r3_inner_product(C[1], atom.pos), 2);
        r = sqrt(r);
        //printf("r [cm] = %f\n", r);
        s_0 = s_0 * exp(-2 * pow((r / (beam.w)), 2));
        //printf("s_0 = %f\n", s_0);
        //exit(0);

        // Doppler shift
        doppler_shift = - 1e4 * (r3_inner_product(C[2], atom.vel)) / lambda; // 2pi kHz
        //printf("doppler_shift [2pi kHz] = %f\n", doppler_shift);
        //r3_print(atom.vel, "vel [cm/s]");
        //printf("doppler_shift [Gamma] = %f\n", doppler_shift / gamma);
        //exit(0);

        // Check each polarization
        for(j = 0; j < 3; j++){
            // Saturation parameter
            s = s_0 * beam.pol_amp[j];

            if(s >= 0.0){
                //printf("s = %f\n", s);

                // Zeeman shift
                zee_shift = zeeman_shift(atom.pos, trans, B_params, pol[j]); // 2pi kHz
                //printf("zeeman_shift(%d) = %f\n", pol[j], zee_shift);

                // Scattering rate
                delta = laser_detuning + doppler_shift + zee_shift; // 2pi kHz
                //printf("delta [2pi kHz] = %f\n", delta);
                R[i] += PI * gamma * 1e3 * s / (1 + s + 4*pow(delta/gamma, 2)); // Hz

                // Sidebands
                for(m = 0; m < beam.sidebands.num; m++){
                    // Right sideband
                    laser_delta = laser_detuning + (m+1) * beam.sidebands.freq; // 2pi kHz
                    delta = laser_delta + zee_shift + doppler_shift;
                    R[i] += (PI*gamma * 1e3) * s / (1 + s + 4*pow(delta/gamma, 2)); // Hz

                    // Left sideband
                    laser_delta = laser_detuning - (m+1) * beam.sidebands.freq; // 2pi kHz
                    delta = laser_delta + zee_shift + doppler_shift;
                    R[i] += (PI*gamma * 1e3) * s / (1 + s + 4*pow(delta/gamma, 2)); // Hz
                }
            }
        }
        //printf("R[%d] = %f\n", i+1, R[i]);
        //printf("\n");

        // Release memory
        for(k = 0; k < 3; k++)
            free(C[k]);

        free(C);
    }

    //exit(0);

    // Release memory
    free(B);

    return R;
}

transition_t get_transition_by_idx(int idx, transitions_set_t trans_set){
    // Variables
    int i;
    transition_t trans;
    trans.idx = -1;

    for(i = 0; i < trans_set.num; i++){
        if(trans_set.trans[i].idx == idx){
            trans = trans_set.trans[i];
            break;
        }
    }

    return trans;
}

beam_t get_beam_by_idx(int idx, beams_setup_t beams_setup){
    // Variables
    int i;
    beam_t beam;
    beam.idx = -1;

    for(i = 0; i < beams_setup.num; i++){
        if(beams_setup.beams[i].idx == idx){
            beam = beams_setup.beams[i];
            break;
        }
    }

    return beam;
}

void set_polarizations_amplitudes(beam_t *beam, double *B, double *eB){
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


    // Release memory
    for(i = 0; i < 3; i++){
        free(M_P[i]);
        free(M_P_dagger[i]);
        free(M_L[i]);
        free(M_B[i]);
        free(R3_basis[i]);
    }

    free(R3_basis);
    free(M_P);
    free(M_P_dagger);
    free(M_L);
    free(M_B);
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

double zeeman_shift(double *pos, transition_t trans, magnetic_field_t B_params, int pol){
    // Variables
    double chi, g_gnd, g_exc, mj_gnd, mj_exc;       // Transition
    double *B, zee_shift;                           // Zeeman shift

    // Transition
    g_gnd = trans.g_gnd;
    g_exc = trans.g_exc;
    mj_gnd = -trans.J_gnd;
    mj_exc = mj_gnd + pol;
    chi = g_gnd * mj_gnd - g_exc * mj_exc;

    // Magnetic field [Gauss]
    B = magnetic_field(B_params, pos);

    // Zeeman shift (units of Gamma)
    zee_shift = 1e3 * (mu_B / h) * r3_mod(B) * chi; // 2pi kHz

    // Release memory
    free(B);

    return zee_shift;
}

double *magnetic_acceleration(atom_t atom, transition_t trans, magnetic_field_t B_params){
    // Variables
    double *a_B;        // Magnetic acceleration
    double *r;          // Position vector
    double *del_B;
    double norm;

    // Allocate memory
    del_B = (double*) calloc(3, sizeof(double));
    a_B = (double*) calloc(3, sizeof(double));

    // Get position considering the axis rotating
    r = r3_apply_operator(B_params.B_basis, atom.pos);
    norm = sqrt(r[0]*r[0] + r[1]*r[1] + 4*r[2]*r[2]);

    if(norm > 0){
        del_B[0] = r[0] / 2;
        del_B[1] = r[1] / 2;
        del_B[2] = 2 * r[2];

        del_B = r3_scalar_product(B_params.B_0 / norm, del_B);
        del_B =  r3_apply_operator(r3_transposed_operator(B_params.B_basis), del_B);

        a_B = r3_scalar_product(0.1 * trans.g_gnd * trans.J_gnd * mu_B / (atom.mass * u), del_B);
    }

    return a_B;
}

beam_t copy_beam(beam_t beam){
    // Variables
    int i;
    beam_t copied_beam;

    // Index
    copied_beam.idx = beam.idx;

    // Wave vector
    copied_beam.k_dir = (double*) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++)
        copied_beam.k_dir[i] = beam.k_dir[i];

    // Polarization
    copied_beam.pol_amp = (double*) calloc(3, sizeof(double));
    for(i = 0; i < 3; i++)
        copied_beam.pol_amp[i] = beam.pol_amp[i];

    // Detuning
    copied_beam.delta = beam.delta;

    // Saturation parameter
    copied_beam.s_0 = beam.s_0;

    // Waist
    copied_beam.w = beam.w;

    // Transition index
    copied_beam.trans_idx = beam.trans_idx;

    // Modulation
    copied_beam.sidebands = beam.sidebands;

    return copied_beam;
}

int is_inside_threshold(atom_t atom, transitions_set_t trans_set, magnetic_field_t B_params, beams_setup_t beams_setup){
    // Variables
    int i;
    int check = 1;
    double rho;

    // Check all beams
    for(i = 0; i < beams_setup.num; i++){
        if(beams_setup.beams[i].s_0 > 0){
            rho = sqrt(r3_inner_product(atom.pos, atom.pos) - pow(r3_inner_product(atom.pos, beams_setup.beams[i].k_dir), 2));
            if(rho > beams_setup.beams[i].w){
                check = 0;
                break;
            }
        }
    }

    return check;
}

int are_all_beams_similars(beams_setup_t beams_setup){
    // Variables
    int i, bool_ret = 1;
    double delta, s_0;

    delta = beams_setup.beams[0].delta;
    s_0 = beams_setup.beams[0].s_0;

    if(beams_setup.num > 1){
        for(i = 1; i < beams_setup.num; i++){
            if((beams_setup.beams[i].delta) > delta*0.99 || (beams_setup.beams[i].delta) < delta*1.01){
                bool_ret = 0;
                break;
            }

            else if((beams_setup.beams[i].s_0) > s_0*1.01 || (beams_setup.beams[i].s_0) < s_0*0.99){
                bool_ret = 0;
                break;
            }
        }
    } else bool_ret = 1;

    return bool_ret;
}
//--

// Utility
//--
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
    char *token, *saveptr;
    double aux_arr[max_size];
    double *arr;

    str = str + 1;
    str[strlen(str)-1] = '\0';
    token = strtok_r(str, " ", &saveptr); 

    //
    // Parse string
    //

    i = 0;
    while(token && (i < max_size)){
        // Py MODULE
        aux_arr[i] = str_to_double(token);
        token = strtok_r(NULL, " ", &saveptr);
        i++;
    }

    arr = (double *) malloc(i * sizeof(double));
    for(j = 0; j < i; j++) arr[j] = aux_arr[j];

    // Release memory
    free(token);

    return arr;
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

double str_to_double(char *str){
    // Variables
    double num;

    // Python module
    if(Py_MODULE) {
        if(str[0] == '-') num = -atof(str_replace(str + 1, ".", ",")); 
        else num = atof(str_replace(str, ".", ","));
    }

    // C program
    else {
        if(str[0] == '-') num = -atof(str + 1); 
        else num = atof(str);
    }

    return num;
}

int str_to_int(char *str){
    // Variables
    int num;

    // Python module
    if(Py_MODULE) {
        if(str[0] == '-') num = -atoi(str_replace(str + 1, ".", ",")); 
        else num = atoi(str_replace(str, ".", ","));
    }

    // C program
    else { 
        if(str[0] == '-') num = -atoi(str + 1); 
        else num = atoi(str);
    }

    return num;
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

double **orthonormal_basis(double *v){
    //
    // Variables
    //

    int i, v3_dir = -1;
    double **B;             // Desired basis
    double *v1, *v2, *v3;   // Auxiliary vectors

    // Allocate memory
    v1 = (double*) calloc(3, sizeof(double));

    // Normalize vector v
    v3 = r3_normalize(v);

    // Verify if the vector v3 is in x, y, or z directions
    for(i = 0;  i < 3; i++){
        if(v3[i] == 1){
            v3_dir = i;
            break;
        }
    }

    // If v3 is in x,y, or z directions, chose x, y, or z for the v1 direction
    if(v3_dir >= 0){
        if(v3_dir == 2 || v3_dir == 1) v1[0] = 1;
        else v1[2] = 1; 
    } else {
        // Generate a random vector  
        for(i = 0; i < 3; i++) v1[i] = ((double) rand()) / ((double) RAND_MAX);

        // Define a orthonormal vector
        v2 = r3_scalar_product(r3_inner_product(v1, v3), v3);
        v1 = r3_normalize(r3_diff(v1, v2));
    }

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

double max_arr(double *arr, int size){
    // Variables
    int i;
    double max_val = arr[0];

    for(i = 1; i < size; i++){
        if(arr[i] > max_val) max_val = arr[i];
    }

    return max_val;
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

double random_exp(double mean){
    return (- mean * log(1 - ((double) rand()) / ((double) RAND_MAX)));
}

void update_hist(histogram_t *hist, double val){
    // Variables
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
}
//--

// Debug
//--
void print_initial_conditions(initial_conditions_t ini_conds){
    printf("Initial Conditions\n--\n");
    printf("T_0 [uK] = %f\n", ini_conds.T_0);
    printf("v_0 [m/s] = %f\n", ini_conds.v_0);
    printf("g_bool = %d\n", ini_conds.g_bool);
    printf("\n");
}

void print_performance(performance_t perform){
    printf("Performance\n--\n");
    printf("recording_time [ms] = %f\n", perform.recording_time);
    printf("max_r [cm] = %f\n", perform.max_r);
    printf("max_v [cm/s] = %f\n", perform.max_v);
    printf("num_bins = %d\n", perform.num_bins);
    printf("wait_time [ms] = %f\n", perform.wait_time);
    printf("time interval [ms] = %f\n", perform.dt);
    printf("\n");
}

void print_magnetic_field(magnetic_field_t B_params){
    printf("Magnetic field\n--\n");
    printf("B_0 [G/cm] = %f\n", B_params.B_0);
    r3_operator_print(B_params.B_basis, "B_basis");
    r3_print(B_params.B_bias, "B_bias [G]");
    r3_print(B_params.B_lin_grad, "B_lin_grad [G/cm]");
    r3_print(B_params.B_origin, "B_origin [G]");
    printf("\n");
}

void print_transitions(transitions_set_t trans_set){
    // Variables
    int i = 0;

    printf("All transitions (n = %d)\n--\n", trans_set.num);
    for(i = 0; i < trans_set.num; i++){
        printf("Index = %d\n", trans_set.trans[i].idx);
        printf("gamma [2pi kHz] = %f\n", trans_set.trans[i].gamma);
        printf("lambda [nm] = %f\n", trans_set.trans[i].lambda);
        printf("J_gnd = %d\n", trans_set.trans[i].J_gnd);
        printf("J_exc = %d\n", trans_set.trans[i].J_exc);
        printf("g_gnd = %f\n", trans_set.trans[i].g_gnd);
        printf("g_exc = %f\n", trans_set.trans[i].g_exc);
        printf("--\n");
    }
    printf("\n");
}

void print_beams_setup(beams_setup_t beams_setup){
    // Variables
    int i = 0;

    printf("All Beams (n = %d)\n--\n", beams_setup.num);

    for(i = 0; i < beams_setup.num; i++){
        printf("Index = %d\n", beams_setup.beams[i].idx);
        r3_print(beams_setup.beams[i].k_dir, "k");
        r3_print(beams_setup.beams[i].pol_amp, "pol");
        printf("s_0 = %f\n", beams_setup.beams[i].s_0);
        printf("w [cm] = %f\n", beams_setup.beams[i].w);
        printf("delta [Gamma] = %f\n", beams_setup.beams[i].delta);
        printf("num_sidebands = %d\n", beams_setup.beams[i].sidebands.num);
        printf("freq_sidebands [kHz] = %f\n", beams_setup.beams[i].sidebands.freq);
        printf("trans_idx = %d\n", beams_setup.beams[i].trans_idx);
        printf("--\n");
    }
    printf("\n");
}

void print_atom(atom_t atom){
    printf("Atom\n--\n");
    printf("symbol = %s\n", atom.symbol);
    printf("Z = %d\n", atom.Z);
    printf("mass [u] = %f\n", atom.mass);
    printf("--\n");
}

void print_progress(atom_t atom, results_t res, magnetic_field_t B_params, int progress){
    // Variables
    double *r;

    r = r3_diff(atom.pos, B_params.B_origin);

    printf("progress = %d\n", progress);
    printf("--\n");
    printf("t [ms] = %f\n", res.time);
    r3_print(r, "pos (from B_origin) [cm]");
    r3_print(atom.vel, "vel [cm / s]");
    printf("|pos| [cm] = %f\n", r3_mod(r));
    printf("|vel| [cm / s] = %f\n", r3_mod(atom.vel));
    printf("--\n");
}

int print_results(results_t res, atom_t atom, int show_hist){
    // Variables
    int i, j;

    printf("Simulation status\n--\n");
    r3_print(res.ini_vel, "Initial Velocity [cm/s]");
    printf("Total time [ms] = %f\n", res.time);
    printf("Trapped atom = %d\n", res.trapped_atom);
    printf("\n");

    if(show_hist){
        //
        // Position frequencies
        //-
        for(i = 0; i < 3;i++){
            printf("pos[%d] = [\n", i+1);
            for(j = 0; j < res.pos_hist[i].num_bins; j++)
                printf("%ld ", res.pos_hist[i].freqs[j]);
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
                printf("%ld ", res.vel_hist[i].freqs[j]);
            printf("]\n\n");
        }
        printf("\n");
        //--
    }

    return 1;
}
//--
