//
// Libraries
//

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "vectors.h"

// 
// Constants
//

#define STRING_BUFFER_SIZE 1024
#define MAX_LINES 124
#define DELIM ","
#define MAX_BEAMS 16
#define Py_MODULE 0

#define h 6.62607004        // Planck constant [10^{-34} J s]
#define e 1.60217662        // Elementary charge [10^{-19} C]s
#define c 2.99792458        // Speed of light [10^{8} m / s]
#define k_B 1.38064852      // Boltzmann constant [10^{-23} J / K]
#define mu_B 9.274009994    // Bohr magneton [10^{-24} J / T]
#define u 1.660539040       // Atomic mass unit [10^{-27} kg]
#define g 980.665           // Gravitational acceleration [cm / s^2]
#define PI 3.14159265358

// Structures
//--
// Transition
typedef struct{
    int idx;        /* Index */
    double gamma;   /* Transition rate [kHz] */
    double lambda;  /* Resonant wave length */
    int J_gnd;      /* Total angular momentum of the ground state */
    int J_exc;      /* Total angular momentum of the excited state */
    double g_gnd;   /* Landè factor of the ground state */
    double g_exc;   /* Landè factor of the excited state */
} transition_t;

// Transitions Set
typedef struct{
    int num;                /* Number of Transitions */
    transition_t *trans;    /* Array of transitions */
} transitions_set_t;

// Atom
typedef struct{
    char *symbol;               /* Atom symbol */
    int Z;                      /* Atomic number */
    double mass;                /* Atomic mass */
    double *pos;                /* Position [cm] */
    double *vel;                /* Velocity [cm / s] */
} atom_t;

// Performance
typedef struct{
    double T_0;                 /* Initial temperature [uK] */
    double recording_time;      /* Recording time [tau] */
    double max_r;               /* Maximum distance (threshold) [cm] */
    double max_v;               /* Maximum speed [cm/s] */
    int num_bins;               /* Number of bins in each histogram */
    double wait_time;           /* Time to reach the equilibrium [tau] */
    double dt;                  /* Time interval [tau] */
} performance_t;

// Initial conditions
typedef struct{
    double T_0;         /* Initial temperature [uK] */
    double v_0;         /* Module of the initial velocity of the atoms */
    double *v_0_dir;    /* Direction of the initial velocity of the atoms */
    double *r_0_dir;    /* Direction of the initial position of the atoms */
    double R;           /* Radius of the atomic beam */
    int g_bool;         /* Gravity (1 - Considered, 0 - Do not consider) */
} initial_conditions_t;

// Magnetic field
typedef struct{
    double B_0;             /* Magnetic Field gradient */
    double **B_basis;       /* Reference basis of the magnetic field */
    double *B_bias;         /* Bias of magnetic field */
    double *B_lin_grad;     /* Linear gradient of magnetic field */
    double *B_origin;       /* Origin */
} magnetic_field_t;

// Beams Modulation
typedef struct {
    int num;        /* Number of sidebands */
    double freq;    /* Resonant frequency */
} beams_modulation_t;

// Beam
typedef struct {
    int idx;                    /* Index */
    double *k_dir;              /* Wave vector direction */
    double *pol_amp;            /* Polarization amplitudes */
    double delta;               /* Laser detuning */
    double s_0;                 /* Resonant saturation parameter */
    double w;                   /* Waist radius */
    beams_modulation_t mod;     /* Beams Modulation */
} beam_t;

// Beams setup
typedef struct {
    int num;                /* Number of beams */
    beam_t *beams;          /* All beams */
} beams_setup_t;

// Histogram
typedef struct {
    double num_bins;    /* Number of bins */
    double bin_size;    /* Bin size */
    double coord0;      /* Initial value */
    int *freqs;         /* Frequencies */
} histogram_t;

// Results
typedef struct{
    histogram_t *pos_hist;          /* Marginal positions histograms */
    histogram_t *vel_hist;          /* Marginal velocities histograms */
    double *ini_vel;                /* Inicial velocity */
    double time;                    /* Total time [s] */
    double escape_square_vel;       /* Escape square velocity [cm/s] */
    int trapped_atom;               /* 1 - Atom was trapped, 0 - Atom was not trapped */
} results_t;
//--

/*
    Flags:
    0 - Simulations has stopped due an error
    1 - Simulations has finished
    2 - Simulation is running
*/

//
// Main functions
//

// Run simulation for a single atom
results_t simulate_atom(char *params_path, int marginals, long seed_time);

// Get parameters of performance
performance_t get_performance(char *params_path);

// Get initial conditions
initial_conditions_t get_initial_conditions(char *params_path);

// Get parameters of the magnetic field
magnetic_field_t get_magnetic_field(char *params_path, performance_t perform);

// Get beams setup
beams_setup_t get_beams(char *params_path);

// Get transition
transition_t get_transition(char *params_path);

// Get atom
atom_t get_atom(initial_conditions_t ini_conds, performance_t perform, beams_setup_t beams_setup, magnetic_field_t B_params, zeeman_slower_t ZS, int opt, char *params_path);

// Get Zeeman slower parameters
zeeman_slower_t get_zeeman_slower(char *params_path);

// Get compression parameters
compression_t get_compression(char *params_path);

// Set initial atom state
int set_ini_atom_state(atom_t *atom, beams_setup_t beams_setup, magnetic_field_t B_params, performance_t perform, initial_conditions_t ini_conds,  zeeman_slower_t ZS, int opt);

// Set histograms
int set_hist(int only_marginals, results_t *res, performance_t conds);

// Apply movement on the atom due to photonic recoil, gravitational force and magnetic force
double move(atom_t *atom, beams_setup_t beams_setup, performance_t conds, initial_conditions_t ini_conds, magnetic_field_t B_params, zeeman_slower_t ZS);

// Get magnetic field vector in the lab frame
double *magnetic_field(magnetic_field_t B_params, double *r);

// Get polarizations amplitudes
//int set_polarizations_amplitudes_old(beam_t *beam, double *eB);

// Get polarizations amplitudes
int set_polarizations_amplitudes(beam_t *beam, double *eB);

// Get scattering rates of each beam
double *get_scatt_rates(beams_setup_t beams_setup, magnetic_field_t B_params, atom_t atom);

// Get scattering rates of the Zeeman slower beam
double get_ZS_scatt_rates(zeeman_slower_t ZS, magnetic_field_t B_params, atom_t atom);

// Get magnetic acceleration
double *magnetic_acceleration(atom_t atom, magnetic_field_t B_params);

// Get Zeeman shift (units of Gamma)
double zeeman_shift(atom_t atom, magnetic_field_t B_params, int pol);

// Check if the atom passes through the threshold
int is_inside_threshold(atom_t atom, magnetic_field_t B_params, beams_setup_t beams_setup, performance_t perform);

// Change-of-matrix basis from polarization frame P to cartisian frame C
complex_t **change_of_basis_from_P_to_C();

/*

// Get the probability to absorb each beam
double *get_probs(beams_setup_t beams_setup, magnetic_field_t B_params, atom_t atom, double dt);

// Get a list of scattering rate for each beam considering all transitions and sidebands
double *get_all_scatt_rate(beams_setup_t beams_setup, magnetic_field_t B_params, atom_t atom);

// Get the sum of scattering rates over all possibilities
double get_total_scatt_rate(beam_t beam, sidebands_t sidebands, double *B, atom_t atom);

// Get all scattering rates for a beam
double *get_scatt_rates(beam_t beam, sidebands_t sidebands, double *B, atom_t atom);

// Get environment
environment_t get_environment(char *params_path);

// Photon absorption event
double photon_absorption(atom_t *atom, beams_setup_t beams_setup, performance_t conds, environment_t env, double dt);

// Photon emission event
double photon_emission(atom_t *atom, environment_t env);

// Get scattering rate
double scattering_rate(atom_t atom, polarized_beam_t beam, performance_t conds, environment_t env, double *B);

// Get components of a vector v on the basis B given the components on basis A
double *change_basis(double *v, double **A, double **B_basis);

*/

//
// Utility functions
//

// Concatenate strings
char *str_concatenate(char *str1, char *str2);

// Read lines from a file
char **read_lines(char *path);

// Replace a character in a string
char* str_replace(char *orig, char *rep, char *with);

// Get double array from string in the format [f1 f2 ... fn]
double *get_double_array(char *str, int *size);

// Generate a orthonormal basis given a vector
double **orthonormal_basis(double *v3);

// Generate a double random number following a Gaussian distribution given a mean and a standard deviation
double random_norm(double mean, double std_dev);

// Generate a double random number following a Exponential distribution given a mean
double random_exp(double mean);

// Generate a double random number following a uniform distribution given by the interval (a, b)
double random_uniform(double a, double b);

// Pick randomly a element of an integer array given an array of probabilities
int random_pick(double *probs, int size);

// Get normalized random vector
double *random_r3_vector(void);

// Update histogram
int update_hist(histogram_t *hist, double val);

// Update multidimensional histogram
int update_hist_3d(histogram_3d_t *hist, double *vals);

// Rotating matrix about
// axis [1 (x), 2 (y), 3 (z)]
// theta [Degree]
double **rotating_matrix(double theta, int axis);

//
// Debug

// Print parameters of performance
int print_performance(performance_t perform);

// Print parameters of initial conditions
int print_initial_conditions(initial_conditions_t ini_conds);

// Print parameters of the magnetic field
int print_magnetic_field(magnetic_field_t B_params);

// Print parameters of the atom
int print_atom(atom_t atom);

// Print parameters of the beams
int print_beams(beams_setup_t beams_setup);

// Print parameters of Zeeman slower
int print_zeeman_slower(zeeman_slower_t ZS);

// Print results
int print_results(results_t res, atom_t atom, int opt, int show_hist);

// Print results
int print_compression(compression_t compression);

// Print simulation status
int print_status(atom_t atom, results_t res, magnetic_field_t B_params, beams_setup_t beams_setup, double progress);

/**

// (Debug) print parameters
int print_params(atom_t atom, performance_t conds, beams_setup_t beams, environment_t env);


// Get int array from string in the format [i1 i2 ... in]
int *get_int_array(char *str, int *size);

// Convert the components (module) of a polarization vector on the basis C to a basis D
double *update_polarization_vector(double *eps, double **r3_C, double **r3_D);

**/