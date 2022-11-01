// Libraries
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "vectors.h"

// Constants
//--
#define STRING_BUFFER_SIZE 1024
#define MAX_LINES 124
#define DELIM ","
#define MAX_SIZE_ARRAY 16
#define Py_MODULE 0

#define h 6.62607004        // Planck constant [10^{-34} J s]
#define e 1.60217662        // Elementary charge [10^{-19} C]s
#define c 2.99792458        // Speed of light [10^{8} m / s]
#define k_B 1.38064852      // Boltzmann constant [10^{-23} J / K]
#define mu_B 9.274009994    // Bohr magneton [10^{-24} J / T]
#define u 1.660539040       // Atomic mass unit [10^{-27} kg]
#define g 9.80665           // Gravitational acceleration [m / s^2]
#define PI 3.14159265358
//--

// Structures
//--
// Transition
typedef struct{
    int idx;        /* Index */
    double gamma;   /* Transition rate [2pi kHz] */
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
    double mass;                /* Atomic mass [Da] */
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
    int g_bool;         /* Gravity (1 - Considered, 0 - Do not consider) */
} initial_conditions_t;

// Magnetic field
typedef struct{
    double B_0;             /* Magnetic Field gradient [G / mm] */
    double **B_basis;       /* Reference basis of the magnetic field */
    double *B_bias;         /* Bias of magnetic field [G] */
    double *B_lin_grad;     /* Linear gradient of magnetic field [G / mm] */
    double *B_origin;       /* Origin [mm] */
} magnetic_field_t;

// Beams Modulation
typedef struct {
    int num;        /* Number of sidebands */
    double freq;    /* Resonant frequency */
} beams_modulation_t;

// Beam
typedef struct {
    int idx;                        /* Positive index */
    double *k_dir;                  /* Wave vector direction */
    double *pol_amp;                /* Polarization amplitudes */
    double delta;                   /* Laser detuning */
    double s_0;                     /* Resonant saturation parameter */
    double w;                       /* Waist radius */
    int trans_idx;                  /* Transition index */
    beams_modulation_t sidebands;   /* Beams Modulation */
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
    long int *freqs;    /* Frequencies */
} histogram_t;

// Results
typedef struct{
    histogram_t *pos_hist;          /* Marginal positions histograms */
    histogram_t *vel_hist;          /* Marginal velocities histograms */
    double *ini_vel;                /* Inicial velocity */
    double time;                    /* Total time [s] */
    double escape_vel;              /* Escape velocity [cm/s] */
    int trapped_atom;               /* 1 - Atom was trapped, 0 - Atom was not trapped */
} results_t;

// Scattering
typedef struct{
    beam_t beam;        /* Beam which will be absorbed */
    transition_t trans; /* Atomic transitions associated with the beam wavelength */
    double dt;          /* [s] Time to happen the scattering event */
} scattering_t;
//--

// Main
//--
// Time evolution of an atom
results_t simulate_atom(char *params_path, int opt, long seed_time);

// Set initial conditions
void set_initial_conditions(char *params_path, initial_conditions_t *ini_conds);

// Set performance parameters
void set_performance(char *params_path, performance_t *perform);

// Set magnetic field parameters
void set_magnetic_field(char *params_path, magnetic_field_t *B_params);

// Set transitions
void set_transitions(char *params_path, transitions_set_t *trans_set);

// Set beams setup
void set_beams_setup(char *params_path, beams_setup_t *beams_setup);

// Set atom
void set_atom(char *params_path, atom_t *atom);

// Set initial atom state
void set_ini_atom_state(atom_t *atom, initial_conditions_t ini_conds, magnetic_field_t B_params, beams_setup_t beams_setup, transition_t trans);

// Set histograms
void set_hist(results_t *res, performance_t perform);

// Get magnetic field vector on the lab frame
double *magnetic_field(magnetic_field_t B_params, double *r);

// Apply movement on the atom due to photonic recoil, gravitational force and magnetic force
double move(atom_t *atom, transitions_set_t trans_set, scattering_t scatt, magnetic_field_t B_params, int g_bool);

// Execute a scattering event
scattering_t scattering_event(atom_t atom, transitions_set_t trans_set, beams_setup_t beams_setup, magnetic_field_t B_params, double max_dt);

// Get scattering rates of each beam
double *get_scatt_rates(beams_setup_t beams_setup, atom_t atom, transitions_set_t trans_set, magnetic_field_t B_params);

// Get transition by index
transition_t get_transition_by_idx(int idx, transitions_set_t trans_set);

// Get beam by index
beam_t get_beam_by_idx(int idx, beams_setup_t beams_setup);

// Get polarizations amplitudes
void set_polarizations_amplitudes(beam_t *beam, double *B, double *eB);

// Change-of-matrix basis from polarization frame P to cartisian frame C
complex_t **change_of_basis_from_P_to_C(void);

// Get Zeeman shift (units of Gamma)
double zeeman_shift(double *pos, transition_t trans, magnetic_field_t B_params, int pol);

// Get magnetic acceleration
double *magnetic_acceleration(atom_t atom, transition_t trans, magnetic_field_t B_params);

// Copy a beam structure to avoid the pointer association
beam_t copy_beam(beam_t beam);

// Check whether the atom is inside a defined threshold related with the Zeeman shift
int is_inside_threshold(atom_t atom, transitions_set_t trans_set, magnetic_field_t B_params, beams_setup_t beams_setup);

// Verify whether all transition have the same detuning and saturation parameter
int are_all_beams_similars(beams_setup_t beams_setup);
//--

// Utility
//--
// Concatenate strings
char *str_concatenate(char *str1, char *str2);

// Read lines from a file
char **read_lines(char *path);

// Replace a character in a string
char* str_replace(char *orig, char *rep, char *with);

// Get double array from string in the format [f1 f2 ... fn]
double *get_double_array(char *str, int *size);

// Rotating matrix which rotates vectors by an angle theta [degree] about axis [1 (x), 2 (y), 3 (z)]
double **rotating_matrix(double theta, int axis);

// Convert string to double
double str_to_double(char *str);

// Convert string to int
int str_to_int(char *str);

// Generate a double random number following a Gaussian distribution
double random_norm(double mean, double std_dev);

// Generate a orthonormal basis given a vector which will be the z-direction
double **orthonormal_basis(double *v3);

// Get the largest value of a float array
double max_arr(double *arr, int size);

// Pick an element of an integer array randomly given an array of probabilities
int random_pick(double *probs, int size);

// Get normalized random vector
double *random_r3_vector(void);

// Sampling of exponential distribution
double random_exp(double mean);

// Update histogram
void update_hist(histogram_t *hist, double val);
//--

// Debug
//--
// Print parameters of initial conditions
void print_initial_conditions(initial_conditions_t ini_conds);

// Print parameters of performance
void print_performance(performance_t perform);

// Print parameters of the magnetic field
void print_magnetic_field(magnetic_field_t B_params);

// Print parameters of the transitions
void print_transitions(transitions_set_t trans_set);

// Print parameters of the transitions
void print_beams_setup(beams_setup_t beams_setup);

// Print parameters of the atom
void print_atom(atom_t atom);

// Print simulation progress
void print_progress(atom_t atom, results_t res, magnetic_field_t B_params, int progress);

// Print results
int print_results(results_t res, atom_t atom, int show_hist);
//--