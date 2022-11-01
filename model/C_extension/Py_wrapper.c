//
// Libraries
//

#include "Py_wrapper.h"

// Read simulation for a single atom
static PyObject* C_simulate_atom(PyObject *self, PyObject *args){
    //
    // Variables
    int i, opt;
    long time;
    PyObject *ret;
    results_t res;
    char *params_path;
    PyObject *pos_freqs, *vel_freqs, *ini_vel;

    //
    // Get parameters path
    if(!PyArg_ParseTuple(args, "sil", &params_path, &opt, &time)){
        return NULL;
    }

    // Results
    //--
    res = simulate_atom(params_path, opt, time);
    
    ini_vel = PyList_New(3);
    for(i = 0; i < 3; i++)
        PyList_SetItem(ini_vel, i, Py_BuildValue("d", res.ini_vel[i]));
    

    // Histograms
    //--
    // Build values
    pos_freqs = build_freqs(res.pos_hist);
    vel_freqs = build_freqs(res.vel_hist);

    //
    // Release memory
    for(i = 0; i < 3; i++){
        free(res.pos_hist[i].freqs);
        free(res.vel_hist[i].freqs);
    }

    free(res.pos_hist);
    free(res.vel_hist);
    //--
    //--

    // Return
    ret = Py_BuildValue("OOdidO", pos_freqs, vel_freqs, res.time, res.trapped_atom, res.escape_vel, ini_vel);

    return ret;
}

// Initialize module
PyMODINIT_FUNC PyInit_mot_sim(void){
    return PyModule_Create(&mot_sim_module);
}

// Convert the results in a PyObject list
PyObject *build_freqs(histogram_t *hist){
    //
    // Variables
    //

    int i, j, dim;
    PyObject *list, *item, *item2;

    //
    // Build list
    //

    list = PyList_New(3);

    for(i = 0; i < 3; i++){
        dim = hist[i].num_bins;
        item = PyList_New(dim);

        for(j = 0; j < dim; j++){
            item2 = Py_BuildValue("i", hist[i].freqs[j]);
            PyList_SetItem(item, j, item2);
        }
        
        PyList_SetItem(list, i, item);
    }

    // Return
    return list;
}