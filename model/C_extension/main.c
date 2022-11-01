//
// Libraries
//

#include "mot_sim.h"
#include <time.h>

void main(){
    int i, j;
    results_t res;
    
    res = simulate_atom("../parameters/", 1, time(0));
    r3_print(res.ini_vel, "Initial Velocity [m/s]");
    printf("trapped_atom = %d\n", res.trapped_atom);
    printf("time [ms] = %f\n", res.time);
}