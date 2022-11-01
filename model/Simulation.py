# Libraries and modules
#--
import sys, os, gc
import numpy as np
import pandas as pd
import mot_sim as C_ext

from multiprocessing import Pool, cpu_count
from datetime import datetime as dt
from model.Results import Results
from random import randint
#--

# Wrapper function for pool
def simulate_atom(args):
    params_dir, opt, s = args
    return C_ext.simulate_atom(params_dir, opt, s)

#
class Simulation:

    #--
    # Attributes
    #--

    # Escape velocity
    @property
    def escape_vel(self):
        return self._escape_vel
    
    # Initial velocity
    @property
    def ini_vel(self):
        return self._ini_vel
    
    # Average escape time
    @property
    def average_escape_time(self):
        return self._average_escape_time

    # Trapped atoms
    @property
    def trapped_atoms(self):
        return self._trapped_atoms

    # (Array) Position histogram array
    @property
    def pos_freqs_arr(self):
        return self._pos_freqs_arr

    # (Array) Velocity histogram array
    @property
    def vel_freqs_arr(self):
        return self._vel_freqs_arr
    
    # Atoms simulated
    @property
    def atoms_simulated(self):
        return self._atoms_simulated

    # Simulation option
    @property
    def option(self):
        return self._option
    
    # Parallel tasks
    @property
    def parallel_tasks(self):
        return self._parallel_tasks

    # (Object) Results object
    @property
    def results(self):
        return self._results    

    #--
    # Operational methods
    #--

    #
    def __init__(self):
        # Set-up initial values
        self._atoms_simulated = -1
        self._results = None
        self._parallel_tasks = 0

    #
    def new(self, group, subgroup = None, shortname = None, opt = 1):
        # Create a results object
        #--
        self._results = Results(int(dt.now().timestamp()), group, subgroup, shortname.strip(), is_new = True)
        self._flip_bins = True
        #--

        # Check simulation option
        #--
        available_opts = {
            1 : "Atoms start from the origin"
        }

        # Set option
        self._option = opt if (opt in available_opts.keys()) else 1

        # Release memory
        del available_opts
        #--

        # Frequencies of positions and velocities
        self._pos_freqs_arr = np.zeros((3, self.results.perform["num_bins"]))
        self._vel_freqs_arr = np.zeros((3, self.results.perform["num_bins"]))

        # Other informations
        #--
        self._trapped_atoms = 0
        self._atoms_simulated = 0
        self._average_escape_time = 0
        self._escape_vel = 0
        self._ini_vel = -1*np.ones((3, int(self.results.perform["num_sim"])))
        #--

        # Parallel tasks
        self._parallel_tasks = int(self.results.perform["parallel_tasks"])

        # Release memory
        gc.collect()

    #
    def run(self):
        # Check simulation status
        #--
        if self.atoms_simulated < self.results.perform["num_sim"]:
            #
            # Check number of executions
            if (self.atoms_simulated + self.parallel_tasks) > self.results.perform["num_sim"]:
                times = 1
            else:
                times = self.parallel_tasks

            # 
            # Arguments to the pool 
            args_pool = []
            for i in range(times):
                seed = int(randint(0, 1e10))
                args_pool.append((self.results.dirs["active_res"] + "parameters/", self.option, seed))

            #
            # Parallel execution
            with Pool(self.parallel_tasks, maxtasksperchild=100) as pool:
                res = pool.map(simulate_atom, args_pool)
                
                # Add frequencies
                for k in range(times):
                    pos_freqs = res[k][0]
                    vel_freqs = res[k][1]
                    self._trapped_atoms += res[k][3]
                    if res[k][3] == 0: self._average_escape_time += res[k][2]
                    self._escape_vel += res[k][4]

                    # Add position and velocity frequencies
                    #--
                    if self.option == 1 or self.option == 2:                
                        for i in range(3):
                            for j in range(self.results.perform["num_bins"]):
                                self._pos_freqs_arr[i][j] += pos_freqs[i][j]
                                self._vel_freqs_arr[i][j] += vel_freqs[i][j]
                    #--

                    # Add initial velocities
                    for i in range(3): 
                        self._ini_vel[i][self.atoms_simulated + k] = res[k][5][i]

                #
                # Release memory
                del pos_freqs
                del vel_freqs
                del args_pool
                del seed

                #
                # Finish pool and release memory
                pool.terminate()

            #
            # Update atoms simulated
            self._atoms_simulated += times

            # Release memory
            gc.collect()

            return times
        #--

    #
    def change_loop(self, idx):
        # Get results
        self.results.loop_idx(idx)

        #
        # Simulate marginal distributions (Lightweight option)
        #--
        if self.option == 1 or self.option == 2:
            #
            # Frequencies (3D-array)
            del self._pos_freqs_arr
            del self._vel_freqs_arr

            self._pos_freqs_arr = np.zeros((3, self.results.perform["num_bins"]))
            self._vel_freqs_arr = np.zeros((3, self.results.perform["num_bins"]))
        #--

        #
        # Simulate 3D distribution (Heavy option)
        #--
        elif self.option == 3:
            #
            # Frequencies (1D-array)
            del self._pos_freqs_arr
            del self._vel_freqs_arr

            self._pos_freqs_arr = np.zeros(self.results.perform["num_bins"]**3)
            self._vel_freqs_arr = np.zeros(self.results.perform["num_bins"]**3)
        #--

        # Other informations
        self._escape_vel = 0
        self._average_escape_time = 0
        self._trapped_atoms = 0
        self._atoms_simulated = 0

        # Parallel tasks
        self._parallel_tasks = int(self.results.perform["parallel_tasks"])

        # Release memory
        gc.collect()

    #
    def save(self):
        # Add position and velocities in marginal histograms
        self.results.add_marginals(self.pos_freqs_arr, self.vel_freqs_arr)

        # Set infos
        if self.trapped_atoms == self.results.perform["num_sim"]:
            self._average_escape_time = -1
            self._escape_vel = -1

        else: 
            self._average_escape_time = self.average_escape_time / (self.results.perform["num_sim"] - self.trapped_atoms)
            self._escape_vel = self.escape_vel / (self.results.perform["num_sim"] - self.trapped_atoms)

        # Add informations
        self.results.add_infos({
                "trapped_atoms": self.trapped_atoms,\
                "average_escape_time": self.average_escape_time,\
                "escape_vel": self.escape_vel
            })

        # Initial velocity histogram
        self.results.add_ini_vel(self.ini_vel)

        # Flip bins
        self._flip_bins = False

        # Release memory
        gc.collect()

    #--
    # Utility methods
    #--

    # Get available groups
    def available_groups(self):
        #
        # Variables
        i = 2
        res = {1:"root"}

        #
        # List all results groups directories
        #--
        groups_dir = os.scandir("model/results/")
        for group_dir in groups_dir:
            # Check if the results group is valid
            if group_dir.is_dir():
                str_splited = group_dir.name.split("_")

                if(str_splited[0] == "group"):
                    name = ""
                    for j in range(1, len(str_splited)):
                        if j == 1: name += str_splited[j]
                        else: name += '_' + str_splited[j]

                    res[i] = name
                    i += 1
        #--

        return res

    # Get available subgroups
    def available_subgroups(self, group):
        # Available subgroups
        i = 1
        subgroups = {i: "root"}


        # Check group
        if (group is None) or (group == "root"): subgroup = None

        # Get subgroups
        else:
            # Check all directories
            #--
            all_subgroup_dir = os.scandir("model/results/group_" + group)

            for subgroup_dir in all_subgroup_dir:
                splited_str = subgroup_dir.name.split("_")

                if splited_str[0] == "subgroup":
                    i += 1
                    subgroups[i] = "_".join(splited_str[1:])
            #--

        return subgroups

    # Get available results
    def available_results(self, results_id):
        # Results
        results = {}

        # Get results
        #--
        # Get results directory
        #---
        path = "model/results/group_" + results_id["group"]
        if (results_id["subgroup"] is not None) and (results_id["subgroup"] != "root"): 
            path += "/subgroup_" + results_id["subgroup"]
        
        results_dir_pter = os.scandir(path)
        #---

        # Check all directories
        for res_dir in results_dir_pter:

            # Check if the results is valid
            if res_dir.is_dir() and (len(res_dir.name) > 6) and res_dir.name[:5].isnumeric():
                
                # Check if the results are complete
                is_valid_res = False
                loop_dir_pter = os.scandir(res_dir.path)
                for loop_dir in loop_dir_pter:
                    is_valid_res = False

                    if loop_dir.is_dir() and loop_dir.name[:3] == "res":
                        loop_file_pter = os.scandir(loop_dir.path)
                        for loop_file in loop_file_pter:
                            if loop_file.name == "log.csv": 
                                is_valid_res = True
                                break

                    if not is_valid_res: break

                # Append valid results
                if is_valid_res:
                    splited_str = res_dir.name.split("_")
                    results[int(splited_str[0])] = "_".join(splited_str[1:])
        #--
        
        return results

    # Check if results code exists
    def check_results_code(self, code):
        #
        # Variables
        dir_path = "model/results/"
        obj_scandir = os.scandir(dir_path)
        ret = False

        for path in obj_scandir:
            str_splited = path.name.split("_")
            sim_code = int(str_splited[0])

            name = ""
            for j in range(1, len(str_splited)):
                if j == 1: name += str_splited[j]
                else: name += '_' + str_splited[j]

            if sim_code == int(code):
                ret = True
                break

        return ret  