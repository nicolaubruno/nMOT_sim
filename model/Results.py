#
# Libraries and modules
#--
import sys, os, gc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from model.BeamsSetup import BeamsSetup

from scipy.optimize import brentq, curve_fit
from scipy.special import erf
#--

#
class Results:

    #--
    # Attributes 
    #--

    # Identification
    @property
    def id(self):
        return self._id

    # Constants in the international system
    @property
    def ctes(self):
        return self._ctes
 
    #
    @property
    def atom(self):
        return self._atom

    # 
    @property
    def transitions(self):
        return self._transitions

    # Initial conditions
    @property
    def ini_conds(self):
        return self._ini_conds

    # Performance
    @property
    def perform(self):
        return self._perform

    # Magnetic field
    @property
    def B_params(self):
        return self._B_params
    
    # Beams setup
    @property
    def beams_setup(self):
        return self._beams_setup

    # Information about the parameters
    @property
    def info(self):
        return self._info
    
    # Marginal histograms of positions
    @property
    def pos_hist(self):
        return self._pos_hist
    
    # Marginal histograms of velocities
    @property
    def vel_hist(self):
        return self._vel_hist
    
    # Trapped atoms
    @property
    def trapped_atoms(self):
        return self._trapped_atoms

    # Escape velocity
    @property
    def escape_vel(self):
        return self._escape_vel

    # Average escape time
    @property
    def average_escape_time(self):
        return self._average_escape_time

    # Looping status
    @property
    def loop(self):
        return self._loop

    # Directories
    @property
    def dirs(self):
        return self._dirs

    # Cut reference for the trap depth calculation
    @property
    def cut_trap_depth(self):
        return 0.5
    
    # Initial Velocity Distribution
    @property
    def ini_vel(self):
        return self._ini_vel
    

    #--
    # Operational Methods 
    #--

    #
    def __init__(self, code, group, subgroup = None, shortname = None, is_new = False):
        # Constants (SI)
        self._ctes = {
            'u': 1.660539040e-27,\
            'k_B': 1.38064852e-23,\
            'hbar': 1.0544718e-34,\
            'h': 6.626070040e-34
        }

        # Loop variable
        self._loop = {
            "var": None,\
            "values":[],\
            "length": -1,\
            "active": 0
        }

        # Directories
        self._dirs = {
            "model":"model/",\
            "root":"model/results/",\
        }

        self._flip_bins = False

        # Identification
        #--
        self._id = {}
        self._id["code"] = int(code)
        self._id["name"] = shortname.strip() if shortname is not None and len(shortname.strip()) > 0 else None
        self._id["group"] = group
        self._id["subgroup"] = subgroup if (len(subgroup.strip()) > 0 and subgroup != "root") else None
        #--

        # Create a blank result object
        if is_new: 
            self.__new()

        # Get an existing results set
        else:
            # Get name
            self.__get_name()

        # Add directories
        #--
        self._dirs["group"] = self.dirs["root"] + "group_" + self.id["group"] + '/'
        if self.id["subgroup"] is not None: self._dirs["group"] += "subgroup_" + self.id["subgroup"] + '/'

        self._dirs["result"] = self.dirs["group"] + str(self.id["code"])
        if self.id["name"] is not None: self._dirs["result"] += "_" + self.id["name"]
        self._dirs["result"] += "/"
        #--

        # Get loop
        if not is_new: self.__get_loop()

    # Looping method
    #--
    #
    def __get_loop_values(self, loop_str):
        # Return variable
        values = []

        # Check string
        if len(loop_str) > 4 and loop_str[0:4] == 'loop':
            #
            # Loop option 1
            if loop_str[4] == '[' and loop_str[-1] == ']':
                opts = loop_str[5:-1].split(' ')

                val = float(opts[0])
                end = float(opts[1])
                step = float(opts[2])
                
                if ((end - val) < 0 and step < 0) or ((end - val) > 0 and step > 0):
                    values = []

                    while val <= end:
                        values.append(val)
                        val += step
                else:
                    raise ValueError('Incorrect looping in the parameters')

            elif loop_str[4] == '{' and loop_str[-1] == '}':
                values = np.array(loop_str[5:-1].split(' '), dtype=float)
                values = sorted(values, key=(lambda x: float(x)))

            else:
                raise ValueError('Invalid loop variable')

        return values
    
    #
    def __get_loop(self):
        # Sorting
        pos = []

        # Check directories
        res_dir_pter = os.scandir(self.dirs["result"])
        for res_dir in res_dir_pter:
            if self.loop["length"] == -1: 
                loop_var = res_dir.name.split("_")
                if len(loop_var) > 1: 
                    self._loop["var"] = "_".join(loop_var[1:])
                    self._loop["length"] = 0

                else: break

            if self.loop["length"] >= 0:
                self.loop_idx(self.loop["length"])
                self.loop["length"] += 1

                for param in [self.B_params, self.ini_conds, self.perform, self.atom, self.beams_setup.general]:
                    if self.loop["var"] in param.index:
                        self.loop["values"].append(param[self.loop["var"]])

        self.loop_idx(0)
    
    #
    def __set_loop(self):
        # First time setting the loop
        if self.loop["length"] == -1:
            # General parameters
            for param in [self.B_params, self.ini_conds, self.perform, self.atom, self.beams_setup.general]:
                for idx in param.index:
                    if type(param[idx]) is not list:
                        values = self.__get_loop_values(str(param[idx]))
                        if len(values) > 0:
                            self._loop["var"] = idx
                            self._loop["values"] = values
                            self._loop["length"] = len(values)
                            break


        # Change loop value
        #--
        for param in [self.B_params, self.ini_conds, self.perform, self.atom, self.beams_setup.general]:
            if self.loop["var"] in param.index:
                param[self.loop["var"]] = self.loop["values"][self.loop["active"]]
        #--
    
    #
    def loop_idx(self, idx):
        # Change active loop index
        self._loop["active"] = idx

        # Change directory
        self._dirs["active_res"] = self.dirs["result"] 
        self._dirs["active_res"] += "res" + str(self.loop["active"] + 1)
        
        if self.loop["var"] is None: 
            self._dirs["active_res"] += "/"
        else:
            self._dirs["active_res"] += "_" + self.loop["var"] + "/"

        self.__get_attr()
        self.__get_dists()
        self.__get_log()
    #--

    #
    def __get_name(self):
        # Status of result directory
        is_res_dir_exists = False

        # Result directory
        res_dir_path = self.dirs["root"]
        res_dir_path += "group_" + self.id["group"] + "/"
        if self.id["subgroup"] is not None:
            res_dir_path += "subgroup_" + self.id["subgroup"] + "/"

        # Get results directory
        res_dir_pter = os.scandir(res_dir_path)
        for res_dir in res_dir_pter:
            splited_res_dir_name = res_dir.name.split("_")

            if int(splited_res_dir_name[0]) == self.id["code"]:
                self._id["name"] = "_".join(splited_res_dir_name[1:])
                if len(self.id["name"].strip()) == 0: self._id["name"] = None
                is_res_dir_exists = True
                break

        # Check result directory
        if not is_res_dir_exists:
            raise ValueError("Result directory does not exist")
    
    # Get attributes
    def __get_attr(self):
        # Directory of parameters
        params_dir = self.dirs["active_res"] + "parameters/"

        # Atom
        path = params_dir + "atom.csv"
        self._atom = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Transitions
        path = params_dir + "transitions.csv"
        self._transitions = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Beams
        self._beams_setup = BeamsSetup(params_dir + "beams/")

        # Initial conditions
        path = params_dir + "initial_conditions.csv"
        self._ini_conds = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Performance
        path = params_dir + "performance.csv"
        self._perform = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)

        # Magnetic field
        path = params_dir + "magnetic_field.csv"
        self._B_params = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
        self._B_params["B_axial"] = self.__str_to_list(self._B_params["B_axial"], size = 3)
        self._B_params["B_bias"] = self.__str_to_list(self._B_params["B_bias"], size = 3)
        self._B_params["B_lin_grad"] = self.__str_to_list(self._B_params["B_lin_grad"], size = 3)

        # Information about the parameters
        path = self._dirs["model"] + "parameters/informations.csv"
        self._info = pd.read_csv(path, header=0)
        self._info.set_index("parameter", inplace=True)
        self._info.fillna("", inplace=True)

        # Cast values
        self.__cast_params_values()

    # Get distributions
    def __get_dists(self):
        # Initial velocity distribution
        self._ini_vel = -1*np.ones((3, int(self.perform["num_sim"])))

        # Marginal histograms
        self._pos_hist = [{"freqs": None, "dens": None, "bins": None} for i in range(3)]
        self._vel_hist = [{"freqs": None, "dens": None, "bins": None} for i in range(3)]

        # Marginal histograms
        #--
        path = self.dirs["active_res"] + "marginals.csv"
        if os.path.exists(path):
            # Read file
            df = pd.read_csv(path, index_col=0)

            # Check if velocities and positions exist
            check_vel = (('vx' in df.columns) and ('vy' in df.columns) and ('vz' in df.columns))
            check_pos = (('x' in df.columns) and ('y' in df.columns) and ('z' in df.columns))

            #
            # Frequencies
            if check_pos:
                self._pos_hist[0]["freqs"] = np.array(df['x'])
                self._pos_hist[1]["freqs"] = np.array(df['y'])
                self._pos_hist[2]["freqs"] = np.array(df['z'])

            if check_vel:
                self._vel_hist[0]["freqs"] = np.array(df['vx'])
                self._vel_hist[1]["freqs"] = np.array(df['vy'])
                self._vel_hist[2]["freqs"] = np.array(df['vz'])

            #
            # Densities and bins of marginal histograms
            #---
            for i in range(3):
                freq_sum_pos = np.sum(self._pos_hist[i]["freqs"])
                freq_sum_vel = np.sum(self._vel_hist[i]["freqs"])

                #
                # Bins
                self._pos_hist[i]["bins"] = - np.ones(int(int(self.perform['num_bins']))) * float(self.perform['max_r'])
                self._vel_hist[i]["bins"] = - np.ones(int(int(self.perform['num_bins']))) * float(self.perform['max_v'])
                pos_delta = 2*float(self.perform['max_r']) / float(int(self.perform['num_bins']))
                vel_delta = 2*float(self.perform['max_v']) / float(int(self.perform['num_bins']))

                for j in range(int(self.perform['num_bins'])):
                    self._pos_hist[i]["bins"][j] += (j+1)*pos_delta - (pos_delta/2)
                    self._vel_hist[i]["bins"][j] += (j+1)*vel_delta - (vel_delta/2)

                # Densities
                # ---
                if freq_sum_pos > 0 and freq_sum_vel > 0:
                    self._pos_hist[i]["dens"] = self._pos_hist[i]["freqs"] / freq_sum_pos
                    self._vel_hist[i]["dens"] = self._vel_hist[i]["freqs"] / freq_sum_vel

                else:
                    self._pos_hist[i]["dens"] = self._pos_hist[i]["freqs"]
                    self._vel_hist[i]["dens"] = self._vel_hist[i]["freqs"]
                # ---
            #---
        #--

        # Initial Velocities
        #--
        path = self.dirs["active_res"] + 'ini_velocities.csv'
        if os.path.exists(path):
            # Read file
            df = pd.read_csv(path, index_col=0)

            # Check if velocity exists
            check_vel = (('vx' in df.columns) and ('vy' in df.columns) and ('vz' in df.columns))

            #
            # Frequencies
            self._ini_vel[0] = np.array(df['vx'])
            self._ini_vel[1] = np.array(df['vy'])
            self._ini_vel[2] = np.array(df['vz'])

            if check_vel:
                self._ini_vel[0] = np.array(df['vx'])
                self._ini_vel[1] = np.array(df['vy'])
                self._ini_vel[2] = np.array(df['vz'])     
        else: self._ini_vel = np.array([])
        #--
    
    #
    def __get_log(self):
        path = self.dirs["active_res"] + 'log.csv'
        if os.path.exists(path):
            log = pd.read_csv(path, header=0, index_col=0, squeeze=True).astype(object)
            
            # Trapped atoms
            if "trapped_atoms" in log:
                self._trapped_atoms = int(log['trapped_atoms'])

            # Average escape time
            if "average_escape_time" in log: 
                self._average_escape_time = float(log['average_escape_time'])

            # Escape velocity
            if "escape_vel" in log:
                self._escape_vel = float(log["escape_vel"])
    
    #
    def __new(self):
        # Create results directory
        #--
        # Group directory
        dir_path = self.dirs["root"] + "group_" + self.id["group"] + "/"

        # Subgroup directory
        if (self.id["subgroup"] is not None):
            dir_path += "subgroup_" + self.id["subgroup"] + "/"

        # Code
        dir_path += str(self.id["code"])
        if self.id["name"] is not None: dir_path += '_' + self.id["name"]
        dir_path += '/'

        if os.path.exists(dir_path):
            raise ValueError("Results already exists!")

        else: os.mkdir(dir_path)
        #--

        # Create directories for each looping
        #--
        # Create new attributes     
        self.__create_attr()

        # Looping
        num_res = self.loop["length"] if self.loop["length"] > 0 else 1
        for i in range(num_res):
            # Result directory
            if self.loop["var"] is not None:
                res_dir = dir_path + "res" + str(i+1) + '_' + self.loop["var"] + '/'

            else:
                res_dir = dir_path + "res1/"

            # Create directory
            os.mkdir(res_dir)

            # Save parameters of the simulation
            #--
            params_dir = res_dir + "parameters/"
            beams_dir = res_dir + "parameters/beams/"
            os.mkdir(params_dir)
            os.mkdir(beams_dir)

            # Change loop variable 
            if self.loop["var"] is not None:
                self.loop["active"] = i
                self.__set_loop()

            self.__list_to_str_all_params()
            self.atom.to_csv(params_dir + "atom.csv", header="atom")
            self.transitions.to_csv(params_dir + "transitions.csv", header="transitions")
            self.ini_conds.to_csv(params_dir + "initial_conditions.csv", header="initial_conditions")
            self.perform.to_csv(params_dir + "performance.csv", header="performance")
            self.B_params.to_csv(params_dir + "magnetic_field.csv", header="magnetic_field")
            self.__str_to_list_all_params()

            # Beams
            general = pd.Series(self.beams_setup.general.to_dict(), name="general")
            general["sidebands"] = "[" + str(general["sidebands"]["num"]) + " " + str(general["sidebands"]["freq"]) + "]"
            general.to_csv(beams_dir + "general.csv", header="general")
            self.beams_setup.get_setup_dataframe().to_csv(beams_dir + "setup.csv", header="setup")
            #--

    #
    def __create_attr(self):
        # Parameters directory
        self._dirs["active_res"] = self._dirs["model"]

        # Get attribute
        self.__get_attr()

        # Check looping values
        self.__set_loop()

    # Check if name exists
    def __check_name(self, name):
        #
        # Variables
        obj_scandir = os.scandir(self.root_dir)
        ret = False

        for path in obj_scandir:
            str_splited = path.name.split("_")
            code = int(str_splited[0])

            check_name = ""
            for j in range(1, len(str_splited)):
                if j == 1: check_name += str_splited[j]
                else: check_name += '_' + str_splited[j]

            if code == self.code:
                if check_name == name:
                    ret = True
                    break

        return ret  

    # Convert a string to a list with "(int) size" float elements
    def __str_to_list(self, str_arr: str, size: int = 3) -> list:
        return [float(val) for val in str_arr[1:-1].split(' ')]

    # Convert a list with float elements to a string
    def __list_to_str(self, arr: list) -> str:
        str_arr = "["

        for val in arr:
            str_arr += str(val) + " "

        str_arr = str_arr[:-1] + "]"

        return str_arr

    # Convert a list to a list with float elements in all parameters
    def __list_to_str_all_params(self):
        self._B_params["B_axial"] = self.__list_to_str(self._B_params["B_axial"])
        self._B_params["B_bias"] = self.__list_to_str(self._B_params["B_bias"])
        self._B_params["B_lin_grad"] = self.__list_to_str(self._B_params["B_lin_grad"])

    # Convert a list to a list with float elements in all parameters
    def __str_to_list_all_params(self):
        self._B_params["B_axial"] = self.__str_to_list(self._B_params["B_axial"], size = 3)
        self._B_params["B_bias"] = self.__str_to_list(self._B_params["B_bias"], size = 3)
        self._B_params["B_lin_grad"] = self.__str_to_list(self._B_params["B_lin_grad"], size = 3)

    # Cast values of the parameters
    def __cast_params_values(self):
        # Atom
        self._atom['Z'] = int(self.atom['Z'])
        self._atom['mass'] = int(self.atom['mass'])

        # Initial conditions
        self._ini_conds['g_bool'] = int(self.ini_conds['g_bool'])
        self._ini_conds['T_0'] = float(self.ini_conds['T_0'])
        self._ini_conds['v_0'] = float(self.ini_conds['v_0'])

        # Magnetic Field
        self._B_params["B_0"] = float(self.B_params["B_0"])

        # Performance
        self._perform['recording_time'] = float(self.perform['recording_time'])
        self._perform['max_r'] = float(self.perform['max_r'])
        self._perform['max_v'] = float(self.perform['max_v'])
        self._perform['num_sim'] = int(self.perform['num_sim'])
        self._perform['num_bins'] = int(self.perform['num_bins'])
        self._perform['wait_time'] = float(self.perform['wait_time'])
        self._perform['dt'] = float(self.perform['dt'])
        self._perform['parallel_tasks'] = int(self.perform['parallel_tasks'])

        # Transitions
        self._transitions["gamma"] = float(self.transitions["gamma"])
        self._transitions["lambda"] = float(self.transitions["lambda"])
        self._transitions["J_gnd"] = int(self.transitions["J_gnd"])
        self._transitions["J_exc"] = int(self.transitions["J_exc"])
        self._transitions["g_gnd"] = float(self.transitions["g_gnd"])
        self._transitions["g_exc"] = float(self.transitions["g_exc"])

    #--
    # Methods to save data
    #--

    # Add frequencies in the marginal histograms of positions
    def add_marginals(self, pos_freqs_arr, vel_freqs_arr):
        data = {
            'x': pos_freqs_arr[0],\
            'y': pos_freqs_arr[1],\
            'z': pos_freqs_arr[2],\
            'vx': vel_freqs_arr[0],\
            'vy': vel_freqs_arr[1],\
            'vz': vel_freqs_arr[2]
        }

        path = self.dirs["active_res"] + "marginals.csv"
        freqs = pd.DataFrame(data).astype("int32")
        freqs.fillna(0, inplace=True)
        freqs.to_csv(path)

        #
        # Release memory
        del freqs
        del data
        del path

        gc.collect()

    # Add histogram of initial velocities
    def add_ini_vel(self, ini_vel):
        data = {
            'vx':ini_vel[0],\
            'vy':ini_vel[1],\
            'vz':ini_vel[2],\
        }

        path = self.dirs["active_res"] + "ini_velocities.csv"
        df = pd.DataFrame(data).astype(float)
        df.fillna(0, inplace=True)
        df.to_csv(path)

        #
        # Release memory
        del df
        del data
        del path

        gc.collect()

    # Add trapped atoms
    def add_infos(self, data):
        path = self.dirs["active_res"] + "log.csv"
        pd.Series(data, name="log").to_csv(path, header="log")

    #--
    # Methods to processing and view data 
    #--

    # Gaussian fitting of position histograms
    def get_fitting_pos_hist(self):
        # Gaussian function
        def gaussian_f(x, amp, mean, std):
            return amp * np.exp(- ((x - mean) / std)**2 / 2)

        # Gumbel distribution
        def gumbel_f(x, amp, mu, beta):
            z = (x - mu) / beta
            return amp*np.exp(-(z + np.exp(-z)))
        
        # Fitting parameters
        params = []

        # Axis perpendicular to gravity (x, y)
        for i in range(3):
            amp_guess = np.max(self.pos_hist[i]["dens"])
            mean_guess = np.dot(self.pos_hist[i]["bins"], self.pos_hist[i]["dens"])
            std_guess = np.sqrt(np.dot(self.pos_hist[i]["bins"]**2,  self.pos_hist[i]["dens"]) - mean_guess**2)
            if i == 2:
                amp_guess *= np.e
                popt, pcov = curve_fit(gumbel_f, self.pos_hist[i]["dens"], self.pos_hist[i]["bins"], p0=[amp_guess, mean_guess, std_guess])

            else:
                popt, pcov = curve_fit(gaussian_f, self.pos_hist[i]["dens"], self.pos_hist[i]["bins"], p0=[amp_guess, mean_guess, std_guess])
            
            amp, mean, std = popt
            #params.append({"amp": amp, "mean": mean, "std": std})
            params.append({"amp":amp_guess, "mean":mean_guess, "std":std_guess})

        #exit(0)
        return params

    # Gaussian fitting
    def gaussian_fitting(self, data="pos_hist"):
        # Parameters
        params = []

        # Gaussian function
        gaussian_f = lambda x, amp, mean, std: amp * np.exp(- ((x - mean) / std)**2 / 2)

        for i in range(3):
            # Data
            if data == "pos_hist":
                X = self.pos_hist[i]["bins"]
                Y = self.pos_hist[i]["dens"]

            elif data == "vel_hist":                
                X = self.vel_hist[i]["bins"]
                Y = self.vel_hist[i]["dens"]

            amp_guess = np.max(Y)
            mean_guess = np.dot(X, Y)
            std_guess = np.sqrt(np.dot(X**2,  Y) - mean_guess**2)
            popt, pcov = curve_fit(gaussian_f, X, Y, p0=[amp_guess, mean_guess, std_guess])
            params.append(popt)

        return params

    # Escape flux of atoms
    def escape_flux_atoms(self):
        data = 0

        # With looping
        if len(self.loop["var"]) > 0:
            Y = np.zeros(self.loop["length"])
            X = np.array(self.loop["values"], dtype="float")

            for i in range(self.loop["length"]):
                self.loop_idx(i)

                # Get escape flux of atoms
                if self.average_escape_time > 0:
                    Y[i] = 2*np.pi*self.transition["gamma"]*1e3 * (1 - self.trapped_atoms / self.perform["num_sim"]) / (self.average_escape_time)

            data = (X, Y)

        # Without looping
        elif self.average_escape_time > 0: 
            data = 2*np.pi*self.transition["gamma"]*1e3 * (1 - self.trapped_atoms / self.perform["num_sim"]) / self.average_escape_time

        return data

    # Normalized trapped atoms
    def normalized_trapped_atoms(self, pin_loop=False):
        # With looping
        #--
        if self.loop["length"] > 0 and (not pin_loop):
            Y = np.zeros(self.loop["length"])
            X = np.array(self.loop["values"], dtype="float")

            for i in range(self.loop["length"]):
                self.loop_idx(i)
                Y[i] = (self.trapped_atoms / self.perform['num_sim'])

            data = (X, Y)
        #--

        # Without loop
        #--
        else:
            data = self.trapped_atoms / self.perform['num_sim']
        #--

        return data

    # Get all escape velocities
    def all_escape_velocities(self):
        # Data
        data = None

        # With loop
        #--
        if self.loop["length"] > 0:
            X, Y, Z = [], [], []

            for i in range(self.loop["length"]):
                # Change loop
                self.loop_idx(i)

                # Check value of the escape velocity
                if self.escape_vel > 0:
                    X.append(float(self.loop["values"][i]))
                    Y.append(self.escape_vel)
                    Z.append(self.normalized_trapped_atoms(pin_loop=True))

            # Loop values
            X = np.array(X, dtype="float")

            # Escape velocities
            Y = np.array(Y, dtype="float")

            # Normalized trapped atoms
            Z = np.array(Z, dtype="float")

            # Data
            data = (X, Y, Z)
        #--

        # Without loop
        #--
        else: data = self.escape_vel
        #--

        return data

    #
    def centre_of_mass(self, axis=[0,1,2], CSV_file=False):
        num_axis = len(axis)
        data = {"loop_values": None, "r0": None}

        # Set data variables
        if num_axis == 1:
            # With loop
            if self.loop["length"] > 0:
                # Set data parameter
                data["loop_values"] = np.array(self.loop["values"], dtype=float)
                data["r0"] = np.zeros(self.loop["length"])

                # Get values
                for i in range(self.loop["length"]):
                    # Load parameters of the loop i
                    self.loop_idx(i)

                    # Get centre of mass
                    x = self.pos_hist[axis[0]]["bins"]
                    p = self.pos_hist[axis[0]]["dens"]

                    data["r0"][i] = np.dot(x.reshape(1, len(x)), p)[0]

            # Without loop
            else:
                # Get centre of mass
                x = self.pos_hist[axis[0]]["bins"]
                p = self.pos_hist[axis[0]]["dens"]

                data = np.dot(x.reshape(1, len(x)), p)[0]

        elif num_axis > 1:
            # With loop
            if self.loop["length"] > 0:
                data["loop_values"] = np.array(self.loop["values"], dtype=float)
                data["r0"] = np.zeros((self.loop["length"], num_axis))

                # Get values
                for i in range(self.loop["length"]):
                    # Load parameters of the loop i
                    self.loop_idx(i)

                    for idx, ax in enumerate(axis):
                        # Get centre of mass
                        x = self.pos_hist[ax]["bins"]
                        p = self.pos_hist[ax]["dens"]

                        data["r0"][i][idx] = np.dot(x.reshape(1, len(x)), p)[0]

            # Without loop
            else:
                # Set data
                data = np.zeros(num_axis)

                # Check axes
                for idx, ax in enumerate(axis):
                    # Get centre of mass
                    x = self.pos_hist[ax]["bins"]
                    p = self.pos_hist[ax]["dens"]

                    data[idx] = np.dot(x.reshape(1, len(x)), p)[0]

        # CSV file
        if CSV_file:
            if self.loop["length"] > 1:
                dict_data = {(self.info.loc[self.loop["var"], "name"] + " [Gamma]"): data["loop_values"]}
                label = ["x0 [cm]", "y0 [cm]", "z0 [cm]"]

                if num_axis > 1:
                    for idx in axis:
                        dict_data[label[idx]] = data["r0"][:,idx]

                else: dict_data[label[axis[0]]] = data["r0"]

            else: dict_data = {"r0 [cm]": data}
            data = dict_data

        return data

    #
    def average_velocity(self, axis=[0,1,2], fixed_loop_idx = False):
        #
        # Returns the best parameters to fit a Gaussian function
        def fit_gaussian(x, y):
            # Convert to numpy array
            x = np.array(x)
            y = np.array(y)

            #
            # Gaussian function
            def gaussian(x, mean, std_dev):
                return np.max(y) * np.exp(-((x - mean)/std_dev)**2 / 2)

            #
            # Guess values
            #--
            guess_values = np.zeros(2)
            guess_values[0] = np.sum(x*y) / np.sum(y)
            guess_values[1] = np.sqrt(np.sum(y * (x - guess_values[0])**2) / np.sum(y))
            #--
            #--

            popt, pcov = curve_fit(gaussian, x, y, guess_values)

            return popt

        #
        # Without looping
        #--
        if len(self.loop["var"]) == 0 or fixed_loop_idx:
            # Variables
            if len(axis) == 1:
                vel_c, std_vel_c = fit_gaussian(self.vel_hist[axis[0]]["bins"], self.vel_hist[axis[0]]["dens"])
            
            else:
                vel_c = np.zeros(len(axis))
                std_vel_c = np.zeros(len(axis))

                #
                # Fit a normal Gaussian function
                for idx, val in enumerate(axis):
                    vel_c[idx], std_vel_c[idx] = fit_gaussian(self.vel_hist[val]["bins"], self.vel_hist[val]["dens"])
        #--

        #
        # With looping
        #--
        else:
            #
            # Set initial variables
            #--
            if len(axis) == 1:
                vel_c = np.zeros(len(self.loop["values"]))
                std_vel_c = np.zeros(len(self.loop["values"]))
            
            else:
                vel_c = np.zeros((len(axis),len(self.loop["values"])))
                std_vel_c = np.zeros((len(axis),len(self.loop["values"])))
            #--

            # Mass centre for each looping value
            for i in range(len(self.loop["values"])):
                #
                # Get looping values
                self.loop_idx(i)

                if len(axis) > 1:
                    for j in axis:
                        vel_c[j][i], std_vel_c[j][i] = fit_gaussian(self.vel_hist[j]["bins"], self.vel_hist[j]["dens"])

                else:
                    x = self.vel_hist[axis[0]]["bins"]
                    p = self.vel_hist[axis[0]]["dens"]

                    vel_c[i], std_vel_c[i] = fit_gaussian(self.vel_hist[axis[0]]["bins"], self.vel_hist[axis[0]]["dens"])
        #--

        return vel_c, std_vel_c      

    # Get temperatures in uK
    def temperature(self, fixed_loop_idx = False, method=0, CSV_file = False):
        # Without looping
        # ---
        if len(self.loop["var"]) == 0 or fixed_loop_idx:
            v_av, v_dev = self.average_velocity(axis=[0, 1, 2], fixed_loop_idx=fixed_loop_idx)
            
            if method == 0:
                temp = ((np.sum(v_dev*1e-2)/3)**2 * float(self.atom['mass']) * self.ctes['u']) / self.ctes['k_B']

            elif method == 1:
                v_av = v_av*1e-2
                v_var = (v_dev*1e-2)**2
                v_square = v_var + v_av**2

                temp = (np.sum(v_square) * float(self.atom['mass']) * self.ctes['u']) / (3*self.ctes['k_B'])
            
            else:
                raise ValueError("Invalid method")
        # ---

        # With looping
        # ---
        else:
            v_av, v_dev = self.average_velocity(axis=[0, 1, 2])

            if method == 0:
                temp = ((np.sum(v_dev*1e-2, axis=0)/3)**2 * float(self.atom['mass']) * self.ctes['u']) / self.ctes['k_B']

            elif method == 1:
                v_av = v_av*1e-2
                v_var = (v_dev*1e-2)**2
                v_square = v_var + v_av**2

                temp = (np.sum(v_square, axis=0) * float(self.atom['mass']) * self.ctes['u']) / (3*self.ctes['k_B'])
        # ---

        # CSV file
        # ---
        if CSV_file:
            if self.loop["length"] > 1:
                loop_label = self.info.loc[self.loop["var"], "name"]
                if len(self.info.loc[self.loop["var"], "unit"]) > 0: 
                    loop_label += self.info.loc[self.loop["var"], "unit"]

                data = {loop_label: np.array(self.loop["values"], dtype=float) , "T [mu K]": temp*1e6}
            else: data = None
        else: data = temp
        # ---

        return data   

    # Doppler temperature
    def doppler_temperature(self, power_broadening=False, fixed_loop_idx = False):
        if power_broadening:
            alpha = np.sqrt(1 + self.beams['main']['s_0'])
        else:
            alpha = 0;

        #
        # Check looping
        if self.loop["var"] == "gamma" and fixed_loop_idx:
            temp = np.zeros(len(self.loop["values"]))
            for i, gamma in enumerate(self.loop["values"]):
                temp[i] = 1e9*(self.ctes['hbar'] * gamma * alpha) / (2 * self.ctes['k_B']) # uK

        else:
            temp = 1e9*(self.ctes['h'] * self.transition['gamma'] * alpha) / (2 * self.ctes['k_B']) # uK

        return temp

    # Trapped atoms
    def all_trapped_atoms(self, fixed_loop_idx = False):
        #
        # Without looping
        #--
        if len(self.loop["var"]) == 0 or fixed_loop_idx:
            res = self.trapped_atoms

        #
        # With looping
        #--
        else:
            res = np.zeros(len(self.loop["values"]))

            for i in range(len(self.loop["values"])):
                self.loop_idx(i)
                res[i] = self.trapped_atoms
        #--

        return res

    # Trap Depth
    def trap_depth(self):
        X, Y = np.zeros(self.loop["length"] * np.ones(2))

        print(X, X.shape)
        print(Y, Y.shape)

        exit(0)

        return (X, Y)

    # Capture velocity
    def capture_velocity(self, fit_func = "erf"):
        # Check loop variable
        if self.loop["var"] == 'v_0' and self.loop["length"] > 1:
            # Get velocities and trapped atoms ratio
            #--
            vel = np.array(self.loop['values'], dtype="float")
            ratio = np.zeros(self.loop['length'])

            for i, val in enumerate(self.loop["values"]):
                self.loop_idx(i)
                ratio[i] = self.trapped_atoms

            ratio = ratio / self.perform["num_sim"]
            max_vel = np.max(vel)
            min_vel = np.min(vel)
            #--

            # Get capture velocity
            #--
            # Polynomial fitting
            if fit_func == "poly":
                # Fit polynomial
                fit_params = np.polyfit(vel, ratio, 10, full = False, cov = False)

                # Polynomial function
                def f(x):
                    y = -self.cut_trap_depth

                    for i in range(11):
                        y += fit_params[10 - i]*x**i

                    return y

                # Get capture velocity
                if f(min_vel) > 0 and f(min_vel)*f(max_vel) < 0: 
                    vel_c = brentq(f, min_vel, max_vel, full_output = False)
                else: vel_c = -1

            # Erf function fitting
            elif fit_func == "erf":
                # General complementary error function
                def general_erfc(t, mean, std_dev, amp):
                    return amp*(1 - (erf((t - mean) / np.sqrt(2 * std_dev**2)) - erf((- mean) / np.sqrt(2 * std_dev**2))) / 2)

                # Get data
                fit_params, covs = curve_fit(general_erfc, vel, ratio, bounds=([min(vel), 0, 0], [max(vel), (max(vel) - min(vel)), 1]))
                f = lambda x: general_erfc(x, fit_params[0], fit_params[1], fit_params[2]) - self.cut_trap_depth
                if f(min_vel) > 0 and f(min_vel)*f(max_vel) < 0:
                    vel_c = brentq(f, min_vel, max_vel, full_output = False)
                else: vel_c = -1
            #--

        else: raise ValueError("The loop variable must be v_0 to calculate the capture velocity")

        return vel_c, fit_params

    # Capture temperature
    def capture_temperature(self, fit_func = "erf"):
        T_mean, T_std_dev = (0,0)

        if self.loop["var"] == 'T_0' and len(self.loop["values"]) > 1:
            # Get temperatures and trapped atoms ratio
            #--
            T = np.array(self.loop['values'], dtype="float")
            ratio = np.zeros(self.loop['length'])

            for i, val in enumerate(self.loop["values"]):
                self.loop_idx(i)
                ratio[i] = self.trapped_atoms

            ratio = ratio / self.perform["num_sim"]
            max_T = np.max(T)
            min_T = np.min(T)
            #--

            # Get capture temperature
            #--
            # Polynomial fitting
            if fit_func == "poly":
                # Fit polynomial
                fit_params = np.polyfit(T, ratio, 10, full = False, cov = False)

                # Polynomial function
                def f(x):
                    y = -self.cut_trap_depth

                    for i in range(11):
                        y += fit_params[10 - i]*x**i

                    return y

                # Get capture velocity
                if f(min_T) > 0 and f(min_T)*f(max_T) < 0: T_c = brentq(f, min_T, max_T, full_output = False)
                else: T_c = -1

            # Erf fitting
            elif fit_func == "erf":
                # General complementary error function
                 # General complementary error function
                def general_erfc(t, mean, std_dev, amp):
                    return amp*(1 - (erf((t - mean) / np.sqrt(2 * std_dev**2)) - erf((- mean) / np.sqrt(2 * std_dev**2))) / 2)

                # Get data
                fit_params, covs = curve_fit(general_erfc, T, ratio, bounds=([min_T, 0, 0], [max_T, (max_T - min_T), 1]))
                f = lambda x: general_erfc(x, fit_params[0], fit_params[1], fit_params[2]) - self.cut_trap_depth
                if f(max_T)*f(0) < 0: T_c = brentq(f, 0, max_T, full_output = False)
                else: T_c = -1
            
            # Gaussian fitting
            elif fit_func == "gaussian":
                f = lambda x, amp, std: amp * np.exp(-x**2 / (2 * std**2))

                # Get fitting parameters
                fit_params, covs = curve_fit(f, T, ratio, bounds=([0, 0], [1, max_T]))

                # Get cut off temperature
                #if f(max_T)*f(0) < 0: T_c = brentq(f, 0, max_T, full_output = False)
                #else: T_c = -1
                T_c = -1

            #--

        return T_c, fit_params

    # Get 2D-histogram of positions removing an axis
    def pos_2Dhist(self, axis = 0, val = 0):
        #
        # Get bin index
        bin_idx = 0
        for idx, bin_size in enumerate(self.pos_3Dhist["bins"][axis]):
            if idx > 0 and (float(val) <= bin_size):
                bin_idx = idx - 1
                break

        #
        # Get bins
        if axis == 0:
            axis_label = {'y', 'z'}
            hist = np.zeros((len(self.pos_3Dhist["bins"][1]), len(self.pos_3Dhist["bins"][2])))

        elif axis == 1:
            axis_label = {'x', 'z'}
            hist = np.zeros((len(self.pos_3Dhist["bins"][0]), len(self.pos_3Dhist["bins"][2])))

        elif axis == 2:
            axis_label = {'x', 'y'}
            hist = np.zeros((len(self.pos_3Dhist["bins"][0]), len(self.pos_3Dhist["bins"][1])))

        #
        # Get densities
        for i in range(len(hist)):
            for j in range(len(hist[i])):
                if axis == 0:
                    hist[i][j] = self.pos_3Dhist["dens"][bin_idx,i,j]

                elif axis == 1:
                    hist[i][j] = self.pos_3Dhist["dens"][i,bin_idx,j]

                elif axis == 2:
                    hist[i][j] = self.pos_3Dhist["dens"][i,j,bin_idx]

        return hist
        path = self.directory + "log.csv"
        pd.Series(data, name="log").to_csv(path, header="log")

    # Get Initial Velocity Histogram
    def get_ini_vel_hist(self, axis=0):
        # Get histogram data
        n, bins = np.histogram(self.ini_vel[axis], bins=int(0.05*self.perform["num_sim"]))

        # Build X and Y values
        #--
        delta = 0
        for j in range(len(bins) - 1):
            delta += bins[j+1] - bins[j] 

        delta = delta / (len(bins) - 1)
        X = np.array(bins[0:-1]) + delta / 2
        Y = np.array(n) / len(self.ini_vel[axis])
        #--

        # Gaussian function
        f = lambda x, amp, mean, std: amp * np.exp(-((x - mean)**2 / (2*std**2)))
        guess_std = np.sqrt(self.ctes["k_B"] * self.ini_conds["T_0"]*1e-3 / (self.atom["mass"] * self.ctes["u"]))*1e2
        fit_params, pcov = curve_fit(f, X, Y, p0=[1 / np.sqrt(2*np.pi*guess_std**2), 0, guess_std])

        return fit_params, X, Y