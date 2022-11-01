#
# Libraries and modules
import os, sys, time
import pandas as pd
from matplotlib import pyplot as plt, ticker, rc
from scipy.special import erf
import seaborn as sns
import numpy as np

from datetime import datetime as dt
from model import Results
from tqdm import tqdm

#
class View:
    
    #--
    # Attributes
    #--

    #
    @property
    def separator(self):
        return self._separator
    
    #
    @property
    def header(self):
        return self._header
    
    #
    @property
    def result_id(self):
        return self._result_id

    #-- 
    # Operational Methods 
    #--

    #
    def __init__(self, simulation):
        # Model objects
        self.__simulation = simulation

        # Terminal separator line
        #--
        self._separator = '\n'
        for i in range(60): self._separator += '-'

        self._header = "\nMonte Carlo simulation of narrow-line magneto-optical traps\n"
        self._header += "Version 2.0, Author: Bruno N. Santos;" + self._separator
        self._header += "\nGlobal options: -1 -> Back | 0 -> Exit" + self.separator + '\n'
        #--

        # Result identification
        self._result_id = {
            "id": None,\
            "name": None,\
            "group": None,\
            "subgroup": None
        }

        # Matplotlib settings
        rc("font", size=14)
        rc("axes", labelsize=14)
        rc("xtick", labelsize=14)
        rc("ytick", labelsize=14)

    # Set result identification
    def set_result_id(self, result_id):
        self._result_id["code"] = result_id["code"]
        self._result_id["name"] = result_id["name"]
        self._result_id["group"] = result_id["group"]
        self._result_id["subgroup"] = result_id["subgroup"]

    #--
    # Data visualization
    #--

    # View histogram of positions
    def pos_histograms(self, res):
        # Message
        self.print_msg('Processing data ...', press_button=False)

        # Fitting
        fit_params = res.gaussian_fitting(data="pos_hist")
        gaussian_f = lambda x, amp, mean, std: amp * np.exp(- ((x - mean) / std)**2 / 2)

        # Set plot
        plt.clf()
        plt.style.use('seaborn-darkgrid') # Theme
        plt.rcParams.update({
                "font.size":14,\
                "axes.titlepad":14
            })

        # Subplots
        fig, ax = plt.subplots(1, 3, sharey=True)
        fig.set_size_inches(13, 4)
        labels = ["v_x", "v_y", "v_z"]
        
        for i in range(3):    
            delta_bins = res.pos_hist[i]["bins"][1] - res.pos_hist[i]["bins"][0] 

            ax[i].bar(res.pos_hist[i]["bins"], res.pos_hist[i]["dens"], width=delta_bins, label="Data")
            x = np.linspace(res.pos_hist[i]["bins"][0], res.pos_hist[i]["bins"][-1], 1000)
            label = (r"$\mu = %f cm$" + "\n" + r"$\sigma = %f cm$") % (fit_params[i][1], fit_params[i][2])
            ax[i].plot(x, [gaussian_f(xi, fit_params[i][0], fit_params[i][1], fit_params[i][2]) for xi in x], linestyle="--", marker="", color="black", label=label)
            ax[i].set_xlabel(r"$" + labels[i] + " [cm / s]" + r"$")
            ax[i].legend(frameon=True, loc="lower right", prop={'size':12})
            ax[i].grid(linestyle="--")
            ax[i].set_aspect('auto')

        plt.close(1)
        plt.tight_layout()
        plt.show()
        self.print_msg('Showing graph ...', press_button=False)

    # View histogram of positions
    def vel_histograms(self, res):
        # Message
        self.print_msg('Processing data ...', press_button=False)

        # Fitting
        fit_params = res.gaussian_fitting(data="vel_hist")
        gaussian_f = lambda x, amp, mean, std: amp * np.exp(- ((x - mean) / std)**2 / 2)

        # Set plot
        plt.clf()
        plt.style.use('seaborn-darkgrid') # Theme
        plt.rcParams.update({
                "font.size":14,\
                "axes.titlepad":14
            })

        # Subplots
        fig, ax = plt.subplots(1, 3, sharey=True)
        fig.set_size_inches(13, 4)
        labels = ["v_x", "v_y", "v_z"]
        
        for i in range(3):    
            delta_bins = res.vel_hist[i]["bins"][1] - res.vel_hist[i]["bins"][0] 

            ax[i].bar(res.vel_hist[i]["bins"], res.vel_hist[i]["dens"], width=delta_bins, label="Data")
            x = np.linspace(res.vel_hist[i]["bins"][0], res.vel_hist[i]["bins"][-1], 1000)
            label = (r"$\mu = %f cm/s$" + "\n" + r"$\sigma = %f cm/s$") % (fit_params[i][1], fit_params[i][2])
            ax[i].plot(x, [gaussian_f(xi, fit_params[i][0], fit_params[i][1], fit_params[i][2]) for xi in x], linestyle="--", marker="", color="black", label=label)
            ax[i].set_xlabel(r"$" + labels[i] + " [cm / s]" + r"$")
            ax[i].legend(frameon=True, loc="lower right", prop={'size':12})
            ax[i].grid(linestyle="--")
            ax[i].set_aspect('auto')

        plt.close(1)
        plt.tight_layout()
        plt.show()
        self.print_msg('Showing graph ...', press_button=False)

    # Initial velocity distribution
    def ini_vel_hist(self, res):
        # Message
        self.print_msg('Processing data ...', press_button=False)

        # Set plot
        plt.clf()
        plt.style.use('seaborn-darkgrid') # Theme
        plt.rcParams.update({
                "font.size":14,\
                "axes.titlepad":14
            })
        plt.suptitle("Histogram of Initial Velocities")

        # Subplots
        fig, ax = plt.subplots(1, 3)
        fig.set_size_inches(15, 5)
        labels = ["x", "y", "z"]
        legend_pos = [0.33, 0.67, 0.985]
        lg = []

        for i in range(3):
            # Get histogram
            fit_params, X, Y = res.get_ini_vel_hist(axis=i)

            # Plot histogram
            ax[i].bar(X, height=Y, width=1.0*(X[1] - X[0]))

            # Plot fitting
            X_fit = np.linspace(min(X)*1.1, max(X)*1.1, 1000)
            f = lambda x: fit_params[0] * np.exp(-((x - fit_params[1])**2 / (2 * fit_params[2]**2)))
            label_fitting = "Gaussian Fitting\n"
            label_fitting += ("mean = %.1f cm/s, std = %.1f cm/s \n" % (fit_params[1], fit_params[2]))
            label_fitting += "Estimated Temperature: %.1f mK\n" % ((float(res.atom["mass"])*res.ctes["u"]*(fit_params[2]*1e-2)**2 / res.ctes["k_B"])*1e3)
            label_fitting += "Expected Temperature: %.1f mK" % res.ini_conds["T_0"]
            ax[i].plot(X_fit, list(map(f, X_fit)), color="black", marker="", linestyle="-", label=label_fitting)

            # Set labels
            ax[i].set_title("Histogram of velocities in the " + labels[i] + "-axis", fontsize=16)
            ax[i].set_xlabel("v [cm / s]")
            #ax[i].set_yticklabels([])
            ax[i].legend(frameon=True, prop={'size':14}, loc="lower left", bbox_to_anchor=(0.05, -0.57), ncol=2)
            ax[i].grid(linestyle="--")
            ax[i].set_aspect('auto')

        plt.close(1)
        plt.tight_layout(rect=[0.0, 0.18, 1.0, 1.0])
        plt.show()
        self.print_msg('Showing graph ...', press_button=False)

    # Heatmap of positions
    def pos_heatmap(self, res):
        self.print_msg('Processing data ...', press_button=False)

        grids = []
        labels = []
        idxs = [[0,1], [0,2]]

        for i in range(2):
            # Centre of Mass and Standard Deviation
            cm, std = res.centre_of_mass(axis=idxs[i], pin_loop = True)

            # Get highest centre of mass
            if(abs(max(cm)) > abs(min(cm))): highest_cm = abs(max(cm))
            else: highest_cm = abs(min(cm))

            # Get highest standard deviation
            if(abs(max(std)) > abs(min(std))): highest_std = abs(max(std))
            else: highest_std = abs(min(std))

            grid_len = len(res.pos_hist[idxs[i][0]]["bins"])
            ini_bin = -1
            final_bin = -1
            pter = (highest_cm + 5*highest_std)

            for j in range(grid_len):
                if (ini_bin < 0) and (pter > abs(res.pos_hist[idxs[i][0]]["bins"][j])):
                    ini_bin = j

                elif (ini_bin > 0) and (pter < abs(res.pos_hist[idxs[i][0]]["bins"][j])): 
                    final_bin = j
                    break

            if ini_bin > final_bin:
                ini_bin = 0
                final_bin = grid_len - 1

            grid_len = (final_bin - ini_bin + 1)
            grids.append(np.zeros([grid_len, grid_len]))
            rows = []
            cols = []

            for j in range(grid_len):
                for k in range(grid_len):
                    grids[-1][j][k] = res.pos_hist[idxs[i][1]]["dens"][final_bin-j]*res.pos_hist[idxs[i][0]]["dens"][ini_bin+k]

                rows.append("%.2f" % (float(res.pos_hist[idxs[i][0]]["bins"][final_bin-j])))
                cols.append("%.2f" % (float(res.pos_hist[idxs[i][0]]["bins"][ini_bin+j])))

            grids[-1] = pd.DataFrame(data = grids[-1], index=rows, columns=cols)

        # Set plot
        plt.clf()
        plt.style.use('seaborn-darkgrid') # Theme

        # Subplots
        fig, ax = plt.subplots(1, 2)
        fig.set_size_inches(8, 4)

        # Grid XY
        #--
        sns.set_style("darkgrid")
        g1 = sns.heatmap(grids[0], ax=ax[0], vmin = 0, cbar = False, cmap="Blues", square=True, xticklabels=int(grids[0].shape[0]/5), yticklabels=int(grids[0].shape[0]/5))
        
        # Settings labels
        g1.set_xlabel("x [cm]")
        g1.set_ylabel("y [cm]")
        g1.set_yticklabels(g1.get_yticklabels(), rotation = 0, fontsize=12)
        g1.set_xticklabels(g1.get_xticklabels(), fontsize=12)

        # Borders
        g1.axhline(y = 0, color='k', linewidth=3)    
        g1.axhline(y = grids[0].shape[1], color='k', linewidth=3)    
        g1.axvline(x = 0, color='k', linewidth=3)    
        g1.axvline(x = grids[0].shape[0], color='k', linewidth=3)    
        #--

        # Grid YZ
        #--
        g2 = sns.heatmap(grids[1], ax=ax[1], vmin = 0, cbar = False, cmap="Blues", square=True, xticklabels=int(grids[1].shape[0]/5), yticklabels=int(grids[1].shape[0]/5))
        g2.set_xlabel("x [cm]")
        g2.set_ylabel("z [cm]")
        g2.set_yticklabels(g2.get_yticklabels(), rotation = 0, fontsize=12)
        g2.set_xticklabels(g2.get_xticklabels(), fontsize=12)

        # Borders
        g2.axhline(y = 0, color='k', linewidth=3)    
        g2.axhline(y = grids[1].shape[1], color='k', linewidth=3)    
        g2.axvline(x = 0, color='k', linewidth=3)    
        g2.axvline(x = grids[1].shape[0], color='k', linewidth=3)   
        #--

        plt.close(1)
        plt.tight_layout()
        #plt.legend(frameon=False)
        plt.show()
        self.print_msg('Showing graph ...', press_button=False)

    # View position marginal histogram
    def pos_marg_hist(self, res, axis=0):
        #
        # Gaussian function
        gaussian = lambda x, mean, std_dev, amp: \
            amp * np.exp(-((x - mean)/std_dev)**2 / 2)

        mean, std_dev = res.mass_centre(axis=[axis], fixed_loop_idx=True)

        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        #plt.tight_layout()
        plt.rcParams.update({
                "figure.figsize": (7,6),\
                "font.size":14,\
                "axes.titlepad":16
            })

        #
        # Set labels
        labels = ['x', 'y', 'z']
        plt.title("Marginal histogram " + labels[axis].upper() + "-position")
        plt.xlabel(labels[axis] + " (cm)")
        plt.ylabel(r"density")

        if axis in [0, 1, 2]:
            style={}
            
            # Plot histogram
            plt.bar(res.pos_hist[axis]["bins"], height=res.pos_hist[axis]["dens"], width=0.1, **style)

            #
            # Plot Gaussian Fit
            #--
            max_dens = np.max(res.pos_hist[axis]["dens"])
            x = res.pos_hist[axis]["bins"]
            y = [gaussian(xi, mean, std_dev, max_dens) for xi in x]

            plt.plot(x, y, label="Gaussian fit", linestyle="--", marker="", color="black")
            #-- 

        #
        # Set plot
        plt.grid(linestyle="--")
        plt.legend(frameon=True)

        #
        # Show
        plt.tight_layout()
        plt.show()

    # View velocity marginal histogram
    def vel_marg_hist(self, res, axis=0):
        # Message
        self.print_msg("Processing data ...", press_button=False)

        #
        # Gaussian function
        gaussian = lambda x, mean, std_dev, amp: \
            amp * np.exp(-((x - mean)/std_dev)**2 / 2)

        mean, std_dev = res.average_velocity(axis=[axis], fixed_loop_idx=True)

        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        #plt.tight_layout()
        plt.rcParams.update({
                "figure.figsize": (7,6),\
                "font.size":14,\
                "axes.titlepad":16
            })

        #
        # Set labels
        labels = ['x', 'y', 'z']
        plt.title("Marginal histogram of velocities (" + labels[axis].upper() + "-axis)")
        plt.xlabel(labels[axis] + " [cm/s]")
        plt.ylabel(r"density")

        if axis in [0, 1, 2]:
            style={}
            
            # Plot histogram
            plt.bar(res.vel_hist[axis]["bins"], height=res.vel_hist[axis]["dens"], width=0.9, **style)

            #
            # Plot Gaussian Fit
            #--
            max_dens = np.max(res.vel_hist[axis]["dens"])
            x = res.vel_hist[axis]["bins"]
            y = [gaussian(xi, mean, std_dev, max_dens) for xi in x]

            plt.plot(x, y, label="Gaussian fit", linestyle="--", marker="", color="black")
            #-- 

        #
        # Set plot
        plt.grid(linestyle="--")
        plt.legend(frameon=True)

        #
        # Show
        plt.tight_layout()
        plt.show()

    # Plot centre of mass
    def centre_of_mass(self, res, axis = [0,1,2], trans_id = 1):
        # Message
        self.print_msg("Processing data ...", press_button=False)

        # With loop
        #--
        if res.loop["length"] > 0:
            # Data
            data = res.centre_of_mass(axis=axis)

            # Set style
            plt.clf()
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15, left=0.17)
            plt.style.use('seaborn-whitegrid')
            #plt.tight_layout()
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })

            # Looping info
            info = res.info.loc[res.loop["var"]]

            # x-axis
            #--
            if info["symbol"] == "T_0":
                x_scale_factor, label = self.__temperature_axis(np.max(data["loop_values"]))

            else:
                x_scale_factor = 1

                if info["unit"]: label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
                else: label = r"$ " + info['symbol'] + r"$"

            plt.xlabel(label)
            #--

            #
            # Set labels
            markers = ['o', '^', 's']
            labels = ['x', 'y', 'z']
            plt.ylabel(r"$r_0\ [mm]$")

            # Plot simulated date
            num_axis = len(axis)

            if num_axis > 1:
                for i in range(num_axis):
                    plt.plot(x_scale_factor * data["loop_values"], 10*data["r0"][:,i], linestyle="--", label=labels[axis[i]], marker=markers[axis[i]])

            else:
                plt.plot(x_scale_factor * data["loop_values"], 10*data["r0"], linestyle="--", label=labels[axis[0]], marker=markers[axis[0]])

            # Set plot
            plt.grid(True, linestyle="--")
            plt.close(1)
            plt.tight_layout()
            plt.legend(frameon=True)
            plt.grid(linestyle="--")
            
            # Show
            self.print_msg('Showing graph ...', press_button=False)
            plt.show()
        #--

        # Without looping
        #--
        else: self.print_msg('Function was not implemented. \nPress any buttom to continue ...', press_button=True)
        #--

    # Plot temperature
    def temperature(self, res, log_scale=0, doppler_temperature=False, method=0):
        #
        # Check looping
        if len(res.loop["var"]) > 0:
            # Data
            temp = res.temperature(method = method)*1e6 # uK

            # Clear stored plots
            plt.clf()

            #
            # Set figure
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15)
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })
            #plt.tight_layout()
            ax = plt.gca()

            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            # Set label
            #--
            # x label
            #--
            if info["symbol"] == "T_0":
                x_scale_factor, label = self.__temperature_axis(np.max(x))

            else:
                x_scale_factor = 1

                if info["unit"]:
                    label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
                else:
                    label = r"$ " + info['symbol'] + r"$"

            plt.xlabel(label)
            #--

            # y label
            y_scale_factor, label = self.__temperature_axis(np.max(temp))
            plt.ylabel(label)
            #--

            # Plot temperature
            plt.plot(x_scale_factor*x, y_scale_factor*temp, label="Simulation", marker='o', linestyle='--')

            # Plot Doppler temperature
            if doppler_temperature:
                plt.plot(x, res.doppler_temperature()*np.ones(len(x)), label="Doppler temperature", linestyle='--', marker='', color='black')

            if log_scale == 0:
                plt.xscale('log')

            elif log_scale == 1:
                plt.yscale('log')

            elif log_scale == 2:
                plt.xscale('log')
                plt.yscale('log')

            # Set plot
            plt.grid(True, linestyle="--", which="both")
            ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
            ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())

            if doppler_temperature:
                plt.legend(frameon=True)

            plt.close(1)
            plt.tight_layout()
            plt.show()

        #
        # Mass centre
        else:
            temp = res.temperature(fixed_loop_idx=True, method = method)*1e6 # uK
            print()
            print("T [uK] = %f" % (temp))
            print("T_{doppler} [uK] = %f" % (res.doppler_temperature()))
            print()

    # Plot trap depth vs detuning
    def trap_depth(self, results, fit_func = "erf"):
        # Variables
        Y = np.zeros(len(results))
        delta = []
        opt = None

        # Message
        self.print_msg('Processing data ...', press_button=False)

        # Get capture values
        i = 0
        for code, name in results.items():
            res = Results(code, self.result_id["group"], self._result_id["subgroup"])

            if len(res.loop["values"]) > 0:
                opt = res.loop["var"]
                delta.append(float(res.beams['main'].delta))

                if opt == "T_0":
                    T_c, fit_params = res.capture_temperature(fit_func = fit_func)
                    Y[i] = T_c

                elif opt == "v_0":
                    vel_c, fit_params = res.capture_velocity(fit_func = fit_func)
                    Y[i] = vel_c

                else:
                    raise ValueError("This visualization option requires looping in T_0 or v_0 variables")

                i += 1

            else:
                raise ValueError("This visualization option requires looping in each results set")

        # Capture velocity
        if opt == "v_0":
            y_label = r"$ v_c\ [cm/s] $"
            y_scale_factor = 1

        # Initial temperature
        if opt == "T_0":          
            y_scale_factor, y_label = self.__temperature_axis(np.max(Y)) 

        # Sort elements
        #--
        for i in range(int(len(delta))):
            for j in range(len(delta) - 1):
                if delta[j] > delta[j + 1]:
                    aux = delta[j]
                    delta[j] = delta[j + 1]
                    delta[j + 1] = aux

                    aux = Y[j]
                    Y[j] = Y[j + 1]
                    Y[j + 1] = aux
        #--

        # Clear stored plots
        plt.clf()

        #
        # Set figure
        plt.figure(figsize=(5,4))
        plt.style.use('seaborn-whitegrid')
        plt.subplots_adjust(top=0.80, bottom=0.15)
        plt.rcParams.update({
                "font.size":14,\
                "axes.titlepad":14
            })
        #plt.tight_layout()
        ax = plt.gca()

        # Looping info
        info = res.info.loc["delta"]

        # Set label
        #--            
        if info["unit"]:
            x_label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
        else:
            x_label = r"$ " + info['symbol'] + r"$"

        # x label
        plt.ylabel(y_label)
        plt.xlabel(x_label)
        #--

        # Plot trapped atoms ratio
        plt.plot(delta, y_scale_factor * Y, label="Simulated Points", marker='o', linestyle='--')

        # Set plot
        plt.grid(True, linestyle="--")

        plt.close(1)
        plt.tight_layout()
        
        #
        # Show
        print('Showing graph ...')
        plt.show()

    # Plot escape flux of atoms
    def escape_flux_atoms(self, res):
        self.print_msg("Processing data ...", press_button=False)

        # Plot graph
        if res.loop["length"] > 0:
            # Data
            X, Y = res.escape_flux_atoms()

            # Set figure
            plt.clf() # Clear previously plots
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })

            # Looping info
            info = res.info.loc[res.loop["var"]]

            # Set labels
            #--
            # x label        
            if info["symbol"] == "T_0":
                x_scale_factor, x_label = self.__temperature_axis(np.max(X))
            
            else:
                x_scale_factor = 1
                x_label = r"$ " + info['symbol'] + r"$"           
                if info["unit"]: x_label += r"$\ [" + info['unit'] + r"] $"

            plt.xlabel(x_label)

            # y label
            plt.ylabel(r"$ \Phi\ [atoms\ /\ s]$")
            #--

            # Plot trapped atoms ratio
            plt.errorbar(x_scale_factor*X, Y, label="MOT On", marker='o', linestyle='')

            # Set plot
            plt.grid(True, linestyle="--")
            plt.close(1)
            plt.tight_layout()
            
            # Show
            self.print_msg('Showing graph ...', press_button=False)
            plt.show()

        # View data
        else:
            data = res.escape_flux_atoms()
            self.print_msg("Escape flux of atoms: %.2f [atoms / s]" % data)

    # Plot trapped atoms ratio
    def trapped_atoms(self, res, fit_func = None):
        # Message
        self.print_msg("Processing data ...", press_button=False)

        # With loop
        #--
        if res.loop["length"] > 0:
            # Data
            X, Y = res.normalized_trapped_atoms()

            # Labels
            #--
            info = res.info.loc[res.loop["var"]]
            if info["symbol"] == "T_0":
                x_scale_factor, x_label = self.__temperature_axis(np.max(X))

            else:
                x_scale_factor = 1
                x_label = r"$ " + info['symbol']
                if info["unit"]: x_label += r"\ [" + info['unit'] + r"] $"
            #--
                

            # Set figure
            plt.clf() # Clear plots
            plt.figure(figsize=(6,5)) # Set figure size
            plt.style.use('seaborn-whitegrid') # Theme
            #plt.subplots_adjust(top=0.80, bottom=0.15)
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })
            ax = plt.gca()

            # Axis labels
            plt.xlabel(x_label)
            plt.ylabel(r"$N / N_{total}$")
            #--

            # Sort elements
            #--
            for i in range(int(len(X))):
                for j in range(len(X) - 1):
                    if X[j] > X[j + 1]:
                        aux = X[j]
                        X[j] = X[j + 1]
                        X[j + 1] = aux

                        aux = Y[j]
                        Y[j] = Y[j + 1]
                        Y[j + 1] = aux
            #--

            # Plot trapped atoms ratio vs looping variable
            plt.plot(x_scale_factor * X, Y, label="Simulated points", marker='o', linestyle="--")

            # Plot trapped atoms ratio vs escape velocity
            #if res.escape_vel >= 0 and info["symbol"] == "v_0":
            #    X2, Y2, Z = np.array(res.all_escape_velocities())
            #    plt.plot(x_scale_factor * Y2, Z, label="Escape velocity", marker='o', linestyle="", color="red")
            #    plt.xlabel(r"$v\ [cm / s]$")

            # Fitting
            #--
            # Initial Velocity
            fit_opts = ["poly", "erf"]

            # Fittings
            if (fit_func in fit_opts) and (res.loop["var"] in ["v_0", "T_0"]):
                if res.loop["var"] == "v_0":
                    vel_c, fit_params = res.capture_velocity(fit_func=fit_func)
                    label_c = r"$v_{cap} = %.2f\ cm/s $" % (x_scale_factor * vel_c)
                    if vel_c > 0: plt.plot(x_scale_factor * vel_c, res.cut_trap_depth, marker="*", markersize=10, linestyle="", label=label_c)

                elif res.loop["var"] == "T_0":
                    T_c, fit_params = res.capture_temperature(fit_func=fit_func)
                    label_c = x_label + (r"$ = %.1f$" % (x_scale_factor * T_c))
                    if T_c > 0: plt.plot(x_scale_factor * T_c, res.cut_trap_depth, marker="*", markersize=10, linestyle="", label=label_c)
                
                # Polynomial function
                if fit_func == "poly":
                    def f(x):
                        y = 0.0

                        for i in range(11):
                            y += fit_params[10 - i]*x**i

                        return y

                elif fit_func == "erf":
                    # Error function
                    def f(t):
                        return fit_params[2] * (1 - (erf((t - fit_params[0]) / np.sqrt(2 * fit_params[1]**2)) - erf((- fit_params[0]) / np.sqrt(2 * fit_params[1]**2))) / 2)

                x_fit = np.linspace(np.min(X), np.max(X), 1000)
                y_fit = np.array(list(map(f, x_fit)))
                
                plt.plot(x_scale_factor * x_fit, y_fit, label="Fitting", marker="", linestyle="-", linewidth=2, color="Black")            
                #--

            # Set plot
            plt.grid(True, linestyle="--")
            plt.legend(frameon = True)

            plt.close(1)
            plt.tight_layout()
            
            # Show
            self.print_msg('Showing graph ...', press_button=False)
            plt.show()
        #--

        # Without loop
        #--
        else:
            self.print_msg("Normalized trapped atoms = %.2f" % res.normalized_trapped_atoms())
        #--

    # Plot escape velocity
    def escape_vel(self, res):
        # With looping
        #--
        if res.loop["values"]:
            # Get X-values
            X = np.array(res.loop["values"])

            # Get Y-values
            #--
            Y = np.zeros(X.size)
            for i in range(Y.size):
                res.loop_idx(i)
                Y[i] = res.escape_vel
            #--

            #
            # Set figure
            plt.clf() # Clear previously plots
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })

            # Looping info
            info = res.info.loc[res.loop["var"]]

            # Set labels
            #--
            # x label        
            if info["symbol"] == "T_0":
                x_scale_factor, x_label = self.__temperature_axis(np.max(X))
            else:
                x_scale_factor = 1
                x_label = r"$ " + info['symbol'] + r"$"           
                if info["unit"]: x_label += r"$\ [" + info['unit'] + r"] $"

            plt.xlabel(x_label)

            # y label
            plt.ylabel(r"$ v_{e}\ [cm / s]$")
            #--

            # Plot escape velocity
            plt.plot(x_scale_factor*X, Y, label="Simulated points", marker='o', linestyle='')

            # Set plot
            plt.grid(True, linestyle="--")
            plt.close(1)
            plt.tight_layout()
            
            # Show
            print('Showing graph ...')
            plt.show()
        #--

        # Without looping
        else:
            print("Escape velocity = %f" % res.escape_vel)

    # Plot r.m.s. cloud size
    def cloud_size(self, res):
        if res.loop["var"]:
            #
            # Get data
            r_c, std_r_c = res.mass_centre()

            #
            # Clear stored plots
            plt.clf()

            # Set figure
            plt.figure(figsize=(5,4))
            plt.style.use('seaborn-whitegrid')
            plt.subplots_adjust(top=0.80, bottom=0.15, left=0.17)
            plt.style.use('seaborn-whitegrid')
            #plt.tight_layout()
            plt.rcParams.update({
                    "font.size":14,\
                    "axes.titlepad":14
                })

            # Looping info
            info = res.info.loc[res.loop["var"]]
            x = np.array(res.loop["values"]).astype(float)

            #
            # Set labels
            #--
            markers = ['o', '^', 's']
            labels = [r'$\sigma_x$', r'$\sigma_y$', r'$\sigma_z$']

            # x label
            #--
            if info["symbol"] == "T_0":
                x_scale_factor, label = self.__temperature_axis(np.max(x))

            else:
                x_scale_factor = 1

                if info["unit"]:
                    label = r"$ " + info['symbol'] + r"\ [" + info['unit'] + r"] $"
                else:
                    label = r"$ " + info['symbol'] + r"$"

            plt.xlabel(label)
            #--

            plt.ylabel("Size [mm]")
            #--

            # Plot simulated date
            for i in range(3):
                plt.plot(x_scale_factor * x, std_r_c[i]*10, label=labels[i], linestyle="--", marker=markers[i])

            # Set plot
            plt.grid(linestyle="--")
            plt.legend(frameon=True)

            # Show
            plt.close(1)
            plt.tight_layout()
            plt.show()

        else:
            print('Visualization not implemented')
    
    # Heat map
    def heatmap(self, res, axis, val):
        #
        # Get data
        hist = res.pos_2Dhist(axis, val)*1e3;

        #
        # Clear stored plots
        plt.clf()

        #
        # Set style
        plt.style.use('seaborn-whitegrid')
        #plt.tight_layout()
        plt.rcParams.update({
                "figure.figsize": (7,6),\
                "font.size":12,\
                "axes.titlepad":12
            })

        #
        # Set labels
        plane_label = ['yz', 'xz', 'xy']
        axis_label = ['x', 'y', 'z']

        plt.title("Position distribution in " + plane_label[axis] + "-axis and " + axis_label[axis] + " = " + ("%.2f" % float(val)))
        plt.xlabel(plane_label[axis][0] + " [cm]")
        plt.ylabel(plane_label[axis][1] + " [cm]")

        #
        # Plot simulated date
        l = float(res.perform['max_r'])
        plt.imshow(hist, cmap="Blues", vmin=np.min(hist), vmax=np.max(hist), extent=[-l, l, -l, l])

        #
        # Show
        plt.tight_layout()
        plt.show()

    # Temperature axis
    def __temperature_axis(self, T_max, symbol="T"):
        if T_max >= 1e3:
            scale_factor = 1e-3
            label = r"$" + symbol + r"\ [K]$"

        elif T_max <= 1:
            scale_factor = 1e-3
            label = r"$" + symbol + r"\ [\mu K]$"

        else:
            scale_factor = 1
            label = r"$" + symbol + r"\ [m K]$"

        return scale_factor, label

    #
    def available_views(self, result):
        # Available options
        available_opts = []

        # All individual options
        all_opts = {
            1: "Histogram of positions",\
            2: "Heatmaps of positions",\
            3: "Centre of mass",\
            4: "Cloud size",\
            5: "Histogram of velocities",\
            6: "Temperature",\
            7: "Escape flux of atoms",\
            8: "Trapped atoms",\
            9: "Trapped atoms (escape velocity)",\
            10: "Escape velocity",\
            11: "Initial velocity distribution"
        }

        # Available options
        #--
        # Position visualizations
        if len(result.pos_hist[0]["freqs"]) > 0:
            available_opts.append(1)
            available_opts.append(2)
            available_opts.append(3)
            available_opts.append(4)

        # Velocity visualizations
        if len(result.vel_hist[0]["freqs"]) > 0:
            available_opts.append(5)
            available_opts.append(6)

        # Trapped atoms
        if result.trapped_atoms >= 0:            
            available_opts.append(7)
            available_opts.append(8)

        # Escape velocity
        if result.escape_vel >= 0:
            #available_opts.append(9)
            available_opts.append(10)

        # Initial velocity distribution
        if result.ini_vel.size > 0:
            available_opts.append(11)
        #--

        return available_opts, all_opts

    #-- 
    # Interface methods
    #--

    # Print the results
    def results_history(self, num = 5, clear_screen = True):
        #
        # Get results
        res = self.__simulation.available_results(num)

        #
        # Clear screen
        if os.name == 'nt': os.system('cls')
        else: os.system('clear')

        #
        # Show header 
        print(self.header)

        #
        # Show results
        print('Results history' + self.separator)

        for i, val in enumerate(res):
            print(str(i + 1), end='')
            print(" - (" + str(val[0]) + ")", end='')
            print(' ' + str(dt.fromtimestamp(val[0])), end='')
            print(' ' + val[1], end='')
            print()
    
    # Print the results groups
    def results_groups(self):
        #
        # Get results
        res = self.__simulation.available_results_groups()

        #
        # Clear screen
        if os.name == 'nt': os.system('cls')
        else: os.system('clear')

        #
        # Show header 
        print(self.header)

        #
        # Show groups
        print('Results groups' + self.separator)

        for i, val in enumerate(res):
            print(str(i + 1), end='')
            print(" - (" + str(val[0]) + ")", end='')
            print(' ' + str(dt.fromtimestamp(val[0])), end='')
            print(' ' + val[1], end='')
            print()
    
    #
    def terminal_menu(self, options, desc = None, msg = None, clear_screen = True, show_main_header = True, enumerated_list = False):
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        # Main header
        if show_main_header: print(self.header, end="")

        # Result identification
        result_id_view = self.__result_id_view()
        if result_id_view is not None: print(result_id_view, end="\n")
        print()

        # Description
        if desc is not None: print(desc + self.separator)

        # Options
        i = 0
        opts_size = len(options)
        for key, val in options.items():
            i += 1
            if enumerated_list: print("[%d] " % i, end='')

            # Print option
            if val: print('%d -> %s;' % (key, val), end='')
            else: print('%d;' % (key), end='')
            if i < opts_size: print()

        print(self.separator, end="\n\n")

        # Message
        if msg is not None: print(msg)

        # Return input
        return input('Option code: ')

    #
    def terminal_input(self, desc = None, msg = None, clear_screen = True, show_header = True):
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

        # Show header 
        if show_header: print(self.header, end="")

        # Result identification
        result_id_view = self.__result_id_view()
        if result_id_view is not None: print(result_id_view)
        print()

        # Message
        if msg is not None: print(msg)

        # Input
        input_str = input(desc + ": ")

        return input_str

    # Print general information of a simulation
    def simulation_status(self, sim_opt = None):
        # Clear screen
        if os.name == 'nt': os.system('cls')
        else: os.system('clear')

        # Header 
        print(self.header)
        if sim_opt is not None: 
            print(sim_opt + self.separator)

        # Group
        if len(self.result_id["group"]) > 0:
            group_view = "(Group) " + self.result_id["group"]
            
            if (self.result_id["subgroup"] is not None) and (self.result_id["subgroup"] != "root"): 
                group_view += " / (Subgroup) " + self.result_id["subgroup"]
            
            group_view += " / " + str(self.__simulation.results.id["code"])
            
            if self.__simulation.results.id["name"] is not None:
                group_view += " " + self.__simulation.results.id["name"]

            print(group_view + self.separator, end="\n")

        # Loop information
        if self.__simulation.results.loop["length"] > 0:
            loop = self.__simulation.results.loop
            info = self.__simulation.results.info
            print('Number of loops: %d' % loop["length"])
            print("Loop variable: (%s) %s [%s]" % (loop["var"], info.loc[loop["var"], "name"], info.loc[loop["var"], "unit"]))

        print()

    #
    def print_msg(self, msg, press_button=True, clear_screen=True):
        # Clear screen
        if clear_screen:
            if os.name == 'nt': os.system('cls')
            else: os.system('clear')

            # Show headers 
            print(self.header, end="")

            # Result identification
            result_id_view = self.__result_id_view()
            if result_id_view is not None: print(result_id_view, end="\n\n")

        print(msg)

        if press_button:
            input("Press any key to continue ...")

    #
    def __result_id_view(self):
        result_id_view = None

        if self.result_id["group"] is not None:
            result_id_view = "(Group) " + self.result_id["group"]
            
            if self.result_id["subgroup"] is not None: 
                result_id_view += " / (Subgroup) " + self.result_id["subgroup"]

            if self.result_id["code"] is not None:
                result_id_view += "\n -- " + str(self.result_id["code"])

                if self.result_id["name"] is not None:
                    result_id_view += " " + self.result_id["name"]

                result_id_view += " --"
            
        
        return result_id_view