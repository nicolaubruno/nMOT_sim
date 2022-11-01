# Libraries and modules
#--
import numpy as np
import time, gc
import pandas as pd

from model import Simulation, Results
from view import View

from datetime import datetime as dt
from tqdm import tqdm
#--

# Controller class
class Controller:
    
    #--
    # Attributes
    #--
    
    #
    @property
    def result_id(self):
        return self._result_id
    
    #--
    # Operational Methods 
    #--

    #
    def __init__(self):
        self.__simulation = Simulation()
        self.__view = View(self.__simulation)

        # Result identification
        self._result_id = {
            "code": None,\
            "name": None,\
            "group": None,\
            "subgroup": None
        }

    #
    def get_group(self):
        group_data = pd.read_csv("controller/active_group.csv", index_col=0, header=None, squeeze=True)
        group_data.fillna(value={"subgroup":None}, inplace=True)
        self.result_id["group"] = group_data["group"]
        self.result_id["subgroup"] = group_data["subgroup"]
        
        self.__view.set_result_id(self.result_id)

    #--
    # Methods to manage the interface
    #--

    #
    def main_menu(self):
        # Get group and subgroup
        self.get_group()

        # Initial menu level
        menu_level = 1

        # Message
        msg = None

        # Loop menu
        while True:
            # Options
            if menu_level == 1:
                # Options
                opts = {
                    1: "Run simulation",\
                    2: "Get results data",\
                    3: "Change group/subgroup"
                }

                # Description
                desc = "Choose an option:"

                # Get group
                chosen_opt, back_menu, msg = self.__call_menu(opts, desc, msg=msg)

                # Simulate a new dataset
                if chosen_opt == 1: self.simulation_menu()

                # Get results data
                elif chosen_opt == 2: self.results_menu()

                # Change group/subgroup
                elif chosen_opt == 3: self.group_menu()

                # Change menu level
                menu_level += -2*back_menu

            # Leave menu
            else:
                print("Exiting ...")
                exit(0)

    # Set the parameters and run a new simulation
    def simulation_menu(self):
        # Message
        msg = None

        # Menu level
        menu_level = 1

        # Loop menu
        while True:
            # Result short name
            if menu_level == 1:
                # Description
                desc = "Run simulation / Short Name"

                # Get shortname
                shortname, back_menu = self.__call_input(desc)
                menu_level += -2*back_menu + 1

            # Execute simulation
            elif menu_level == 2:
                # Create a new simulation
                self.__simulation.new(self.result_id["group"], self.result_id["subgroup"], shortname)

                # Print simulation status
                self.__view.simulation_status()
                
                # Check looping
                loop_num = self.__simulation.results.loop["length"]
                if loop_num <= 0: loop_num = 1

                # Set progress bars
                #---
                pbars = []
                desc = None

                for i in range(loop_num):
                    if loop_num > 1:
                        desc = str(self.__simulation.results.loop["values"][i])

                    pbars.append(tqdm(total=self.__simulation.results.perform["num_sim"], desc=desc, position=i))
                #---


                # Run simulation
                #--
                for i in range(loop_num):
                    # Change loop
                    self.__simulation.change_loop(i)
                    
                    # Simulate atoms
                    #--
                    while self.__simulation.atoms_simulated < self.__simulation.results.perform["num_sim"]:
                        # Simulate atoms
                        times = self.__simulation.run()
                        
                        pbars[i].update(times)
                    #--

                    # Save simulation
                    self.__simulation.save()
                #--

                # Close progress bars
                for i in range(loop_num): pbars[i].close()

                #
                # Release memory
                gc.collect()

                # Information about the simulation
                self.__view.print_msg("\nThe simulation has finished!", clear_screen=True)
                menu_level = 0

            # Leave menu
            else: break

            '''

            # Create a new simulation
            self.__simulation.new(shortname, opt, results_group)

            # Check looping
            loop_num = len(self.__simulation.results.loop["values"])
            if loop_num == 0: loop_num = 1

            # Update time
            check_time = dt.now().timestamp()

            # Print simulation status
            self.__view.simulation_header(group=available_groups[results_group], sim_opt=opts[opt], clear_screen=True)

            # Set progress bars
            pbars = []
            for i in range(loop_num):
                desc = "Atoms simulated" if loop_num == 1 else self.__simulation.results.loop["var"] + " = " + ("%.2f" % self.__simulation.results.loop["values"][i])
                pbars.append(tqdm(total=self.__simulation.results.perform["num_sim"], desc=desc, position=i))

            #
            # Run simulation
            #--
            for i in range(loop_num):
                # Open new simulation for each looping value
                if i > 0: self.__simulation.open(self.__simulation.results.code, i, opt)
                
                #
                # Simulate atoms
                #--
                while self.__simulation.atoms_simulated < self.__simulation.results.perform["num_sim"]:
                    # Simulate atoms
                    times = self.__simulation.run()
                    
                    pbars[i].update(times)
                #--

                # Save simulation
                self.__simulation.save()
            #--

            # Close progress bars
            for i in range(loop_num): pbars[i].close()

            #
            # Release memory
            gc.collect()

            # Information about the simulation
            #self.__view.simulation_header(group=available_groups[results_group], sim_opt=opts[opt], clear_screen=True)
            input("Enter any key to continue ... ")

            #
            # Set menu level
            self._menu_level = 0
            '''

    # Plot/Get results
    def results_menu(self):
        # Message
        msg = None

        # Set menu level
        menu_level = 1

        # Loop menu
        while True:
            # Check collective or individual results option
            #--
            if menu_level == 1:
                # Description
                desc = "Plot/Get Results / Choose an option:"

                # Options
                opts = {
                    1: "Individual results",\
                    2: "Collective results"
                }

                # Get group
                collective_results, back_menu, msg = self.__call_menu(opts, desc, msg)
                collective_results -= 1

                # Set menu level
                menu_level += -2*back_menu + 1
            #--

            # Results option
            elif menu_level > 1:
                # Collective results
                if collective_results:
                    # Get result option
                    if menu_level == 2:
                        # Description
                        desc = "Collective Results / Choose an option:"

                        # Options
                        #--
                        opts = {
                            1: "Trap depth"
                        }
                        #--

                        # Get group
                        view_opt, back_menu, msg = self.__call_menu(opts, desc, msg)
                        menu_level += -2*back_menu + 1

                    # Graph or CSV file option
                    elif menu_level == 3:
                        # Description
                        #--
                        desc = "... / " + opts[view_opt]
                        desc += " / Choose an option"
                        #--

                        # Options
                        #--
                        opts = {
                            1: "Visualization",\
                            2: "Get CSV file"
                        }
                        #--

                        # Get group
                        view_or_CSV, back_menu, msg = self.__call_menu(opts, desc, msg)

                        # Set menu level
                        menu_level += -2*back_menu + 1

                    # Visualization
                    elif menu_level == 4:
                        # Trap depth
                        if view_opt == 1:
                            if view_or_CSV == 1:
                                results = self.__simulation.available_results(self.result_id)
                                self.__view.trap_depth(results)
                                menu_level -= 1

                            else:
                                self.__view.print_msg("Creating CSV file ...", press_button = False)

                                results = self.__simulation.available_results(self.result_id)
                                X = np.zeros(len(results))
                                Y = np.zeros(len(results))
                                
                                i = 0
                                for code, name in results.items():
                                    res = Results(code, self.result_id["group"], self.result_id["subgroup"])

                                    X[i] = res.beams["main"].delta

                                    if res.loop["var"] == "T_0":
                                        T_c, fit_params = res.capture_temperature()
                                        Y[i] = T_c if T_c > 0 else 0

                                    elif res.loop["var"] == "v_0":
                                        vel_c, fit_params = res.capture_velocity()
                                        Y[i] = vel_c if vel_c > 0 else 0

                                    i += 1

                                # Save CSV file
                                #--
                                data = {}
                                data["delta"] = X

                                if res.loop["var"] == "T_0 [mK]": data["T_c"] = Y
                                else: data["v_c [cm/s]"] = Y

                                path = "temp/trap_depth.csv"
                                pd.DataFrame(data).to_csv(path, index=False)
                                #--

                                self.__view.print_msg("The CSV file has been created! Path: %s" % path)

                                # Set menu level
                                menu_level -= 2

                # Individual results
                else:
                    # Get code of a simulated data
                    if menu_level == 2:
                        # Description
                        desc = "Individual Results / Choose a set of simulated data:"

                        # Available results
                        #--
                        results_opts = self.__simulation.available_results(self.result_id)
                        #--

                        # Get result object
                        results_code, back_menu, msg = self.__call_menu(results_opts, desc, msg, enumerated_list=True)
                        
                        if results_code > 0:
                            self.__view.print_msg("Loading result ...", press_button=False)
                            results = Results(results_code, self.result_id["group"], self.result_id["subgroup"])
                            self._result_id['code'] = results.id["code"]
                            self.result_id["name"] = results.id["name"]
                            self.__view.set_result_id(self.result_id)

                        # Set menu level
                        menu_level += -2*back_menu + 1

                    # View option lvl 1
                    elif menu_level == 3:
                        # Description
                        #--
                        desc = "... / Individual results / Choose an option"
                        #--

                        # Options
                        available_view_opts, all_view_opts = self.__view.available_views(results)
                        opts = {(k+1):all_view_opts[v] for (k,v) in enumerate(available_view_opts)}

                        # Get visualization
                        view_opt, back_menu, msg = self.__call_menu(opts, desc, msg)
                        view_opt = available_view_opts[view_opt-1]
                        menu_level += -2*back_menu + 1

                    # View option lvl 2 (View of CSV File)
                    elif menu_level == 4:
                        # Set view as default
                        view_or_CSV = 1

                        # Check loop
                        if results.loop["length"] > 0:
                            # Description
                            #--
                            desc = "... / " + all_view_opts[view_opt]
                            desc += " / Choose an option"
                            #--

                            # Options
                            #--
                            opts = {
                                1: "Visualization",\
                                2: "Get CSV file"
                            }
                            #--

                            # Get group
                            view_or_CSV, back_menu, msg = self.__call_menu(opts, desc, msg)

                            # Set menu level
                            menu_level += -2*back_menu + 1
                        
                        # Visualization
                        else: menu_level += 1

                    # View option lvl 3
                    elif menu_level == 5:
                        # Histogram of positions (or velocities)
                        if view_opt == 1 or view_opt == 5:
                            # Description
                            desc = "... / " + all_view_opts[view_opt] + "/ Visualization / "
                            desc += "\nChoose a loop value (" + results.info.loc[results.loop["var"], "name"] + "):"

                            # Options
                            keys = np.arange(1, results.loop["length"] + 1, 1)
                            values = results.loop["values"]
                            opts = {key:values[key - 1] for key in keys}
                            loop_idx, back_menu, msg = self.__call_menu(opts, desc, msg)
                            loop_idx -= 1
                            menu_level += -2*back_menu + 1

                        else: menu_level = 6

                    # View data
                    elif menu_level == 6:
                        # Histogram of positions
                        if view_opt == 1:
                            # Set loop value
                            results.loop_idx(loop_idx)

                            # Plot graph
                            if view_or_CSV == 1:
                                self.__view.pos_histograms(results)
                                menu_level -= 2

                            # Get CSV file
                            elif view_or_CSV == 2:
                                self.__view.print_msg("Option has not implemented!", press_button = False)

                                # Set menu level
                                menu_level -= 2

                        # Position Heatmap
                        elif view_opt == 2:
                            # Plot graph
                            if view_or_CSV == 1:
                                self.__view.pos_heatmap(results)
                                menu_level -= 2

                            # Get CSV file
                            elif view_or_CSV == 2:
                                self.__view.print_msg("Option has not implemented!", press_button = False)

                                # Set menu level
                                menu_level -= 2

                        # Centre of mass
                        elif view_opt == 3:
                            # Plot graph
                            if view_or_CSV == 1:
                                self.__view.centre_of_mass(results)
                                menu_level -= 2

                            # Get CSV file
                            elif view_or_CSV == 2:
                                self.__view.print_msg("Creating CSV file ...", press_button = False)
                                data = results.centre_of_mass(CSV_file = True)

                                # Save CSV file
                                #--
                                path = "temp/"
                                path += str(self.result_id["code"])
                                path += "_centre_of_mass.csv"
                                pd.DataFrame(data).to_csv(path, index=False)
                                #--

                                self.__view.print_msg("The CSV file has been created! Path: %s" % path)

                                # Set menu level
                                menu_level -= 2

                        # Cloud size
                        elif view_opt == 4:
                            self.__view.print_msg("Function has not implemented!")
                            menu_level -= 2
                                
                        # Histogram of velocities
                        elif view_opt == 5:
                            # Set loop value
                            results.loop_idx(loop_idx)

                            # Plot graph
                            if view_or_CSV == 1:
                                self.__view.vel_histograms(results)
                                menu_level -= 2

                            # Get CSV file
                            elif view_or_CSV == 2:
                                self.__view.print_msg("Option has not implemented!", press_button = False)

                                # Set menu level
                                menu_level -= 2

                        # Temperature
                        elif view_opt == 6:
                            # Plot graph
                            if view_or_CSV == 1:
                                self.__view.temperature(results)
                                menu_level -= 2

                            # Get CSV file
                            elif view_or_CSV == 2:
                                self.__view.print_msg("Creating CSV file ...", press_button = False)
                                data = results.temperature(CSV_file = True)

                                # Save CSV file
                                #--
                                path = "temp/"
                                path += str(self.result_id["code"])
                                path += "_temperature.csv"
                                pd.DataFrame(data).to_csv(path, index=False)
                                #--

                                self.__view.print_msg("The CSV file has been created! Path: %s" % path)

                                # Set menu level
                                menu_level -= 2

                        # Escape flux of atoms
                        elif view_opt == 7:
                            # Plot graph
                            if view_or_CSV == 1:
                                self.__view.escape_flux_atoms(results)
                                menu_level -= 2

                            # Get CSV file
                            elif view_or_CSV == 2:
                                self.__view.print_msg("Creating CSV file ...", press_button = False)
                                X, Y = results.escape_flux_atoms()

                                # Save CSV file
                                #--
                                data = {
                                    results.loop["var"]: X,\
                                    "escape_flux_atoms": Y
                                }

                                path = "temp/"
                                path += str(self.result_id["code"])
                                path += "_escape_flux_atoms.csv"
                                pd.DataFrame(data).to_csv(path, index=False)
                                #--

                                self.__view.print_msg("The CSV file has been created! Path: %s" % path)

                                # Set menu level
                                menu_level -= 2

                        # Trapped atoms
                        elif view_opt == 8:
                            # Plot graph
                            if view_or_CSV == 1:
                                self.__view.trapped_atoms(results, fit_func="erf")
                                menu_level -= 2

                            # Get CSV file
                            elif view_or_CSV == 2:
                                self.__view.print_msg("Creating CSV file ...", press_button = False)
                                X, Y = results.normalized_trapped_atoms()

                                # Save CSV file
                                #--
                                data = {
                                    results.loop["var"]: X,\
                                    "trapped_atoms": Y
                                }

                                path = "temp/"
                                path += str(self.result_id["code"])
                                path += "_trapped_atoms.csv"
                                pd.DataFrame(data).to_csv(path, index=False)
                                #--

                                self.__view.print_msg("The CSV file has been created! Path: %s" % path)

                                # Set menu level
                                menu_level -= 2

                        # Trapped atoms (Escape velocity)
                        elif view_opt == 9:
                            self.__view.print_msg("Function has not implemented!")
                            menu_level -= 2

                        # Escape velocity
                        elif view_opt == 10:
                            self.__view.print_msg("Function has not implemented!")
                            menu_level -= 2
                        
                        # Initial velocity distribution
                        elif view_opt == 11:
                            if results.loop["length"] > 1:
                                # Description
                                #--
                                desc = "... / " + all_view_opts[view_opt]
                                desc += " / Choose an option:"
                                #--

                                # Options
                                #--
                                opts = {}
                                for i in range(results.loop["length"]):
                                    opts[i+1] = results.loop["var"] + " = " + str(results.loop["values"][i])
                                #--

                                # Get loop index
                                loop_idx, back_menu, msg = self.__call_menu(opts, desc, msg)

                                # Set menu level
                                menu_level += -2*back_menu + 1
                                if back_menu == 0: results.loop_idx(loop_idx-1)

                            else: menu_level += 1

                        # Option has not implemented
                        else:
                            self.__view.print_msg("Option has not implemented!", press_button = False)
                            menu_level -= 1

                    # View data lvl 2
                    elif menu_level == 6:
                        # Initial velocity distribution
                        if view_opt == 11:
                            # Plot graph
                            if view_or_CSV == 1:
                                self.__view.ini_vel_hist(results)
                                menu_level -= 3

                            # Get CSV file
                            elif view_or_CSV == 2:
                                self.__view.print_msg("Option has not implemented!", press_button = False)

                                # Set menu level
                                menu_level -= 3
                        
                        # Option has not implemented
                        else:
                            self.__view.print_msg("Option has not implemented!", press_button = False)
                            menu_level -= 3

                    # Leave menu
                    else: 
                        print("Exiting ...")
                        exit(0)

            # Leave menu
            else: break

    # Set the active group and subgroup
    def group_menu(self):
        # Message
        msg = None

        # Set menu level
        menu_level = 1

        # Loop menu
        while True:
            # Get group
            if menu_level == 1:
                # Available groups
                opts = self.__simulation.available_groups()
                desc = "Choose a group:"

                # Get group
                chosen_opt, back_menu, msg = self.__call_menu(opts, desc, msg)
                if back_menu == 0: self._result_id["group"] = opts[chosen_opt]

                # Change menu level
                menu_level += -2*back_menu + 1

            # Get subgroup:
            elif menu_level == 2:
                # Available subgroups
                opts = self.__simulation.available_subgroups(self.result_id["group"])

                # Check available subgroups
                if len(opts) > 0:
                    # Get subgroup
                    desc = ("(Group) %s / Choose a subgroup:" % self.result_id["group"])
                    chosen_opt, back_menu, msg = self.__call_menu(opts, desc, msg)
                    if back_menu == 0: self._result_id["subgroup"] = opts[chosen_opt]

                    # Change menu level
                    menu_level += -2*back_menu + 1

                # Leave looping
                else: 
                    menu_level = 0
                    self._result_id["subgroup"] = None

            # Leave looping
            else: break

        # Change CSV file
        #--
        group_data = {
            "group":self.result_id["group"],\
            "subgroup":self.result_id["subgroup"]
        }

        pd.Series(group_data, name="active_group").to_csv("controller/active_group.csv")
        #--

        # Change view
        self.__view.set_result_id(self.result_id)

    #
    def __call_menu(self, opts, desc = None, msg = None, clear_screen = True, enumerated_list = False):
        # Variables
        back_menu = 0

        # Call menu
        #--
        while True:
            # Get option
            opt = self.__view.terminal_menu(opts, desc=desc, msg=msg, clear_screen=clear_screen, enumerated_list = enumerated_list)
            msg = None

            # Exit
            if opt == '0':
                print('\nExiting ...', end='\n\n')
                exit(0)

            # Back
            elif opt == "-1":
                opt = int(opt)
                back_menu = 1
                break

            # Get in order
            elif opt[:2] == "_o" and len(opt) > 2 and opt[2:].isdigit():
                opts_keys = list(opts.keys())
                if ((int(opt[2:]) - 1) < len(opts_keys)) and int(opt[2:]) > 0:
                    opt = opts_keys[int(opt[2:]) - 1]
                    break

                else: msg = "Option not found"

            # Check value
            elif len(opt) > 0 and (opt.isdigit() or (opt[0] == "-" and opt[1:].isdigit())):
                opt = int(opt)

                if (opt in opts.keys()): break
                else: msg = "Invalid option"

            else: msg = "Invalid option"
        #--

        return opt, back_menu, msg

    #
    def __call_input(self, desc, validation=None, back_opt=True):
        # Variables
        back_menu = 0
        msg = None

        # Get value
        opt = self.__view.terminal_input(desc)

        # Exit
        if opt == '0':
            print('\nExiting ...\n')
            exit(0)

        # Back
        elif opt == "-1" and back_opt:
            back_menu = 1

        # Validation
        elif (validation is not None) and (opt != "-1"):
            while (not validation(opt)) and (opt != "-1"):
                msg = "Invalid value!"
                opt = self.__view.terminal_input(description, header=header, clear_screen=clear_screen, footer=msg)

                # Exit
                if opt == '0':
                    print('\nExiting ...\n')
                    exit(0)

                # Back
                elif opt == "-1":
                    opt = int(opt)
                    back_menu = 1
                    break

                # Empty
                elif len(opt.strip()) == 0:
                    opt = None

        return opt, back_menu
