# Library
import numpy as np
import pandas as pd

class BeamsSetup:
    #--
    # Getters
    #--

    # Number of beams
    @property
    def num_beams(self) -> int:
        return self._num_beams

    # General values
    @property
    def general(self):
        return self._general
    
    # All Beams
    @property
    def beams(self) -> dict:
        return self._beams
    
    #--
    # Operational methods
    #--

    #
    def __init__(self, dir_path: str):
        # General values
        #--
        self._general = pd.read_csv(dir_path + "general.csv", header=0, index_col=0, squeeze=True).astype(object)

        # Casting sidebands
        self._general["sidebands"] = self.general["sidebands"][1:-1].split(" ")
        self._general["sidebands"] = {"num": int(self.general["sidebands"][0]), "freq": float(self.general["sidebands"][1])}
        #--

        # Get all beams
        #--
        beams_setup = pd.read_csv(dir_path + "setup.csv", header=0, index_col=0).astype(object)
        self._num_beams = beams_setup.shape[0]
        self._beams = {}

        # Casting values
        #--
        for i in beams_setup.index:
            beam = pd.Series(beams_setup.loc[i])

            # Wave vector direction
            beam["k_dir"] = self.__str_to_list(beam["k_dir"], size = 3)

            # Polarization
            beam["pol_amp"] = self.__str_to_list(beam["pol_amp"], size = 3)

            # Saturation parameter
            if beam["s_0"] != "--": beam["s_0"] = float(beam["s_0"])

            # Waist
            if beam["w"] != "--": beam["w"] = float(beam["w"])

            # Detuning
            if beam["delta"] != "--": beam["delta"] = float(beam["delta"])

            if beam["trans_id"] != "--": beam["trans_id"] = int(beam["trans_id"])

            if beam["sidebands"] != "--": 
                sidebands = self.__str_to_list(beam["sidebands"], size = 2)
                beam["sidebands"] = {"num": int(sidebands[0]), "freq": float(sidebands[1])}

            self._beams[i] = beam
        #--
        #--

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

    # Print all attributes
    def print_all(self):
        print("General Values")
        print("--")
        print(self.general)
        print("--")
        print("Number of beams: %d" % self.num_beams)
        print("--")
        for idx, beam in self.beams.items():
            print("Beam %d" % idx)
            print("--")
            print(beam)
            print("--", end="\n\n")

    # Get Dataframe of beams setup
    def get_setup_dataframe(self):
        indexes = []
        cols = {"k_dir":[], "pol_amp":[], "s_0":[], "w":[], "delta":[], "sidebands":[], "trans_id":[]}

        for idx, beam in self.beams.items():
            indexes.append(idx)
            cols["k_dir"].append(self.__list_to_str(beam["k_dir"]))
            cols["pol_amp"].append(self.__list_to_str(beam["pol_amp"]))
            cols["s_0"].append(beam["s_0"])
            cols["w"].append(beam["w"])
            cols["delta"].append(beam["delta"])
            cols["trans_id"].append(beam["trans_id"])

            if beam["sidebands"] != "--":
                sidebands = "["
                sidebands += str(beam["sidebands"]["num"]) + " "
                sidebands += str(beam["sidebands"]["freq"]) + "]"

            else: sidebands = "--"

            cols["sidebands"].append(sidebands)

        df = pd.DataFrame(cols, index=indexes)
        df.index.name = "idx"
        
        return df[list(cols.keys())]