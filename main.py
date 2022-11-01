# Libraries and Modules
from controller import Controller
from datetime import datetime as dt
from model import Results
from model import Simulation
from view import View
import numpy as np
import sys

# Controller Object
ctrler = Controller()
  
# Main
if __name__ == '__main__':  
    # Initial menu
    ctrler.main_menu()

