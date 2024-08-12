#!/usr/bin/env python3

#### >>>> Description of Code <<<< ####
__description__= """
  This code recalculate openMX first principle calculation until Magnetic moment is observed.
  Firstly, collinear case with strategy: change initial magnetic moment, change kgrid, change mixing is considered
"""
__author__ = "Cha"
__email__ = "caca518@kaist.ac.kr"
__version__ = 0.0
#### >>>> Description of Code End <<<< ####

#### >>>> Modules and constant <<<< ####
import os, sys, subprocess
from itertools import combinations
import numpy as np
import glob

from janus.dft import openmx as omx

NCORE = 2
COMMAND = f"mpirun -np {NCORE} openmx *.dat"
THRESHOLD_MAGMOM = 0.05
#### >>>> Modules and constant End <<<< ####



#### >>>> Function Definition <<<< ####
def is_magmom_observed():
    files = glob.glob("*.out")
    if not files:
        raise FileNotFoundError("No output file is found")
    if len(files) > 1:
        raise ValueError("Multiple output files are found")
    if omx.OmxOutput(files[0]).magmom > THRESHOLD_MAGMOM:
        return True
    return False

def strategy_initial_magmom():
    # define searching range usign initial magnetic moments
    searching_range: list = get_magmom_searching_range(omx.OmxInput("*.dat").structure.spin[0])
    # change initial magnetic moments using according to searching range
    for magmom in searching_range:
        change_initial_magmom(magmom)
        # run calculation
        subprocess.run(COMMAND, shell=True, check=True)
        # check magnetic moment
        if is_magmom_observed():
            return True
    # return True if magnetic moment is observed
    return False


def change_initial_magmom(magmom):
    input = omx.OmxInput("*.dat")
    input.structure.spin[0] = magmom
    input.write_datfile()
    return

def get_magmom_searching_range(magmom):
    # below 5 over 0.5. No smaller steps than 0.1
    # start from vicinity of original magmom

    max_magmom = magmom*1.5 if magmom*1.5 < 5 else 5
    min_magmom = magmom*0.5 if magmom*0.5 > 0.5 else 0.5

    # Generate a list of 5 numbers equally distributed from min_magmom to max_magmom
    searching_range = np.linspace(min_magmom, max_magmom, 5)

    # Sort the list by the absolute difference from magmom
    searching_range = sorted(searching_range, key=lambda x: abs(x - magmom))
    return searching_range


def strategy_kgrid():
    # define searching range using kgrid
    searching_range: list = get_kgrid_searching_range(omx.OmxInput("*.dat").options.kpoints[0])
    # change kgrid using according to searching range
    for kgrid in searching_range:
        change_kgrid(kgrid)
        # run calculation
        subprocess.run(COMMAND, shell=True, check=True)
        # check magnetic moment
        if is_magmom_observed():
            return True
        # return True if magnetic moment is observed
    return False

def get_kgrid_searching_range(kgrid):
    # just do the few steps more of it.
    searching_range = [kgrid, kgrid + 1, kgrid + 2]
    return searching_range

def change_kgrid(kgrid):
    input = omx.OmxInput("*.dat")
    input.options.kpoints = [kgrid, kgrid, 1]
    return

def strategy_mixing():
    # mixing type
    # RMM-DIISK and Kerker (good for metals)
    #
    # Generally, larger mixinghistory(30-50) and everypulay 1, larger mixing_start_pulay(10 -30), larger electronic temperature (1000)
    # Also, Kerker factor supress the charge sloshing and slowing convergence. 
    # But note that we are not dealing with the charge sloshing problem.
    #
    changes = [(change_mixing_type, "RMM-DIISK"), (change_mixing_history, 30), (change_mixing_start_pulay, 10), (change_every_pulay, 1), (change_electronic_temperature, 1000), (change_kerker_factor, 0.5)]
    iterate_changes(changes)


    pass

def iterate_changes(functions):
    # Takes a list of tuples (functions, arguments) and iterate all possible combinations of the functions
    # Change options separately
    for func in functions:
        func[0](*func[1:])
        subprocess.run(COMMAND, shell=True, check=True)

    # Change options simultaneously
    for i in range(2, len(functions) + 1):
        for combo in combinations(functions, i):
            for func in combo:
                func[0](*func[1:])
                subprocess.run(COMMAND, shell=True, check=True)

def change_mixing_type():
    pass
def change_mixing_history():
    pass
def change_mixing_start_pulay():
    pass
def change_every_pulay():
    pass
def change_electronic_temperature():
    pass
def change_kerker_factor():
    pass
    

def change_mixing():
    pass


#### >>>> Function Definition End <<<< ####

#### >>>> Main Function <<<< ####
def main():
    # run calculation for first time
    subprocess.run(COMMAND, shell=True, check=True)
    if is_magmom_observed():
        print("Magnetic moment is observed")
        return
    if strategy_initial_magmom():
        print("Magnetic moment is observed by changing initial magnetic moment")
        return
    if strategy_kgrid():
        print("Magnetic moment is observed by changing kgrid")
        return
    #if strategy_mixing():
        #print("Magnetic moment is observed by changing mixing")
        #return
    print("Magnetic moment is not observed")
    
#### >>>> Main Function End <<<< ####

if __name__ == '__main__':
    main()