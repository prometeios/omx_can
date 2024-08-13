#!/usr/bin/env python3

#### >>>> Description of Code <<<< ####
__description__= """
@brief: Canonical algorithm for get magnetic local minima
@author: Cha
@email: caca518@kaist.ac.kr
@date: 2024.08.03
@version: 0.0
"""
#### >>>> Description of Code End <<<< ####


#### >>>> Modules and constant <<<< ####
import os

from . import omx_io

THRESHOLD_MAGMOM = 0.05
#### >>>> Modules and constant End <<<< ####


#### >>>> Class Definition <<<< ####

#### >>>> Class Definition End <<<< ####


#### >>>> Function Definition <<<< ####

def is_mag(outfile):
    # Check if result has magnetic moments
    output = omx_io.OmxOutput(outfile)
    is_mag = output.magmom > THRESHOLD_MAGMOM
    return is_mag
    

def alter_init_magmom(inputfile, numcores):
    # Change the initial magnetic moment until converge to magnetic local minima
    
    # Check initial guess of magnetic moment
    input = omx_io.OmxInput(inputfile)
    INIT_MAG = input.structure.spin[0][0]
    sysname = input.options.sysname
    
    # Set workdir to the input file directory
    workpath = os.path.dirname(inputfile)
    input.options.workdir = workpath

    
    # TODO: it assumes magnetic materials are placed first. Change it to be more general
    
    # TODO: Guess the initial magnetic moment based on stoichiometry
    
    ## Generate an list of initial magnetic moments to check in likely order
    MAX_SEARCH = 10 # Maximum number of search
    MIN_STEP = 0.2 # Minimum steps to change magnetic moment
    MAX_MAG = 5.0 # Maximum magnetic moment
    
    max_mag = 1.5 * INIT_MAG
    if max_mag > MAX_MAG:
        max_mag = MAX_MAG
    min_mag = 0.5 * INIT_MAG
    nsearch = MAX_SEARCH
    step = (max_mag - min_mag) / (MAX_SEARCH-1)
    if step < MIN_STEP:
        step = MIN_STEP
        nsearch = int((max_mag - min_mag) / MIN_STEP) + 1
    
    magmom_list = [min_mag + i * step for i in range(nsearch)]
    
    # Sort the list in the order of difference from INIT_MAG
    magmom_list = sorted(magmom_list, key=lambda x: abs(x - INIT_MAG))
    
    
    
    ## Loop over the list of magnetic moments
    for magmom in magmom_list:
        input.structure.spin[0][0] = magmom
        input.write_datfile(inputfile) # This will write the dat file in users working directory. Not exactly overwriting the input file
        # Run the calculation
        os.system(f"mpirun -np {numcores} openmx {inputfile}")
        if is_mag(f"{sysname}.out"):
            print(f"Initial magnetic moment {magmom} is magnetic")
            break
    print("Altering initial magnetic moment is done. No magnetic moment is found.")
    

def alter_mixing():
    # Change the mixing parameters and algorithms
    pass

def alter_U(inputfile, numcores):
    # Start from the DFT+U calculation
    input = omx_io.OmxInput(inputfile)
 
    # Set workdir to the input file directory
    workpath = os.path.dirname(inputfile)
    input.options.workdir = workpath
    transition_metal = input.structure.positions[0][0]
    # TODO: Currently it assumes the first atom is the transition metal. Change it to be more general

    # Check the U value that makes magnetic local minima
    MAX_U = 5
    for U in range(MAX_U):
        # Set the U value
        input.options.is_restart = False
        input.options.is_dft_u_correction = True
        input.options.hubbard_u_value[transition_metal]['1d'] = U
        input.write_datfile(inputfile)

        # Run the calculation
        os.system(f"mpirun -np {numcores} openmx {inputfile}")
        # Check if magnetic local minima is found
        if is_mag(f"{input.options.sysname}.out"):
            print(f"System is magnetic at hubbard U value {U}")
            break
        # If found, break the loop
    # If not found, finish the function
    if U == MAX_U:
        if not is_mag(f"{input.options.sysname}.out"):
            print("No magnetic local minima is found.")
            return
    

    # Diminish the U value until 0, escape when magnetic local miniam vanished
    USTEP = 0.5
    for U in range(U, -USTEP, -USTEP):
        # Set the U value
        input.options.is_restart = True
        input.options.hubbard_u_value[transition_metal]['1d'] = U
        input.write_datfile(inputfile)

        # Run the calculation
        os.system(f"mpirun -np {numcores} openmx {inputfile}")
        # Check if magnetic local minima is found
        if not is_mag(f"{input.options.sysname}.out"):
            print(f"System is non-magnetic at hubbard U value {U}")
            break
        if U == 0:
            print("magnetic local minima is found.")
    return

def alter_cons():
    # Start from the constrained calculation
    pass

#### >>>> Function Definition End <<<< ####


#### >>>> Main Function <<<< ####
def main():
    pass
#### >>>> Main Function End <<<< ####

if __name__ == '__main__':
    main()
