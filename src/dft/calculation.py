#!/usr/bin/env python3

#### >>>> Description of Code <<<< ####
__description__= """
  Methods to excute OpenMX calculation. Contain automatic recalculation with modified option
"""
__author__ = "Cha"
__email__ = "caca518@kaist.ac.kr"
__version__ = 0.0
#### >>>> Description of Code End <<<< ####

#### >>>> Modules and constant <<<< ####
import os

from janus.dft import openmx as omx


#### >>>> Modules and constant End <<<< ####

#### >>>> Methods <<<< ####
def run_openmx(datfile, core_number):
    """ Run OpenMX calculation"""
    # mpirun with core number
    os.system("mpirun -np %d openmx %s" % (core_number, datfile))

def modify_init_magmom(datfile, magmom):
    """ Modify initial magmom"""
    # Read datfile
    omx_input = omx.OmxInput(datfile)
    # Modify magmom
    omx_input.options.init_magmom = magmom
    # Write datfile
    omx_input.write_datfile()
    
def check_

#### >>>> Class <<<< ####