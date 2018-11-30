"""
Script for writting transfer models (.inp) for the transfer program in fortran ./recub
"""

#==============================================================================
#  Python script to run radiative transfer models
#==============================================================================

import subprocess
import os
import numpy as np
import string
import glob
from copy import deepcopy
from ferpy import radtransf_modelpars_v6 as radtransf_modelpars
import re
import gc

# =============================================================================
# Model Runner
# =============================================================================
def model_runner(path, model_name):
    """
    To run .inp models
    Path ending with /
    model_name without .inp
    """
    # Starting workdir
    origdir =  os.getcwd()

    # Radiative transfer code directory
    inpdir = '/Users/frico/Documents/Ed/transf/'
    # Program directory
    progdir = inpdir+'program/1drec/'
    # Status file
    statusdir = progdir+'STATUS.REC'

    # Model input directory
    modeldir = path
    if not os.path.exists(modeldir):
        os.makedirs(modeldir)
    # Input model
    modell = model_name+'.inp'
    if not os.path.exists(modeldir+modell+'_.res'):
        # Running model
        os.chdir(progdir)
        ### Modify status file
        statusfile = open(statusdir, 'w')
        statusfile.write("'"+modeldir+modell+"'") # modeldir/modell.inp
        statusfile.close()
        ### To avoid problem with dylib and macOS
        subprocess.call('install_name_tool -change @rpath/libiomp5.dylib /opt/intel/compilers_and_libraries_2017.2.163/mac/compiler/lib/libiomp5.dylib ./recub', shell='True')
        ### Running the program
        subprocess.call('./recub', shell='True')
        os.chdir(origdir)
    else:
        print('\nModel already ran: %s (change model name or delete out files)' % modell)
        

    