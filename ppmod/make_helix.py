#! /usr/bin/env python
"""
Creates an alpha-helix given an input sequence (or a json config file). 

Chimera needs to be installed and available on the path
"""

#names of columns, that have the fraction of native contacts for each cc segment
cc_names = ['cc_0','cc_1','cc_2','cc_3','cc_4','cc_5']

import pandas as pd
import numpy as np
#command line arguments handling 
import argparse 

parser = argparse.ArgumentParser(description=__doc__, #uses the scripts docstring
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a','--sequence', 
          help='Aminoacid sequence of the helix to build.', 
          type =str)

parser.add_argument('-o','--out-file', 
          help='Output file name.', 
          default='helix.pdb', type=str)

a = parser.parse_args()
print a