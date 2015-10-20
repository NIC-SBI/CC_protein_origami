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
import os
import subprocess

parser = argparse.ArgumentParser(description=__doc__, #uses the scripts docstring
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a','--seq', 
          help='Aminoacid sequence of the helix to build.', 
          type =str)

parser.add_argument('-o','--out-file', 
          help='Output file name.', 
          default='helix.pdb', type=str)

parser.add_argument('-d','--debug', 
          help='Debug mode (shows chimera debug output).', 
          default=False, action='store_true')

parser.add_argument('-f','--fake', 
          help='Only print command but do not execute it.', 
          default=False, action='store_true')

a = parser.parse_args()

chimera_path = 'chimera'
script_dir = os.path.dirname(os.path.realpath(__file__))
chimera_flags = "--nogui"
if a.debug:
  chimera_flags += " --debug"

cmd = "bash -c \"{path} {flags} --script '{script_dir}/make_helix_chimera.py --out-file={out_file} --seq={seq}'\"".format(
       path=chimera_path, flags=chimera_flags, script_dir=script_dir, out_file=a.out_file, seq=a.seq)
if a.fake:
    print cmd
else:
    subprocess.call(cmd, shell=True)