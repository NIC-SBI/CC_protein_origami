#! /usr/bin/env python
"""
Creates an alpha-helix given an input sequence (or a json config file). 

Chimera needs to be installed and available on the path
"""

#names of columns, that have the fraction of native contacts for each cc segment

from __future__ import print_function, division, absolute_import, unicode_literals
import ppmod.utils as u

#command line arguments handling 
import argparse 
import os
import subprocess

cc_names = ['cc_0','cc_1','cc_2','cc_3','cc_4','cc_5']

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

parser.add_argument('-j','--json', 
          help='Json file from which to load aminoacid sequnce.', 
          type =str)

parser.add_argument('-e','--auto-exit', 
          help='Exit at the end of script', 
          type =bool, default=True)

a = parser.parse_args()


if not a.json is None:
    d = u.load_json_data(a.json)
    a.seq = d.entire_sequence

chimera_path = 'chimera'
script_dir = os.path.dirname(os.path.realpath(__file__))
chimera_flags = ""#"--nogui"




if a.debug:
  chimera_flags += " --debug"

cmd = "bash -c \"{path} {flags} --script '{script_dir}/make_helix_chimera.py --out-file={out_file} --seq={seq} --auto-exit={ae} '\"".format(
       path=chimera_path, flags=chimera_flags, script_dir=script_dir, out_file=a.out_file, seq=a.seq, ae=a.auto_exit)

if a.fake:
    print(cmd)
else:
    subprocess.call(cmd, shell=True)
