#! /usr/bin/env python
"""
Creates an alpha-helix given an input sequence 
"""
import chimera
import Midas
#command line arguments handling 
import argparse 

parser = argparse.ArgumentParser(description=__doc__, #uses the scripts docstring
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-s','--seq', 
          help='Aminoacid sequence of the helix to build.', 
          type =str)

parser.add_argument('-o','--out-file', 
          help='Output file name.', 
          default='helix.pdb', type=str)

a = parser.parse_args()

#trim quotes and spaces
a.seq=a.seq.replace('"','')
a.seq=a.seq.replace("'",'')
a.seq=a.seq.replace(" ",'')


from BuildStructure import placePeptide

#rotlib Richardson.mode is faster
placePeptide(a.seq, [(-57, -47)] * len(a.seq), rotlib="Richardson.mode")
#get the model object (it's the last object created)

mol = chimera.openModels.list()[-1]
Midas.write(mol, None, a.out_file)

