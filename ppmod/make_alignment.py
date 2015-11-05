#!/usr/bin/env python
"""
Generates an alignment file to be used with Modeller

The target sequence is read from .json file. Json file must also contain the information about pdf files with template
structure. Program aligns dtemplate sequences to target sequence. Any residual residues from template sequence are removed
and the new structure is written to a pdb file which is used in homology modelling. Finally the alignment file is generated.
It is advised to check the file before use in subsequent modelling.
"""

from __future__ import print_function, absolute_import
import mdtraj as md
import argparse
import re
import utils as u
import os
import numpy as np
from modeller import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-j', '--json', help="json file name", default='data.json')
    parser.add_argument('-aln', '--alignment', help="alignment file name", default='alignment-file.ali')
    parser.add_argument('-p', '--path', help="path to template pdb files", default='../../building_blocks/')
    args = parser.parse_args()

    d = u.load_json_data(args.json)         #read json file
    aln_str = d.entire_sequence

    with open(args.alignment, 'w') as f1:       
        for pair, pair2 in d.pairs:  #go through all CC pairs
            p1, p2, pdbname = u.find_pair(pair, d.segments)
            topology = md.load(os.path.join(args.path, pdbname)).topology   #read topology
            position = md.load(os.path.join(args.path, pdbname)).xyz        #and position from pdb file  
            p1f = topology.to_fasta(0)
            p2f = topology.to_fasta(1)             #convert topology to fasta sequence
            
            #align template structures to target check weather the template sequence is to long and determine alignemnt positions by checking the quality of different alignments
            i, ii = u.align(p1f, d.segments[p1]['sequence'], 6, 6)            
            
 

            #shorten the template sequence if needed and write the topology and the coordinates to a new pdb file
            l = len(min((p1f[i:], d.segments[p1]['sequence'][ii:]), key=len)) #compare the length of aligned template and target sequence and get the length of teh shorter one
            path = (os.path.join(args.path, pair+'-new.pdb')) #path to new pdb files      
            u.writepdb(i, l-1, topology, position, path)
           
            #write the alignment file taking into account previously determined alignment position
            count = 0
            f1.write('>P1;{}\n'.format(pair))    
            f1.write('structureX:{}::A:::::-1.00:-1.00\n'.format(pair + '-new.pdb'))
            
            print('>P1;{}'.format(pair))    
            print('structureX:{}::A:::::-1.00:-1.00'.format(pair + '-new.pdb'))
            while count < len(aln_str):
                if count == d.segments[p1]['start']-1 + ii:
                    f1.write(p1f[i : i + l])
                    print(p1f[i : i + l])
                    count = count+len(p1f[i : i + l])
                elif count == d.segments[p2]['start']-1 + ii:
                    f1.write(p2f[i : i + l])
                    print(p2f[i : i + l])
                    count = count+len(p2f[i : i + l])
                else:
                    f1.write("-")
                    print("-", end="")
                    count = count +1
            print('*')
            print("")
            f1.write('*\n')
            f1.write('\n')
        f1.write('>P1;{}\n'.format(d.name))
        f1.write('sequence:{}:1:A::::: 0.00:0.00\n'.format(d.name))
        f1.write(aln_str)
        f1.write('*\n')
        print('>P1;{}\n'.format(d.name), end="")
        print('sequence:{}:1:A::::: 0.00:0.00\n'.format(d.name), end="")
        print(aln_str, end="")
        print('*\n', end="")
