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
        for pair in d.pairs:  #go through all CC pairs

            pair_1_id, pair_2_id, pdbname = u.find_pair(pair['pair'][0], d.segments)
            pair_name = d.segments[pair_1_id].name + "_" + d.segments[pair_2_id].name 

            md_obj = md.load(os.path.join(args.path, pdbname)) 
            topology = md_obj.topology   #read topology
            position = md_obj.xyz        #and position from pdb file  

            pair_1_seq = u.mdtraj_to_fasta(topology,0)
            pair_2_seq = u.mdtraj_to_fasta(topology,1)           #convert topology to fasta sequence
            
            #align template structures to target check weather the template sequence is to long and determine alignemnt positions by checking the quality of different alignments
            template_start, target_start = u.align(pair_1_seq, d.segments[pair_1_id]['sequence'], 6, 6)            
            #compare the length of aligned template and target sequence and get the length of teh shorter one
            min_length = len(min((pair_1_seq[template_start:], d.segments[pair_1_id]['sequence'][target_start:]), key=len)) 
            if len(pair_1_seq) > min_length:
                #shorten the template sequence if needed and write the topology and the coordinates to a new pdb file
                pdbname = pdbname.replace('.pdb','') + '_' + str(template_start) + '_' + str(template_start + min_length) + '.pdb'
                path = (os.path.join(args.path, pdbname)) #path to new pdb files      
                u.writepdb(template_start, min_length-1, topology, position, path)
           
            #write the alignment file taking into account previously determined alignment position
            count = 0
            f1.write('>P1;{}\n'.format(pair_name))    
            f1.write('structureX:{}::A::B:::-1.00:-1.00\n'.format(pdbname))
            
            print('>P1;{}'.format(pair_name))    
            print('structureX:{}::A::B:::-1.00:-1.00'.format(pdbname))
            while count < len(aln_str):
                if count == d.segments[pair_1_id]['start']-1 + target_start:
                    f1.write(pair_1_seq[template_start : template_start + min_length] + "/")
                    print(pair_1_seq[template_start : template_start + min_length] + "/")
                    count = count+len(pair_1_seq[template_start : template_start + min_length]) + 1
                elif count == d.segments[pair_2_id]['start']-1 + target_start:
                    f1.write(pair_2_seq[template_start : template_start + min_length])
                    print(pair_2_seq[template_start : template_start + min_length])
                    count = count+len(pair_2_seq[template_start : template_start + min_length])
                else:
                    f1.write("-")
                    print("-", end="")
                    count = count +1
            print('*')
            print("")
            f1.write('*\n')
            f1.write('\n')
        f1.write('>P1;{}\n'.format(d.name))
        f1.write('sequence:{}:1:A::B::: 0.00:0.00\n'.format(d.name))
        f1.write(aln_str)
        f1.write('*\n')
        print('>P1;{}\n'.format(d.name), end="")
        print('sequence:{}:1:A::B::: 0.00:0.00\n'.format(d.name), end="")
        print(aln_str, end="")
        print('*\n', end="")
