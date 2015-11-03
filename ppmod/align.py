#!/usr/bin/env python
"""
Generates an alignment file to be used with Modeller

The target sequence is read from .json file. Json file must also contain the information about pdf files with template
structure. Program aligns dtemplate sequences to target sequence. Any residual residues from template sequence are removed
and the new structure is written to a pdb file which is used in homology modelling. Finally the alignment file is generated.
It is advised to check the file before use in subsequent modelling.
"""

from __future__ import print_function
import mdtraj as md
import argparse
import re
import utils as u
import os
import numpy as np
from modeller import *

parser = argparse.ArgumentParser('Generate alignment file from information in JSON file',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-j', '--json', help="json file name",default='data.json')
parser.add_argument('-aln', '--alignment', help="alignment file name",default='alignment-file.ali')
parser.add_argument('-p', '--path', help="path to template pdb files",default='building_blocks/')
args = parser.parse_args()

def score(p1f,seq,flag):
    "This function scores the quality of the alignment. Higher score means better alignment"
    total = 0
    if flag == 0:
        for i in list(range(8)):  #check the first ten residues
            if p1f[i] == seq[i]:   #if p1f and seq have a matching residue on position i add a point
                total = total+1  
    else:
        for i in list(range(-1,-9,-1)):  #check the last ten residues
            if p1f[i] == seq[i]:          #if p1f and seq have a matching residue on position i add a point
                total = total+1 
    return total                          #return final score



d = u.load_json_data(args.json)         #read json file
aln_str = d.entire_sequence

with open(args.alignment,'w') as f1:
    for i in range(len(d.pairs)):  #go through all CC pairs
        for n in range(len(d.segments)):
            if d.segments[n]['name'] == d.pairs[i][0]:
                p1 = d.segments[n]['id']-1         #get segment name
                p2 = d.segments[n]['pair_id']-1
                pdbname=d.segments[n]['pdb_template']   #get template pdb file name
                break
        topology = md.load(os.path.join(args.path, pdbname)).topology   #read topology
        position = md.load(os.path.join(args.path, pdbname)).xyz        #and position from pdb file  
        p1f = topology.to_fasta(0)
        p2f = topology.to_fasta(1)             #convert topology to fasta sequence
        
        #align template structures to target
        i = 0            
        j = -1
        k = 0
        l = -1
        ii = 0
        jj = -1
        kk = 0
        ll = -1
        #check weather the template sequence is to long and determine alignemnt positions by checking the quality
        #of different alignments
        scoreold = score(p1f[i:],d.segments[p1]['sequence'][ii:],0)
        for t in list(range(6)):   
            for tt in list(range(6)): 
                scorenew = score(p1f[t:],d.segments[p1]['sequence'][tt:],0)  
                if scorenew > scoreold:
                    scoreold = scorenew
                    i = t
                    ii = tt
                    
        scoreold = score(p1f[:],d.segments[p1]['sequence'][:],1)
        for t in list(range(-1,-7,-1)):
            for tt in list(range(-1,-7,-1)):
                scorenew = score(p1f[:t],d.segments[p1]['sequence'][:tt],1)
                if scorenew > scoreold:
                    scoreold = scorenew
                    j = t
                    jj = tt
       
        #shorten the template sequence if needed and write the topology and the coordinates to a new pdb file        
        topsub1 = topology.subset(list(range(topology.select("resid {:d}".format(i))[0],topology.select("resid {:d}".format(topology.chain(0).n_residues+j))[-1])))
        topsub2 = topology.subset(list(range(topology.select("resid {:d}".format(topology.chain(0).n_residues+i))[0],topology.select("resid {:d}".format(topology.n_residues+j))[-1])))   
        topsubj = topsub1.join(topsub2)
        coord = np.concatenate((position[0,topology.select("resid {:d}".format(i))[0]:topology.select("resid {:d}".format(topology.chain(0).n_residues+j))[-1],:]*10,position[0,topology.select("resid {:d}".format(topology.chain(0).n_residues+i))[0]:topology.select("resid {:d}".format(topology.n_residues+j))[-1],:]*10),axis=0)
        fpdb = md.formats.PDBTrajectoryFile(os.path.join(args.path, pdbname.rpartition('.')[0]+'-new.pdb'), mode='w', force_overwrite=True)
        fpdb.write(coord,topsubj)
        
        #write the alignment file taking account previously determined alignemnt position
        count = 0
        f1.write('>P1;{}\n'.format(pdbname.rpartition('.')[0]))    
        f1.write('structureX:{}::A:::::-1.00:-1.00\n'.format(pdbname.rpartition('.')[0]+'-new.pdb'))
        
        print('>P1;{}'.format(pdbname.rpartition('.')[0]))    
        print('structureX:{}::A:::::-1.00:-1.00'.format(pdbname.rpartition('.')[0]+'-new.pdb'))
        while count < len(aln_str):
            if count == d.segments[p1]['start']-1+ii:
                f1.write(p1f[i:len(p1f)+1+j])
                print(p1f[i:len(p1f)+1+j])
                count=count+len(p1f[i:len(p1f)+1+j])
            elif count == d.segments[p2]['start']-1+ii:
                f1.write(p2f[i:len(p2f)+1+j])
                print(p2f[i:len(p2f)+1+j])
                count = count+len(p2f[i:len(p2f)+1+j])
            else:
                f1.write("-")
                print("-",end="")
                count = count +1
        print('*')
        print("")
        f1.write('*\n')
        f1.write('\n')
    f1.write('>P1;{}\n'.format(d.name))
    f1.write('sequence:{}:1:A::::: 0.00:0.00\n'.format(d.name))
    f1.write(aln_str)
    f1.write('*\n')
    print('>P1;{}\n'.format(d.name),end="")
    print('sequence:{}:1:A::::: 0.00:0.00\n'.format(d.name),end="")
    print(aln_str,end="")
    print('*\n',end="")
