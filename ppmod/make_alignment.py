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

parser = argparse.ArgumentParser(__doc__,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-j', '--json', help="json file name", default='data.json')
parser.add_argument('-aln', '--alignment', help="alignment file name", default='alignment-file.ali')
parser.add_argument('-p', '--path', help="path to template pdb files", default='../../building_blocks/')
args = parser.parse_args()

def align(p1f, seq, j, jj):
    """Align the template sequence to the target sequence

    Parameters
    ----------
    p1f : str
        Template sequence string
    seq : str
        Target sequence string
    j : int
        How many different template sequences are checked for optimum alignment (different sequences are obtained by sequentially shortening p1f by one AA)  
    jj : int  
        How many different target sequences are checked for optimum alignment (different sequences are obtained by sequentially shortening seq by one AA)

    Returns
    -------
    i : int
        Start of the alignment on template sequence 
    ii : int
        Start of the alignment on target sequence
    """    
    i = 0
    ii = 0
    scoreold = score(p1f, seq)
    for t in range(j):   
       for tt in range(jj): 
           scorenew = score(p1f[t:], seq[tt:])  
           if scorenew > scoreold:
               scoreold = scorenew
               i = t
               ii = tt
    return (i, ii)

def score(p1f, seq):
    """This function scores the quality of the alignment. Higher score means better alignment

    Parameters
    ----------
    p1f : str
        Template sequence string
    seq : str
        Target sequence string
    
    Returns
    -------
    total : int
        Final score of the alignmnet
    """
    total = 0
    for i in range(8):  #check the first ten residues
        if p1f[i] == seq[i]:   #if p1f and seq have a matching residue on position i add a point
            total = total+1  
    return total                          #return final score

def find_pair(pair):
     """This function scores the quality of the alignment. Higher score means better alignment

    Parameters
    ----------
    pair : str
        Name of CC pair
    
    Returns
    -------
    p1 : int
        Id number of pair in 
    p2 : int
    pdbname : str
    """
    for s in d.segments:
        if s['name'] == pair:
            p1 = s[n]['id']-1         #get segment name
            p2 = s[n]['pair_id']-1
            pdbname = s['pdb_template']   #get template pdb file name
            break     
    return p1, p2, pdbname

def selres(i, t):
    """Select i-th residue from topology

    Parameters
    ----------
    i : int
        Residue index
    t : mdtraj topology
        Topology
    
    Returns
    -------
    res : mdtraj residue
        i-th residue
    """
    res = t.select("resid {:d}".format(i))
    return res

def writepdb(i, l, t, p, path):
    """Write the aligned part of the template sequence to pdb

    Parameters
    ----------
    i : int
        Start of the aligning part on the template sequence
    l : int
        Length of the aligned template sequence
    t : mdtraj topology
        topology of the template 
    p : md traj positions
        atom positions in template    
    """
    cl0 = t.chain(0).n_residues  #cl0 length of the first chain of a CC pair
    topsub1 = t.subset(list(range(selres(i, t)[0],selres(i+l, t)[-1])))
    topsub2 = t.subset(list(range(selres(cl0+i, t)[0], selres(cl0+i+l, t)[-1])))   
    topsubjoin = topsub1.join(topsub2)
    coord = np.concatenate((p[0, selres(i, t)[0]:selres(i+l, t)[-1], :]*10, p[0, selres(cl0+i, t)[0]:selres(cl0+i, t)[-1], :]*10), axis=0)
    fpdb = md.formats.PDBTrajectoryFile(path, mode='w')
    fpdb.write(coord, topsubjoin)
    return

d = u.load_json_data(args.json)         #read json file
aln_str = d.entire_sequence

with open(args.alignment, 'w') as f1:
    for pair in d.pairs:  #go through all CC pairs
        p1, p2, pdbname=find_pair(pair[0])
        topology = md.load(os.path.join(args.path, pdbname)).topology   #read topology
        position = md.load(os.path.join(args.path, pdbname)).xyz        #and position from pdb file  
        p1f = topology.to_fasta(0)
        p2f = topology.to_fasta(1)             #convert topology to fasta sequence
        
        #align template structures to target
        #check weather the template sequence is to long and determine alignemnt positions by checking the quality
        #of different alignments
             
        i, ii = align(p1f, d.segments[p1]['sequence'], 6, 6)            
        l = len(min((p1f[i:], d.segments[p1]['sequence'][ii:]), key=len)) 

        #shorten the template sequence if needed and write the topology and the coordinates to a new pdb file
        path = (os.path.join(args.path, pair[0]+'-new.pdb'))        
        writepdb(i, l-1, topology, positions, path)
       
        #write the alignment file taking account previously determined alignemnt position
        count = 0
        f1.write('>P1;{}\n'.format(pdbname.rpartition('.')[0]))    
        f1.write('structureX:{}::A:::::-1.00:-1.00\n'.format(pair[0]+'-new.pdb'))
        
        print('>P1;{}'.format(pdbname.rpartition('.')[0]))    
        print('structureX:{}::A:::::-1.00:-1.00'.format(pair[0]+'-new.pdb'))
        while count < len(aln_str):
            if count == s[p1]['start']-1+ii:
                f1.write(p1f[i:i+l])
                print(p1f[i:i+l])
                count = count+len(p1f[i:i+l])
            elif count == s[p2]['start']-1+ii:
                f1.write(p2f[i:i+l])
                print(p2f[i:i+l])
                count = count+len(p2f[i:i+l])
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
