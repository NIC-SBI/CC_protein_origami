#taken from http://hayd.github.io/2013/dotable-dictionaries/
import os

class Dotable(dict):
    """Generates nested dotable dicts from a json-like object. This makes is possible to write
    a.b[0].c
    """
    __getattr__= dict.__getitem__

    def __init__(self, d):
        self.update(**dict((k, self.parse(v))
                           for k, v in d.items()))

    @classmethod
    def parse(cls, v):
        if isinstance(v, dict):
            return cls(v)
        elif isinstance(v, list):
            return [cls.parse(i) for i in v]
        else:
            return v
        
        
#taken from http://stackoverflow.com/a/13105359/952600
def byteify(input):
    """Takes a jason-like structure and convets unicode str to str"""
    if isinstance(input, dict):
        return {byteify(key):byteify(value) for key,value in input.items()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input        
    
def load_json_data(file_name):
    """Loads data from json and returns a dotable dict with ascii strings"""
    import json
    with open(file_name, 'r') as infile:
        data=json.load(infile)
#    data=byteify(data)
    return Dotable.parse(data)

def pair_ids_from_segments(segments):
    """Return the pair IDs of the segments. The IDs are ZERO based. The matching is based on naming"""
    pair_ids = []
    for n in range(len(segments)):
        for n1 in range(n+1,len(segments)):
            if segments[n]['name']==segments[n1]['pair_name']:
                pair_ids.append((n, n1))
    return pair_ids    
    

import string
import random
import numpy as np
import mdtraj as md
#taken from http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))    


def sequnce_and_knowns(alnfile):
    """Parses the IDs of the known structures and the sequences from the alignment file.
    Returns a tuple:
      (seq, (known1, known2 ...))

    """
    knowns = ()
    with open(alnfile,'r') as f:
        while True:
            line = f.readline()
            if not line: break
            if line[:4] == '>P1;':
                seqname = line[4:].rstrip()
                line = f.readline()            
                if line.split(':')[0] == 'sequence':
                    sequence = seqname
                else:
                     knowns = knowns+(seqname,)
    return (sequence, knowns)                 

def relative_to(file_dir, a_path):
    """Returns the directory of file_dir and appends a_path"""
    #os.path.realpath(
    #os.path.join(os.getcwd(), os.path.dirname(__file__)))
    a_dir = os.path.dirname(file_dir)
    return os.path.join(a_dir,a_path)

def align(template, target, n_template, n_target):
    """Align the template sequence to the target sequence

    Parameters
    ----------
    template : str
        Template sequence string
    target : str
        Target sequence string
    n_template : int
        How many different template sequences are checked for optimum alignment (different sequences are obtained by sequentially shortening p1f by one AA)  
    n_target : int  
        How many different target sequences are checked for optimum alignment (different sequences are obtained by sequentially shortening seq by one AA)

    Returns
    -------
    template_start : int
        Start of the alignment on template sequence 
    target_start : int
        Start of the alignment on target sequence
    """    
    template_start = 0
    target_start = 0
    scoreold = score(template, target)
    for t in range(n_template):   
       for tt in range(n_target): 
           if t == tt: continue
           scorenew = score(template[t:], target[tt:])  
           if scorenew > scoreold:
               scoreold = scorenew               
               template_start = t
               target_start = tt
    return (template_start, target_start)

def score(template, target):
    """This function scores the quality of the alignment. Higher score means better alignment

    Parameters
    ----------
    template : str
        Template sequence string
    target : str
        Target sequence string
    
    Returns
    -------
    total : int
        Final score of the alignmnet
    """
    total = 0
    for i in range(8):  #check the first ten residues
        if template[i] == target[i]:   #if p1f and seq have a matching residue on position i add a point
            total = total+1  
    return total                          #return final score

def find_pair(pair, segments):
    """Find pair's id, pair_id, and its pdb file

    Parameters
    ----------
    pair : str
        Name of CC segments
    
    Returns
    -------
    pair1 : int
        Id number of the CC segment 
    pair2 : int
        Pair_id number of the CC segment
    pdbname : str
        Name of the pdb file where the CC segment is written
    """
    for s in segments:
        if s['name'] == pair:
            pair_1 = s['id']-1         #get segment name
            pair_2 = s['pair_id']-1
            pdbname = s['pdb_template']   #get template pdb file name
            break     
    return pair_1, pair_2, pdbname

def selres(index, topology):
    """Select i-th residue from topology

    Parameters
    ----------
    index : int
        Residue index
    topology : mdtraj topology
        Topology
    
    Returns
    -------
    residue : mdtraj residue
        i-th residue
    """
    residue = topology.select("resid {:d}".format(index))
    return residue

def writepdb(temp_start, temp_len, top, positions, path):
    """Write the aligned part of the template sequence to pdb

    Parameters
    ----------
    temp_start : int
        Start of the aligning part on the template sequence
    templ_len : int
        Length of the aligned template sequence
    top : mdtraj topology
        topology of the template 
    positions : md traj positions
        atom positions in template 
    path : str
        path to folder where the pdb file will be saved   
    """
    chainlength = top.chain(0).n_residues  #length of the first chain of a CC pair
    topsubset1 = top.subset(list(range(selres(temp_start, top)[0],selres(temp_start + temp_len, top)[-1]))) #topology of atoms in the first aligned chain of the CC segment
    topsubset2 = top.subset(list(range(selres(chainlength + temp_start, top)[0], selres(chainlength +temp_start + temp_len, top)[-1]))) #topology of atoms in the first aligned chain of the CC segment
    topsubjoin = topsubset1.join(topsubset2) # join topology subsets

    coordinate1 = positions[0, selres(temp_start, top)[0]:selres(temp_start + temp_len, top)[-1], :]*10 #position of atoms in the first subset
    coordinate2 = positions[0, selres(chainlength + temp_start, top)[0]:selres(chainlength + temp_start + temp_len, top)[-1], :]*10 #position of atoms in the second subset
    coordinate = np.concatenate((coordinate1, coordinate2), axis=0) #join position subsets

    fpdb = md.formats.PDBTrajectoryFile(path, mode='w') 
    fpdb.write(coordinate, topsubjoin) #write to pdb
    return
