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
           if t == tt: continue
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

def find_pair(pair, segments):
    """Find pair's id, pair_id, and its pdb file

    Parameters
    ----------
    pair : str
        Name of CC segments
    
    Returns
    -------
    p1 : int
        Id number of the CC segment 
    p2 : int
        Pair_id number of the CC segment
    pdbname : str
        Name of the pdb file where the CC segment is written
    """
    for s in segments:
        if s['name'] == pair:
            p1 = s['id']-1         #get segment name
            p2 = s['pair_id']-1
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
    topsub1 = t.subset(list(range(selres(i, t)[0],selres(i + l, t)[-1]))) #topology of atoms in the first aligned chain of the CC segment
    topsub2 = t.subset(list(range(selres(cl0 + i, t)[0], selres(cl0 +i + l, t)[-1]))) #topology of atoms in the first aligned chain of the CC segment
    topsubjoin = topsub1.join(topsub2) # join topology subsets
    coord1 = p[0, selres(i, t)[0]:selres(i + l, t)[-1], :]*10 #position of atoms in the first subset
    coord2 = p[0, selres(cl0 + i, t)[0]:selres(cl0 + i + l, t)[-1], :]*10 #position of atoms in the second subset
    coord = np.concatenate((coord1, coord2), axis=0) #join position subsets
    fpdb = md.formats.PDBTrajectoryFile(path, mode='w') 
    fpdb.write(coord, topsubjoin) #write to pdb
    return
