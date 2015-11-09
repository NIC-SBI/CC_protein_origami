
# -*- coding: utf-8 -*-
"""
Tools to generate a json file with infomration about segments.
"""
from __future__ import print_function
import utils as u
import yaml
import six

def load_pairs(yaml_str):
    """Loads a yaml string and returns a list of pairs"""
    pairs = yaml.load(yaml_str)
    
    #split the pairs and chains if they are seperated by a colon    
    for n in range(len(pairs)):
        #if it's a calumn split it        
        if isinstance(pairs[n]['pair'], six.string_types):    
            pairs[n]['pair']= pairs[n]['pair'].split(':')
        #strip whitespace    
        pairs[n]['pair'] = [p.strip() for p in pairs[n]['pair']]
        
        #do the same for chains
        if isinstance(pairs[n]['chains'], six.string_types):                            
            pairs[n]['chains']= pairs[n]['chains'].split(':')
            
        #strip whitespace    
        pairs[n]['chains'] = [c.strip() for c in pairs[n]['chains']]
    return pairs
    
    
def find_pair_by_segment_name(pairs, segment_name):
    """Given a segment name (for example p3) return the matching pair structue. """
            #find the pair name
    for n_pair, pair in enumerate(pairs):
        if segment_name in pair['pair']:
            return pair 

def get_other_segment_name(pair, segment_name):
    """returns the name of the other segment"""
    pair_name = pair['pair'][:]
    pair_name.remove(segment_name)                        
    pair_name=pair_name[0]
    return pair_name
        
def generate_json(name, entire_sequence, segments_str, pairs, out_name=None):
    """generates json file with information about the polypeptide polyhedra.
    Segment parings are calculated, etc...

    Parameters
    ----------
    name : str
        The name of the polyhedral design.
    entire_sequence : str
        The whole sequnce, including his tags etc..
    segments_str: str
        A string of segment sequnces and segment names. The names must be tab separated.
    pairs : dict
        

    Returns
    -------
      Does not return anything. By default generates an file with name {name}.json
    """  
    entire_sequence=entire_sequence.replace(' ','')
    segments = segments_str.split('\n')    
    #ignore the linkers
    segments = list(filter(lambda s: len(s.strip())>10, segments)) 
    
    
    # split into seq name pairs
    for n, seg in enumerate(segments):
        spl = seg.split('\t')
        seq = spl[0].replace(" ","")
        name = spl[-1].strip()
    
        seq = seq.replace(" ","")
        
        #has a segment allready been found?
        pos_in_seq=0
        for n1 in range(n):
            print(n1, n-1,name)
            print(name)
            print(segments[n1]['start'],segments[n1]['len'])
            if segments[n1]['name']==name:
                pos_in_seq=segments[n1]['start']+segments[n1]['len']
                

        start = entire_sequence.find(seq,pos_in_seq)
        
        pair = find_pair_by_segment_name(pairs, name)
        pair_name = get_other_segment_name(pair, name)
        
        
        segments[n]={"id":n+1, "name": name, 
                     "sequence" :seq, 
                     "start" : start+1,
                     "len" : len(seq),
                     "end" : start+len(seq),
                     "pair_id" : -1,
                     "pair_name": pair_name,
                     "pair_type": pair['type'],
                     "pdb_template" : pair['template'],
                     "pdb_chain" : ""
                    }
       
    #find tha pair IDs
    for n in range(len(segments)):
        for n1 in range(n+1,len(segments)):
            if segments[n]['name']==segments[n1]['pair_name']:
                segments[n]['pair_id']=segments[n1]['id']      
                segments[n1]['pair_id']=segments[n]['id']


if __name__ == "main":
    import argparse    
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument("-j","--json",default='out/data.json',help='specify path to json file')
    parser.add_argument('-i','--initial-model', help='Initial model in pdb format')
    parser.add_argument('-a','--alnfile', help='secify path to alignemnt file')
