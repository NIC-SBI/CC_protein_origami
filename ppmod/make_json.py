
# -*- coding: utf-8 -*-
"""
Tools to generate a json file with infomration about segments.
"""
from __future__ import print_function
import utils as u
import yaml
import six
import re

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


def get_segment_chain(pair, segment_name):
    """Returns the chain belonging to the pair. 
       If homodimers, returns the first chain, but this must be fixed later."""
    ind = pair['pair'].index(segment_name)    
    return pair['chains'][ind]

        
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
    segments = list(filter(lambda s: len(s.strip())>15, segments)) 
    pos_in_seq=0
    # split into seq name pairs
    for n, seg in enumerate(segments):

        #split on ' ', ',' ';','|'    
        spl = re.split('\s*[ |;,]\s*',seg)        
        seq = spl[0].replace(" ","")
        seg_name = spl[-1].strip()
        seq = seq.replace(" ","")        

        start = entire_sequence.find(seq, pos_in_seq)
        pos_in_seq = start + len(seq)

        pair = find_pair_by_segment_name(pairs, seg_name)
        pair_name = get_other_segment_name(pair, seg_name)
        chain = get_segment_chain(pair, seg_name)
        
        segments[n]={"id":n+1, "name": seg_name, 
                     "sequence" :seq, 
                     "start" : start+1,
                     "len" : len(seq),
                     "end" : start+len(seq),
                     "pair_id" : -1,
                     "pair_name": pair_name,
                     "pair_type": pair['type'],
                     "pdb_template" : pair['template'],
                     "pdb_chain" : chain
                    }
       
    #find tha pair IDs
    for n1 in range(len(segments)):
        for n2 in range(n1+1,len(segments)):
            if segments[n1]['name'] == segments[n2]['pair_name']:
                segments[n1]['pair_id'] = segments[n2]['id']      
                segments[n2]['pair_id'] = segments[n1]['id']
                
                #correct chains for homodimers
                if segments[n1]['name'] == segments[n2]['name']:
                    #TODO: This is a hack! Now homodimers only support chains A and B                    
                    segments[n1]['pdb_chain'] = "A"
                    segments[n2]['pdb_chain'] = "B"
                
    import json            
    if out_name is None:    
        out_name = name+".json"
    

    #don't remember if pair_types are used anywhere
    pair_types = [p['type'] for p in pairs]    
    
    with open(out_name, 'w') as outfile:
        json.dump({"name":name,
                   "entire_sequence": entire_sequence, 
                   "pairs": pairs, 
                   "pair_types": pair_types, 
                   "segments":segments}, outfile, indent=True, sort_keys=True)            


if __name__ == "main":
    import argparse    
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument("-j","--json",default='out/data.json',help='specify path to json file')
    parser.add_argument('-i','--initial-model', help='Initial model in pdb format')
    parser.add_argument('-a','--alnfile', help='secify path to alignemnt file')
