#taken from http://hayd.github.io/2013/dotable-dictionaries/
class Dotable(dict):
    """Generates nested dotable dicts from a json-like object. This makes is possible to write
    a.b[0].c
    """
    __getattr__= dict.__getitem__

    def __init__(self, d):
        self.update(**dict((k, self.parse(v))
                           for k, v in d.iteritems()))

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
        return {byteify(key):byteify(value) for key,value in input.iteritems()}
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
    data=byteify(data)
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
#taken from http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))    


def sequnce_and_knowns(alnfile):
    """Parses the IDs of the known strucutes and the sequnces from the alignment file.
    Returns a tuple:
      (seq, (known1, known2 ...))

    """
    knowns = ()
    with open(alnfile,'r') as f:
        while True:
            line = f.readline()
            if not line: break
            if line[:4] == '>P1;':
                seqname = line[4:].rstrip('\n')
                line = f.readline()            
                if line.split(':')[0] == 'sequence':
                    sequence = seqname
                else:
                     knowns = knowns+(seqname,)
    return (sequnce, knowns)                 
