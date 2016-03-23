
# -*- coding: utf-8 -*-
"""
Given a JSON file generate scripts for coloring the representation in Chimera and VMD.
"""
from __future__ import print_function
import utils as u
import yaml
import six
import re

def chimera_color(json_file, out_file, model_number="", verbose=True):
    d=u.load_json_data("pyr-p1b-SN-N5.json")

    with open(out_file+'.chimera', 'w') as f:        
        for s in d.segments:
            template = """color {color} #{model}:{start}-{end}"""
            s['model']=''
            line = template.format(**s)
            if verbose:
                print(line)
            f.write(line)


if __name__ == "__main__":
    import argparse    
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument("-j","--json", default='data.json', help='specify path to json file')
    parser.add_argument("-o","--out-file", default='color', help='Output name (without extension) of the coloring script')
    parser.add_argument("-v","--verbose", default=True, help='print scripts to screen')
    

    a = parser.parse_args()
    chimera_color(a.json, a.out_file, verbose=a.verbose)
