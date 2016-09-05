# -*- coding: utf-8 -*-
"""
Given a JSON file generate scripts for coloring the representation in Chimera and VMD.
"""
from __future__ import print_function, division, absolute_import
import cocopod.utils as u

def print_and_write(line, file, verbose=True):
    """prints to std out and also writes to file. Only printed if verbose us True"""
    if verbose:
        print(line)
    file.write(line+"\n")

 

def chimera_color(json_file, out_file, model_number="", verbose=True, color_map=None, 
                  default_color="#d2d2b4b48c8c", add_surface=True, thick_vdw_size=2.5, default_vdw_size=1.7, print_thick_vdw=False):
    """Print scipt for coloring in chimera"""

    d=u.load_json_data(json_file)
    
    with open(out_file+'.chimera', 'w') as f:  
    
        def print_segment(s, vdw_size=None):
            template = """color {color} #{model}:{start}-{end}"""
            s['model']=model_number
            s['vdw_size']=vdw_size
            #override the color of segments if map is present
            if not color_map is None:
                s['color']=color_map.get(s.name, s.get('color', default_color))
            else:
                s['color']=s.get('color', default_color)

            line = template.format(**s)
            print_and_write("#"+s.name, f, verbose)
            print_and_write(line, f, verbose)       

            if vdw_size:
                #define the VDW scale
                line = """vdwdefine {vdw_size} #{model}:{start}-{end}""".format(**s)
                print_and_write(line, f, verbose) 



        #color the model with the default color (the color for the linkers)
        print_and_write("#Default model color and vdw size", f, verbose)
        line = "color {color} #{model}".format(color=default_color,model=model_number)
        print_and_write(line, f, verbose)
        line = "vdwdefine {vdw_size} #{model}".format(vdw_size=default_vdw_size,model=model_number)
        print_and_write(line, f, verbose)

        #group the commands by pairs of segments
        printed_segments_ids = []
        for s in d.segments:
            if not s.id in printed_segments_ids:
                print_segment(s)    
                print_segment(d.segments[s.pair_id-1])
                printed_segments_ids.extend([s.id, s.pair_id])

        if add_surface:    
            print_and_write("", f, verbose)
            line = "surface probeRadius 3   vertexDensity 2 allComponents false  protein"
            print_and_write(line, f, verbose)
            line = "surface probeRadius 3   vertexDensity 2 allComponents false  protein"
            line = "transparency 80,s #{model}:".format(model=model_number)
            print_and_write(line, f, verbose)

        if print_thick_vdw:
            print_and_write("", f, verbose)
            print_and_write("#Thicker segments", f, verbose)
            printed_segments_ids = []
            for s in d.segments:
                if not s.id in printed_segments_ids:
                    print_segment(s, thick_vdw_size)    
                    print_segment(d.segments[s.pair_id-1], thick_vdw_size)
                    printed_segments_ids.extend([s.id, s.pair_id])            

if __name__ == "__main__":


    import argparse    
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument("-j","--json", default='data.json', help='specify path to json file')
    parser.add_argument("-o","--out-file", default=None, help='Output name (without extension) of the coloring script. Default taken from json')
    parser.add_argument("-v","--verbose", default=True, help='print scripts to screen')
    parser.add_argument("-m","--color-map", default=None, help='colormap to override colors given in the json. Works even if no colors are present in the Json')
    parser.add_argument("-t","--thick_vdw", default=True, help='Print specification for thicker segments (for movie creation).')

    a = parser.parse_args()
    
    if a.out_file is None:
        a.out_file = a.json.replace('.json','')

    if not a.color_map is None:
        a.color_map = u.load_json_data(a.color_map)
    chimera_color(a.json, a.out_file, verbose=a.verbose, color_map=a.color_map, print_thick_vdw=a.thick_vdw)

    

 


