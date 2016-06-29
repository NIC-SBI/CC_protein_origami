#!/usr/bin/env python

"""
Constructs a homology refined model from JSON config file and an initial polyhedral peptide model (see make_initial_model.py).

The raw model is then polished by homology modeling. In this step structures of individual segment pairs are used
to generate the final model. Multiple models are generated and evaluated according to their DOPE score and modeller
objective function. It is assumed that the model with the lowest DOPE score is the best.
"""
from __future__ import print_function
import argparse
import os
import utils as u
from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, quasi_newton, actions
from modeller.automodel import *

parser=argparse.ArgumentParser(__doc__)
parser.add_argument("-j","--json",default='out/data.json',help='specify path to json file')
parser.add_argument('-i','--initial-model', help='Initial model in pdb format')
parser.add_argument('-a','--alnfile', help='secify path to alignemnt file')
parser.add_argument('-mvi','--max-var-iterations', help='set max number of iterations for conjugate gradiens optimization method', default=500)
parser.add_argument('-nr','--repeat', help='number of optimization repeats', default=1)
parser.add_argument('-md','--md-level', help='How much MD optimization is performed', default = 'slow', 
                 choices=['very_fast', 'fast', 'slow', 'very_slow', 'slow_large'])
parser.add_argument('-smi','--start-index', help='set starting model index',default=1)
parser.add_argument('-emi','--end-index', help='set ending model index',default=1)
parser.add_argument('-o','--out-dir',help='output directory. Default is name+random-ppostfix', default=None, type=str)
parser.add_argument('-r','--rand-seed',help='Random seed. For modeler it mist be between -2 and -50000. If none, a random number will be chosen', 
                                       default=None, type=int)

args = parser.parse_args()

if args.rand_seed is None:
    import random
    args.rand_seed = random.randint(-50000, -2) 

d = u.load_json_data(args.json)
env = environ(rand_seed=args.rand_seed)

log.verbose()    # request verbose output

if args.out_dir is None:
    args.out_dir=d.name + u.id_generator()

if not os.path.exists(args.out_dir):
    os.makedirs(args.out_dir)

#convert to absolute path before chdir
args.alnfile =  os.path.abspath(args.alnfile)

os.chdir(args.out_dir)

print("################################")
print("Blocks_file after chdir", os.path.abspath(u.relative_to(__file__,'../../building_blocks')))
    
# Read parameters (needed to build models from internal coordinates)
env.libs.topology.read('${LIB}/top_heav.lib') 
env.libs.parameters.read('${LIB}/par.lib')


env.io.atom_files_directory = ['.','..','building_blocks', 
                               os.path.abspath(u.relative_to(__file__,'../../building_blocks'))] #where to read atom files
env.edat.dynamic_sphere = True

md_opt_dict = \
    {'very_fast'  : refine.very_fast,
     'fast'       : refine.fast,
     'slow'       : refine.slow,
     'very_slow'  : refine.very_slow,
     'slow_large' : refine.slow_large,
    }

#print dir(automodel.get_model_filename)

import shutil

class AlphaModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
        #add secondary alpha-helical restraints for each segement
        for seg in d.segments:
            #print seg.name, seg.start, seg.end
            rsr.add(secondary_structure.alpha(self.residue_range(str(seg.start), str(seg.end))))

    def get_model_filename(self, sequence, id1, id2, file_ext):
        print(self, sequence, id1, id2, file_ext)

        return "03-homology-model-{id2:02}{file_ext}".format(
            outdir=args.out_dir, id1=id1, id2=id2, file_ext=file_ext)

#determine knowns and target sequence            
sequence, knowns = u.sequnce_and_knowns(args.alnfile)
                    
a = AlphaModel(env,
              alnfile=args.alnfile, 
              knowns=knowns,     
              sequence=sequence,        
              inifile=args.initial_model,
              assess_methods=[assess.DOPE, assess.GA341, assess.normalized_dope])              

a.max_var_iterations = int(args.max_var_iterations)         
a.md_level = md_opt_dict[args.md_level]
#a.md_level = refine.slow_large     
a.repeat_optimization = int(args.repeat)
a.initial_malign3d = True

a.starting_model = int(args.start_index)                 
a.ending_model = int(args.end_index)
                                   
a.make()                           
