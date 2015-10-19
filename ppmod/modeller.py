#!/usr/bin/env python

"""
Constructs a homology model for different TET variants from sequence

Json, alignment files and an initial model of the protein should be prepared beforehand. Initial model should be the protein
in the form of an alpha helix (easily achieved by chimera). Initially the program treats segments as rigid bodies. 
Several md runs are carried out in combination with harmonic restraints to bring the protein in tetrahedron shape.
The raw model is then polished by homology modeling. In this step structures of individual segment pairs are used
to generate the final model. Multiple models are generated and evaluated according to their DOPE score and modeller
objective function. It is assumed that the model with the lowest DOPE score is the best.
"""
import argparse
import os
import utils as u
from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, quasi_newton, actions
from modeller.automodel import *

parser=argparse.ArgumentParser('Program for making homology models for TET.')
parser.add_argument("-j","--json",default='out/data.json',help='specify path to json file')
parser.add_argument("-hx","--helix",help='specify path to pdb file with protein as alpha helix')
parser.add_argument('-ms','--mdsteps',help='specify number of md steps',default=50000)
parser.add_argument('-m','--mean',help='specify desired distance between segments',default=9.0)
parser.add_argument('-sd','--stdev',help='specify desired standard deviation of distance between segments',default=1.0)
parser.add_argument('-aln','--alnfile', help='secify path to alignemnt file')
parser.add_argument('-mvi','--max_var_iterations', help='set max number of iterations for conjugate gradiens optimization method', default=100)
parser.add_argument('-nr','--repeat', help='number of optimization repeats', default=10)
parser.add_argument('-smi','--startindex', help='set starting model index',default=1)
parser.add_argument('-emi','--endindex', help='set ending model index',default=10)
parser.add_argument('-fri','--friction',help='set friction factor', default=0.0)
parser.add_argument('-t','--temp',help='set temperature during md run', default = 300)
parser.add_argument('-s','--shift',help='set limit of atomic shifts along an axis',default=0.39)
parser.add_argument('--time',help='set time step size in md runs (in fs)', default=4.0)

args = parser.parse_args()
d = u.load_json_data(args.json)
env = environ()
os.getcwd()
log.verbose()    # request verbose output

# Read parameters (needed to build models from internal coordinates)
env.libs.topology.read('${LIB}/top_heav.lib') 
env.libs.parameters.read('${LIB}/par.lib')

env.io.atom_files_directory = ['./out','./building_blocks'] #where to read atom files
env.edat.dynamic_sphere = True

# read model file
mdl = complete_pdb(env, args.helix) 

# Select all atoms:
atmsel = selection(mdl)
at = mdl.chains[0].atoms
# Generate the restraints:
mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
chain = ":A"   #namechain
#make rigid bodies
for seg in d.segments:
    print seg.name, seg.start, seg.end
    r = rigid_body(mdl.residue_range(str(seg.start)+chain, str(seg.end)+chain))
    mdl.restraints.rigid_bodies.append(r)
#%%time
#--------------------------------------------------------------------------------
for i in range(len(d.pairs)):  #pairs are formed in the same order as written in .json file
    out_name = 'out/{:02d}-{:.3s}'.format((i+1),d.pairs[i][0]) 
    
    for n in range(len(d.segments)):
        if d.segments[n]['name'] == d.pairs[i][0]:
            p1 = d.segments[n]['id']
            p2 = d.segments[n]['pair_id']
            break
    run_md_steps = int(args.mdsteps)    
             
#convert to 0 based indices
    p1 = p1 -1
    p2 = p2 -1
    p1_start = str(d.segments[p1].start)
    p1_end = str(d.segments[p1].end)
    p2_start = str(d.segments[p2].start)
    p2_end = str(d.segments[p2].end)

    print "PAIR: {} and {}".format(d.segments[p1].name, d.segments[p2].name),d.pair_types[i]
    print p1_start,p1_end
    print p2_start,p2_end
    
#add restrain between segments belonging to the same pair
    if d.pair_types[i] == 'A':
        mdl.restraints.add(forms.gaussian(group=physical.xy_distance,     
                               feature=features.distance(at['CA:'+str(p1_start)],
                                                        at['CA:'+str(p2_end)]),
                               mean=float(args.mean), stdev=float(args.stdev)))            
        mdl.restraints.add(forms.gaussian(group=physical.xy_distance,
                               feature=features.distance(at['CA:'+str(p1_end)],
                                                        at['CA:'+str(p2_start)]),
                               mean=float(args.mean), stdev=float(args.stdev))) 
    else: 
        mdl.restraints.add(forms.gaussian(group=physical.xy_distance,     
                               feature=features.distance(at['CA:'+str(p1_start)],
                                                        at['CA:'+str(p2_start)]),
                               mean=float(args.mean), stdev=float(args.stdev)))            
        mdl.restraints.add(forms.gaussian(group=physical.xy_distance,
                               feature=features.distance(at['CA:'+str(p1_end)],
                                                        at['CA:'+str(p2_end)]),
                               mean=float(args.mean), stdev=float(args.stdev)))    

# Create optimizer objects and set defaults for all further optimizations
    cg = conjugate_gradients(output='REPORT')  
    md = molecular_dynamics(output='REPORT')

# Open a file to get basic stats on each optimization
    trcfil = open(out_name+'.log', 'w')
    cg.optimize(atmsel, max_iterations=100, actions=[actions.trace(10, trcfil)])       
    md.optimize(atmsel, max_iterations=run_md_steps, friction=float(args.friction), temperature=float(args.temp), init_velocities=True,            
            cap_atom_shift=float(args.shift), md_time_step=float(args.time),
            #guide_time=40, guide_factor=1,
            actions=[actions.charmm_trajectory(100, filename=out_name+'-md.dcd'),
                     actions.trace(10, trcfil)])
    mdl.write(file=out_name+'-md.pdb')
    mpdf = atmsel.energy()
    mdl.write(file=out_name+'.pdb')


#----------------------------------------------------
#one final optimization run
out_name="out/07a-final-min"

trcfil = open(out_name+'.log', 'w')
md.optimize(atmsel, max_iterations=run_md_steps, friction=float(args.friction), temperature=float(args.temp), init_velocities = True,            
            cap_atom_shift=float(args.shift), md_time_step=float(args.time),
            #guide_time=40, guide_factor=1,
            actions=[actions.charmm_trajectory(100, filename=out_name+'-md.dcd'),
                     actions.trace(10, trcfil)])
mdl.write(file=out_name+'-md.pdb')
cg.optimize(atmsel, max_iterations=500, actions=[actions.trace(10, trcfil)])

# do the actual homology modelling
mdl.write(file=out_name+'.pdb')

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
        for seg in d.segments:
            print seg.name, seg.start, seg.end
            rsr.add(secondary_structure.alpha(self.residue_range(str(seg.start), str(seg.end))))

#determine knowns and target sequence            
knowns = ()
with open(args.alnfile,'r') as f:
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
                    
a = MyModel(env,
              alnfile=args.alnfile, 
              knowns=knowns,     
              sequence=sequence,        
	      inifile='07a-final-min.pdb',
              assess_methods=assess.DOPE)              

a.max_var_iterations = int(args.max_var_iterations)         
a.md_level = refine.slow 	
#a.md_level = refine.slow_large 	
a.repeat_optimization = int(args.repeat)
a.initial_malign3d = True
a.get_refine_actions()
a.starting_model = int(args.startindex)                 
a.ending_model = int(args.endindex)
                                   
a.make()                           