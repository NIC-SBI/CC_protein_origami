#!/usr/bin/env python

"""
Constructs an initial rough polyhedral model from JSON config file.

Json, alignment files and an initial model of the protein should be prepared beforehand. Initial model should be the protein
in the form of an alpha helix. Initially the program treats segments as rigid bodies. 
Several md runs are carried out in combination with harmonic restraints to bring the protein in roughly correct polyhedral shape.
"""
from __future__ import print_function, division, absolute_import
import ppmod.utils as u
import argparse
import os
from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, quasi_newton, actions
from modeller.automodel import *

parser=argparse.ArgumentParser('Program for making initial homology models of protein polyhedra.')
parser.add_argument("-j","--json",default='out/data.json',help='specify path to json file')
parser.add_argument("-hx","--helix",help='specify path to pdb file with protein as alpha helix')
parser.add_argument('-ms','--mdsteps',help='specify number of md steps',default=50000)
parser.add_argument('-m','--mean',help='specify desired distance between segments',default=9.0)
parser.add_argument('-sd','--stdev',help='specify desired standard deviation (strength of harominc constraint) of distance between segments',default=1.0)
parser.add_argument('-fri','--friction',help='set friction factor', default=0.0)
parser.add_argument('-t','--temp',help='set temperature during md run', default = 300)
parser.add_argument('-s','--shift',help='set limit of atomic shifts along an axis',default=0.39)
parser.add_argument('--time-step',help='set time step size in md runs (in fs)', default=4.0)
parser.add_argument('-o','--out-dir',help='output directory. Default is name+random-ppostfix', default=None, type=str)
parser.add_argument('-r','--rand-seed',help='Random seed. For modeler it mist be between -2 and -50000. If none, a random number will be chosen', 
                                       default=None, type=int)


args = parser.parse_args()

if args.rand_seed is None:
    import random
    args.rand_seed = random.randint(-50000, -2) 

d = u.load_json_data(args.json)
env = environ(rand_seed=args.rand_seed)
os.getcwd()
log.verbose()    # request verbose output

if args.out_dir is None:
    args.out_dir=d.name + u.id_generator()

# Read parameters (needed to build models from internal coordinates)
env.libs.topology.read('${LIB}/top_heav.lib') 
env.libs.parameters.read('${LIB}/par.lib')

env.io.atom_files_directory = ['./out','./building_blocks'] #where to read atom files
env.edat.dynamic_sphere = True

# read model file
mdl = complete_pdb(env, args.helix) 

#write the psf topology
mdl.write_psf(args.out_dir+'/00.psf')

# Select all atoms:
atmsel = selection(mdl)
at = mdl.chains[0].atoms
# Generate the restraints:
mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
chain = ":A"   #name of chain

#make rigid bodies
for seg in d.segments:
    print(seg.name, seg.start, seg.end)
    r = rigid_body(mdl.residue_range(str(seg.start)+chain, str(seg.end)+chain))
    mdl.restraints.rigid_bodies.append(r)


#--------------------------------------------------------------------------------
for i in range(len(d.pairs)):  #pairs are formed in the same order as written in .json file
    out_name = args.out_dir+'/01-{:02d}-{}-{}-initial'.format((i+1),d.pairs[i]['pair'][0], d.pairs[i]['pair'][1]) 
    
    for n in range(len(d.segments)):
        if d.segments[n]['name'] == d.pairs[i]['pair'][0]:
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

    print("PAIR: {} and {}".format(d.segments[p1].name, d.segments[p2].name),d.pair_types[i])
    print(p1_start,p1_end)
    print(p2_start,p2_end)
    
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
    trcfil = open(out_name+'.ener', 'w')
    cg.optimize(atmsel, max_iterations=100, actions=[actions.trace(10, trcfil)])       
    md.optimize(atmsel, max_iterations=run_md_steps, friction=float(args.friction), temperature=float(args.temp), init_velocities=True,            
            cap_atom_shift=float(args.shift), md_time_step=float(args.time_step),
            #guide_time=40, guide_factor=1,
            actions=[actions.charmm_trajectory(100, filename=out_name+'.dcd'),
                     actions.trace(100, trcfil)])
    mpdf = atmsel.energy()
    mdl.write(file=out_name+'.pdb')


#----------------------------------------------------
#one final optimization run
out_name=args.out_dir+"/02-final-initial-min"

trcfil = open(out_name+'.ener', 'w')
md.optimize(atmsel, max_iterations=run_md_steps, friction=float(args.friction), temperature=float(args.temp), init_velocities = True,            
            cap_atom_shift=float(args.shift), md_time_step=float(args.time_step),
            #guide_time=40, guide_factor=1,
            actions=[actions.charmm_trajectory(100, filename=out_name+'-md.dcd'),
                     actions.trace(100, trcfil)])

cg.optimize(atmsel, max_iterations=500, actions=[actions.trace(10, trcfil)])


#move the model to the origin
mdl.orient()
mdl.write(file=out_name+'.pdb')
#write the psf topology (it's the same as 00.psf)
mdl.write_psf(file=out_name+'.psf')
