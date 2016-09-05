#! /usr/bin/python
"""Script for calculating approximate volume of coiled-coil polyhedra

Requires: SciPy v0.17.0
"""
from __future__ import print_function, division, absolute_import
import cocopod.utils as u
from scipy.spatial import ConvexHull
import numpy as np
import mdtraj as md
import argparse

parser = argparse.ArgumentParser(description="Script for calculating approximate volume of the polyhedra")
parser.add_argument('-m', '--model', help='Input file name (pdb)', type=str)
parser.add_argument('-j', '--json', help='JSON file defining segments configuration', type=str)
parser.add_argument('-o', '--out', help='Name of the outputted Bild file', type=str)
args = parser.parse_args()

#read input files
d = u.load_json_data(args.json)
md_obj = md.load(args.model)

#read topology and position
topology = md_obj.topology
positions = md_obj.xyz 
#converting positions from 3D np array into a nested list
positions = np.reshape(positions, (topology.n_atoms, 3))
positions = positions.tolist()

points = [] #points that will be used for ConvexHull calculation
for pair in d.pairs:  #go through all CC pairs
    pair_1_id, pair_2_id, pdbname = u.find_pair(pair['pair'][0], d.segments)
    p1_start = d.segments[pair_1_id]['start'] 
    p1_end = d.segments[pair_1_id]['end']
    p2_start = d.segments[pair_2_id]['start'] 
    p2_end = d.segments[pair_2_id]['end']
    #select CA atoms in the starting and ending residues of both pair segments 
    atom_1 = positions[topology.select("residue " + str(p1_start) + " and name CA")]
    atom_2 = positions[topology.select("residue " + str(p1_end) + " and name CA")]
    atom_3 = positions[topology.select("residue " + str(p2_start) + " and name CA")]
    atom_4 = positions[topology.select("residue " + str(p2_end) + " and name CA")]
    #from atom positions specify the aproximate coordinates of the start and end of the CC pair
    if d.pair_types[d.pairs.index(pair)] == 'A':
            avg1 = [(atom_1[i] + atom_4[i])/2 for i in range(len(atom_1))]
            avg2 = [(atom_2[i] + atom_3[i])/2 for i in range(len(atom_1))]
    else:
        avg1 = [(atom_1[i] + atom_3[i])/2 for i in range(len(atom_1))]
        avg2 = [(atom_2[i] + atom_4[i])/2 for i in range(len(atom_1))]
    #add the coordinates of the start and end of the CC pair to list of points
    points.append(avg1)
    points.append(avg2)
points = np.array(points)
hull = ConvexHull(points, qhull_options='Qt')
print("{} nm³".format(hull.volume))
#for simplex in hull.simplices:
#    print(points[simplex, 0], points[simplex, 1], points[simplex, 2])     
f = open(args.out, 'w')
f.write(".color red\n")
for simplex in hull.vertices:
    f.write(".dot {:6.3f} {:6.3f} {:6.3f}\n".format(points[simplex, 0]*10, points[simplex, 1]*10, points[simplex, 2]*10))   
