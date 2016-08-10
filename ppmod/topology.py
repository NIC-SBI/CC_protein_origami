#!/usr/bin/python
from __future__ import print_function, division

#####################################
###Ajasja Ljubetic (ajasja.ljubetic@gmail.com) ##################
###utilities for dealing with topologies
from __future__ import print_function, division, absolute_import
import ppmod.utils as u

import collections

def permute_segment_left(segs, n=1):
    """Permutes a segment to the left for n movements"""
    q = collections.deque(segs)
    q.rotate(-n)
    return list(q)


def permute_vertices_left(verts, n=1):
    """Permutes a segment to the left for n movements"""
    q = collections.deque(verts)
    #last vertex is same as first
    last = q.pop()
    q.rotate(-n)
    q.append(q[0])
    return list(q)

def get_permutation_name(basename, permutation=1):
    """gets the name for the permutation with a, b, c...
    Currently does not work for negative"""
    return str(basename)+"."+str(permutation)    


def get_segment_distances_dict(segment_topology):
    """given an array of char as segments find the distance between pairs in the sequnce.
    Antipralalel are given a shorter distance than parallel."""
    stl = [c.upper() for c in segment_topology]
    #delete duplicates
    blocks = sorted(list(set(stl)))
    
    dist = collections.OrderedDict()
    for block in blocks:
        inds = [i for i, x in enumerate(segment_topology) if x.upper() == block]
        #print(inds)
        assert len(inds) == 2, block+ " does not have a pair: "+ str(len(inds))
        
        #if they are difftent case they are antiprallel
        #print(segment_topology[inds[0]],segment_topology[inds[1]]) 
        if segment_topology[inds[0]]!=segment_topology[inds[1]]:
            #antiparallel are shorter distance
            dist[block]=(abs(inds[1]-inds[0]-1))
        else:
            dist[block]=(abs(inds[1]-inds[0]))
    return dist    

def get_segment_distances(segment_topology):
    """given an array of char as segments find the distance between pairs in the sequnce.
    Antipralalel are given a shorter distance than parallel."""
    return list(get_segment_distances_dict(segment_topology).values())


def name_topologies_and_permutations(top_dataframe):
    """Given a dataframe of topologies, enumerate all topologies"""
    import pandas
    import numpy as np
    allpnames = []
    for n, (i, row) in enumerate(top_dataframe.iterrows()):
        
        for p in range(len(row['topo'])):
            allpnames.append(get_permutation_name(n+1,p+1))
    
    df = pandas.DataFrame(index=allpnames,columns="segments num_AP num_P num_cross reflected".split())
    for n, (i, row) in enumerate(top_dataframe.iterrows()):
        s = list(row['topo'])
        #print(s)
        for p in range(len(s)):
            name = get_permutation_name(n+1,p+1)
            ps = permute_segment_left(s, p)

            df.loc[name, 'segments'] = "".join(ps)
            df.loc[name, 'num_AP']    = row['num_AP']
            df.loc[name, 'num_P']     = row['num_P']
            df.loc[name, 'num_cross'] = row['num_cross']
            df.loc[name, 'reflected'] = row['reflected']
    return df

def calculate_TCO(top_dataframe):
    """Given a dataframe of topologies, enumerate all topologies"""
    import pandas
    import numpy as np
    
    df = pandas.DataFrame(index=top_dataframe.index,columns="min max TCO stdTCO segments num_AP num_P num_cross ref".split())
    for n, (name, row) in enumerate(top_dataframe.iterrows()):
        s = list(row['segments'])
        dist = get_segment_distances(s)
        
        df.loc[name, 'min'] = np.min(dist)
        df.loc[name, 'max'] = np.max(dist)
        df.loc[name, 'TCO'] = np.mean(dist)
        df.loc[name, 'stdTCO'] = np.std(dist)
        df.loc[name, 'segments'] = row['segments']
        df.loc[name, 'num_AP']    = row['num_AP']
        df.loc[name, 'num_P']     = row['num_P']
        df.loc[name, 'num_cross'] = row['num_cross']
        df.loc[name, 'ref'] = row['reflected']
    return df    


def name_of_topology(top, top_list, verbose=False):
    """Given a list of topologies (as dict or dataframe), return the "name" of the topology from the list"""
    res = []
    s_top = standard(top)
    
    for index, row in top_list.iterrows():
        #print(index)
        s_form = standard(row['segments'])
        if verbose:
            print (s_form, s_top)
        if s_form  == s_top:
            res.append(index)
    
    #check for reflected
    for index, row in top_list.iterrows():
        s_form =  standard(row['segments'][::-1])           
        if s_form  == s_top:
            res.append(str(index)+'R')


    
    assert len(res) >= 1, str(res)
    #due to symmetry there might be more than one solution. Just return the first.
    if len(res) > 1:
        print(res)
    return res[0]    
#####################################
#####################################

########################################
#Load ply files
def load_vfaces(ply_file):
    """Loads a list of vFaces from a ply file."""
    import plyfile
    plydata = plyfile.PlyData.read(ply_file)
    #print(plydata)
    #convert to list of lists. List of numpy arrays does not work.
    #Trouble with vertex_index vs vertex_indicies. Just take the first property     
    vfaces = [list(f[0]) for f in plydata['face']]
    return vfaces
    

def convert_vface_to_efaces(vfaces):
   """Convert faces given by vertices to faces given by strings (edge faces).
It works for up to 26 edges."""
   #list of edges   
   eds = []
   nA = ord("A")
   na = ord("a")
   faces = []
   for ff in vfaces:
      #print("ff",ff) 
      f = ff+[ff[0]]
      n = len(f)
      e = ""
      for i in range(n-1):
         aa = f[i]
         bb = f[i+1]
         a,b = aa,bb 
         if a > b:
            a,b = b,a
         if [a,b] not in eds:
            eds.append([a,b])
         i = eds.index([a,b])
         if aa > bb:
            e += chr(na+i)
         else:
            e += chr(nA+i)
      faces.append(e)
   return faces
############################################         
         

# Some additions for the Biomathematics afternoon in Ljubljana 2016
# February 17,2016.
#bipyramid = ["ahF","Abc","egf","giB","dic","eDh"]
bipyramid = ["ahf","Abc","egF","gbI","dic","eDh"]
triprism = ["abc","def","bifH","cgeI","ahdG"] #Triangular prism
tetrahedron = ["abc","adF","beD","cfE"] #Tetrahedron in the sphere
tetra = tetrahedron
bipyramid = ["adF","beD","cfE","agI","bhG","ciH"] #Double Tertahedron (triangular bipyramid in the sphere)
octahedron = ["aeF","bfG","cgH","dhE","aiJ","bjK","ckL","dlI"] #Octahedron
square_pyramid = ["aeF","bfG","cgH","dhE","ABCD"] #Pyramid
pyramid = square_pyramid
doubletetra = bipyramid 
cube = ["abcd","ijkl","aeIH","bfJE","cgKF","dhLG"]

# Modified for all embeddings.
#
# Slight improvements. Ferburay 1, 2016.
#
# We are given a polyhedron, compute its graph(1-skeleton)
# and compute all local rotations.
# for each local rotation we compute back a set of faces.
# T.P. January 29, 2016.

from time import *
from itertools import permutations

def explore(poly,verbose=False, return_dataframe = True):
    """Returns all topologies of a given polyhedron. 
    If returned as a Pandas dataframe it is sorted by 'num_cross', 'num_AP', 'topo', 'reflected'"""
    #a list of hashes returned as result
    top_data = []
    (sup,crs) = represent(poly,verbose)
    for i in range(len(sup)):
        if verbose:
            print("id,num_AP,num_P,num_cross,topo,reflected")                
        else:
            print("id,num_cross, super canonical")
            print(i,crs[i],sup[i],others(sup[i]))       
        same_crossings([sup[i]], verbose=verbose, crossings=crs[i], top_data=top_data)
    if return_dataframe:
        import pandas
        top_data = pandas.DataFrame(data=top_data)
        top_data.sort_values(by=['num_cross', 'num_AP', 'topo', 'reflected'], inplace=True)

    return top_data

def represent(poly,verbose = False):
    """should return all crossing non-equvalent representatives of a polyhedron."""
    tt = super_canonical(makesingle(poly))
    if not verbose:
        print("super canonical",tt,others(tt))
    ae = all_embeddings(tt)
    ttl = [canonical(makesingle(stringsfromfaces(facesfromrotation(r)))) for r in ae]
    tts = list(set(ttl))
    if not verbose:
        print("ttl",len(ttl),"tts",len(tts))
    canons = []
    sup= []
    rots = []
    k = 0
    for n,ctt in enumerate(tts):
        if ctt in canons:
            continue
        newcan = generate_canonicals(ctt)
        supercan= min(newcan)
        if not verbose:
            print(n,k, len(newcan),supercan)
            if supercan == tt:
                print("------------")
            print(ae[n])
        k += 1
        sup.append(min(newcan))
        rots.append(ae[n])
        canons.extend(newcan)
    i0 = sup.index(tt)
    rot0 = rots[i0]
    crs = []
    for j,t in enumerate(sup):
        crs.append(compare(rot0,rots[j]))
    return (sup,crs)

def same_crossings(poly,verbose=False, crossings=None, top_data=None):
    """all topologies with the same crossings."""
    ga = generateall(poly)
    gar = [oriented_canonical(x[::-1]) for x in ga]
    dic = {}
    n = len(ga)
    k = 0
    for i in range(n):
        (na,np) = noparallel(ga[i])
        dic[na] = dic.get(na,0)+1 
        result_dict = {"num_AP":na, "num_P":np, "num_cross":crossings,
                        "topo":ga[i], "reflected":ga[i] == gar[i]}

        top_data.append(result_dict)
        if ga[i] == gar[i]:
            if verbose:
                print(i,na,np,crossings,ga[i], "", sep=",")
        else:
            if verbose:
                print(i,na,np,crossings,ga[i], "", sep=",")
                print(i,na,np,crossings,gar[i],"R", sep=",")
            k += 1
    if not verbose:
        print("Number of topologies ", n)
        print("Number of reflexive topologies ",n-k)
        print("Number of non-reflexive pairs ",k)
        print()
        print("Number of topologies by antiparallel dimers. \n (num AP, num topo)")
        for k in dic:
            print(k,dic[k])
    return (n,k,dic,ga)        



def compare(rot0,rot):
    tot = 0
    for x in rot0:
        if rot0[x] != rot[x]:
              tot += 1
    return tot

# T.P. February 17, 2016.
def representatives(poly):
    """should return all crossing non-equvalent representatives of a polyhedron."""
    tt = canonical(makesingle(poly))
    print("**")
    ae = all_embeddings(tt)
    print("ae ")
    ff = [facesfromrotation(r) for r in ae]
    print("ff")
    tts = [stringsfromfaces(r) for r in ff]
    print("tts",len(tts),tts)
    sstrs = [super_canonical(makesingle(x)) for x in tts]
    return sstrs

def all_embeddings(tt):
    """should return all rotations for a given strand."""
    edges = skeletonedges(tt)
#    print("skeleton edges",edges)
    rot = all_rotations(edges)
    return rot

def skeletonedges(tt):
#    print("skeleton edges input tt :",tt)
    skel = skeleton(tt)
#    print("skel ",skel)
    edgs = skel[2]
    verts = skel[1]
    newv = {verts[i]:i for i in range(len(verts))}
    edges = []
    for (xx,yy) in edgs:
        x = newv[xx]
        y = newv[yy]
        if x > y:
            x,y = y,x
        if not (x,y) in edges:
            edges.append((x,y))
            
    return edges    

# TODO: DEBUG HERE!
def skeleton(tt):
    """Determine the skeleton of the fundamental polygon."""
    (gp,cn,prt) = partition(tt)
##    print("gp = ",gp)
##    print("cn = ",cn)
##    print("prt = ",prt)
    edgs = edges(tt)
    # print(tt)
    # print(edgs)
    vts = []
    for i in  cn:
        if i == cn[i]:
            vts.append(i)
    nedgs = []
    for (u,v) in edgs:
        nedgs.append((cn[u],cn[v]))
##    print("skeleton output :")
##    print("sk 0 :", cn)
##    print("sk 1 :", vts)
##    print("sk 2 :", nedgs)
##    print("sk 3 :", gp)
##    print("sk 4 :", prt)
    return (cn,vts,nedgs,gp,prt)

def partition(tt):
    """Partition of the set of vertices into equivalence classes"""
    n = len(tt)
    gp = gluev(tt)
    prt = {}
    for i in range(n):
        prt[i] = {i}
    goon = True
    while goon:
        goon = False
        for (i,j) in gp:
            np = prt[i]|prt[j]
            if len(np) > min(len(prt[i]),len(prt[j])):
                goon = True
            prt[i] = prt[j] = np
    pr = {}
    for i in range(n):
        pr[i] = min(prt[i])
    return (gp,pr,prt)

def gluev(tt):
    """ Glue vertices for each pair of letters."""
    n = len(tt)
#    print(n,tt)
    ttl = tt.lower()
    p = pairs(ttl)
    gp = []
    for (i,j) in p:
        if tt[i] == tt[j]:
            gp.append((i,j))
            gp.append(((i+1)%n,(j+1)%n))
        else:
            gp.append((i,(j+1)%n))
            gp.append(((i+1)%n,j))
    return gp

def pairs(t):
    """Collection of indices of pairs of letters (regardless of orientation)"""
    tt = t.lower()
    n = len(tt)
    p = []
    for i,x in enumerate(tt):
        j = tt.find(x)
        if j == i:
            continue
        p.append((j,i))
    return p    
        
def stringsfromfaces(faces):
    """Faces that are given in terms of closed walks using vertices,
are transformed into strings."""
    arcs = []
    for f in faces:
        x = f[0]
        for y in f[1:]:
            arcs.append((x,y))
            x = y
    arcs.sort()
    verts = []
    for (x,y) in arcs:
        if not x in verts:
            verts.append(x)
        if not y in verts:
            verts.append(y)
    verts.sort()
    edges = [(x,y) for (x,y) in arcs if x < y]
    edges.sort()
    syml = "a"
    symu = "A"
    sdic = {}
    if len(edges) > 26:
        print("alarm")
    for (x,y) in edges:
        sdic[(x,y)] = syml
        sdic[(y,x)] = symu
        syml = chr(ord(syml)+1)
        symu = chr(ord(symu)+1)
    sfaces = []
    for f in faces:
        sf = ""
        x = f[0]
        for y in f[1:]:
            sf += sdic[(x,y)]
            x = y
        sfaces.append(sf)
    return sfaces
        
def generateall(poly):
    """ Generate all non-intersecting polypeptides with a given polyhedron."""
    res = []
    some = makesingle(poly)
    ns = set(some.lower())
    nn = ""
    for x in ns:
        nn += x
    d = len(nn)
    D = 2**d
    for chi in range(D):
        samp = ""
        t = chi
        i = 0
        while chi > 0:
            if chi % 2 == 1:
                samp += nn[i]
            i += 1
            chi //= 2
        pol = glueorsplit(poly,samp)
        if len(pol) == 1:
            sx = canonical(pol[0])
            found = False
            for x in res:
                if isomorphic(x,sx):
                    found = True
    #                print(isomorphism(x,sx))
                    break
            if not found:
                res.append(sx)
    #            print(i)
    return res


def makesingle(llf,samp=""):
    """Given a polyhedron described by a list of faces, find a surface with one face.
first apply some gluesplitting along samp"""
    lf = glueorsplit(llf,samp)
    while len(lf)>1:
#        print(lf)
        lf = gluetwo(lf)
    return lf[0]

        







# T.P. January 31, 2016.
# Canonical re-labeling of a fundamental polygon with respect to dihadral equivalence.

def standard(tt):
    """Transform a string to an isomorphic that starts withn an 'a'
    without rotation or reflection. """
    su = "a"
    sl = "A"
    dic = {}
    res = ""
    for s in tt:
        if not s in dic:
            dic[s] = sl
            res = res+sl
            sl = chr(ord(sl)+1)
            dic[other(s)] = su
            su = chr(ord(su)+1)
        else:
            res = res + dic[s]
    return res
    
def all_equivalents(tt):
    """a list of all rotated and reflected strings. """
    ae = []
    ttr = tt[::-1]
    for i in range(len(tt)):
        ttt = tt[i:]+tt[:i]
        ae.append(ttt)
        tttr = ttr[i:]+ttr[:i]
        ae.append(tttr)
    return ae

def oriented_equivalents(tt):
    """a list of all rotated strings. """
    ae = []
    for i in range(len(tt)):
        ttt = tt[i:]+tt[:i]
        ae.append(ttt)
    return ae

#
# canonical can be speeded up by reversing the order of operation??
#

def generate_canonicals(tt):
    ga = generateall([tt])
#    print("ga = ",ga)
    ca = [canonical(x) for x in ga]
    return ca           
             


def super_canonical(tt):
    ga = generateall([tt])
#    print("ga = ",ga)
    ca = [canonical(x) for x in ga]
    return min(ca)           

             
             
def canonical(tt):
    """ canonical value of a string. """
    ae = all_equivalents(tt)
    aes = [standard(x) for x in ae]
    return min(aes)

def all_standard(tt):
    """all standard values of a string. """
    ae = all_equivalents(tt)
    aes = [standard(x) for x in ae]
    return aes


def oriented_canonical(tt):
    """oriented canonical value of a string. """
    ae = oriented_equivalents(tt)
    aes = [standard(x) for x in ae]
    return min(aes)

# T.P. January 29, 2016.

    
def facesfromrotation(rot):
    """ Typical orientable face-tracing algorithm (see Pisanski, Potocnik in
    Handbook of Graph Theory)."""
#    print("faces from rotation.")
    eds = []
    for x in rot:
        for y in rot[x]:
            if (x,y) not in eds:
                eds.append((x,y))
##    for x in rot:
##        print(x,rot[x])
    m = len(eds)
#    print ("m",m)
    faces = []
    start = True
    while len(eds) > 0:
        if start:
            (x,y) = eds[0]
#            print(x,y)
            eds.remove((x,y))
#            print(x,y)
            f = [x,y]
            xx,yy = x,y
            start = False
        i = rot[y].index(x)+1
        if i >= len(rot[y]):
            i = 0
        z = rot[y][i]
        x,y = y,z
        if xx == x and yy == y:
            start = True
            faces.append(f)
        else:
            eds.remove((x,y))
            f.append(z)
    faces.append(f)
    return faces
            
        
    

def rotation(edges):
    """ From ordered pairs to an ordered neighborhood
    dictionary. """
    rot = {}
#    print("rotation input edges ",edges)
    for (x,y) in edges:
        rot[x] = rot.get(x,[])
        rot[x].append(y)
        rot[y] = rot.get(y,[])
        rot[y].append(x)
    for z in rot:
        li = rot[z]
        if len(li) > 1:
            mi = min(li)
            i = li.index(mi)
            li = li[i:]+li[:i]
            if li[1]> li[-1]:
                lir = li[1:]
                lir = lir[::-1]
                li = li[0:1]+lir
            rot[z] = li
#    print ("rotation ouptut rot :",rot)
    return rot

def all_rotations(edges):
    """ Determine all rotations and the number of crosssings for a given graph. """
    rot = rotation(edges)
#    print("rot",rot)
    arot = [] 
    for x in rot:
        arot.append((x,allrot(rot[x])))
    alrot = allr(arot)
    dics = []
    for r in alrot:
        tmp = {}
        for (x,v) in r:
            tmp[x] = v
        dics.append(tmp)
    return dics

def allr(arot):
    if len(arot) == 0:
        return [[]]
    (x,rts) = arot[0]
    dics = allr(arot[1:])
    res = []
    for di in dics:
        for r in rts:
            tmp = [(x,r)]+di
            res.append(tmp[:])
    return res

def allrot(li):
#    li.sort()
    if len(li) < 4:
        return [(0,li)]
    gi = list(permutations(li[1:]))
    ti = [list(x) for x in gi]
    fi = [x for x in ti if x[0] < x[-1]]
    t = li[0]
    res = [[t]+x for x in fi]
    k = 1
    if li == res:
        k = 0
    return (k,res)   

def allrot(li):
    li.sort()
    if len(li) < 4:
        return [li]
    gi = list(permutations(li[1:]))
    ti = [list(x) for x in gi]
    fi = [x for x in ti if x[0] < x[-1]]
    t = li[0]
    res = [[t]+x for x in fi]
    return res   
        
        


# Three realizations of the tetrahedron.
tet1 = "abcadEcfEBdF"
tet2 = "abcadEBdFCeF"
tet3 = "adFabeFCBdEc"

# Experiments with self-assembling surfaces.
# T. Pisanski June 4,2013 

# Isomorphisms and automorphisms
# T. Pisanski May 20,2013.

def isomorphic(tt,ttt):
    """True if isomorphic."""
    (a,b) = isomorphism(tt,ttt)
    return (len(b) > 0)


def simpleiso(tt,ttt):
    """Simple isomorphism - change of letters only.""" 
    dic = {}
    for j,x in enumerate(tt):
#            print(j,x)
        dic[x] = dic.get(x,"")
        y = dic[x]
        if y == "":
            y = ttt[j]
            dic[x] = y
            yy = other(y)
            xx = other(x)
            dic[xx] = yy
        elif y == ttt[j]:
            continue
        else:
            return {}
    return dic

def isomorphism(tt,ttt):
    cyciso = cyclicisomorphism(tt,ttt)
    return cyciso

def dihedralisomorphism(tt,ttt):
    """Determine the dictionary of some dihedral isomorphism between tt and ttt"""
    dihiso = cyclicisomorphism(tt,ttt[::-1])
    return dihiso
        
        
def cyclicisomorphism(tt,ttt):
    """Determine the dictionary of some cyclic isomorphism between tt and ttt"""
    n = len(tt)
    if n != len(ttt):
#        print("Not the same length")
        return dic
    for i in range(n):
#        print(i)
        ttr = ttt[i:]+ttt[:i]
        dic = simpleiso(tt,ttr)
        if len(dic) > 0:
            return(i,dic) 
    return (-1,{})
            
def automorphisms(tt):
    """Return the list of shifts that will result in an automorphism of the string."""
    aut = []
    for i in range(len(tt)):
        ttt = tt[i:]+tt[:i]
        dic = simpleiso(tt,ttt)
        if len(dic) > 0:
#            print(dic)
            aut.append(i)
    return aut
                



# Single face embedding procedure
# T. Pisanski, May 19,2013
from random import *


smallpolyhedra = {"tetrahedron":tetra,"4-sided pyramid":pyramid,"triangular prism": triprism}
polyhedra = {"tetrahedron":tetra,"cube":cube,"4-sided pyramid":pyramid,"triangular prism": triprism,"triangular double pyramid":doubletetra}
denpolyhedra = {"triangular double pyramid":doubletetra}
enpolyhedra = {"4-sided pyramid":pyramid}
polyhedra = {"tetrahedron":tetra,"cube":cube,"4-sided pyramid":pyramid,"triangular prism": triprism,"triangular double pyramid":doubletetra}
octa = {"octahedron":octahedron}

def combine(lf):
    """Combine a list of faces into a single string."""
    s = ""
    for x in lf:
        s += x
    return s

def others(tt):
    ntt = ""
    for e in tt:
        ne = other(e)
        ntt += ne
    return ntt

def other(e):
    """Reverse the edge"""
    if e == e.lower():
        return e.upper()
    else:
        return e.lower()

def noparallel(sequ):
    """Number of parallel vs antiparallel edges"""
    n = 0
    for x in sequ:
        if other(x) in sequ:
            n += 1
    return (n//2, (len(sequ)-n)//2)
    
#def np(lf):
#    """Same as noparallel except that for multiple faces"""
#    return noparallel(combine(lf))

def opposite(face):
    """Reverse the face orientation.""" 
    nf = ""
    for x in face:
        nf += other(x)
    return nf[::-1]

def findsingles(face):
    """Find all letters in the face with a sigle occurence."""
    res =""
    for e in face:
        if face.count(e) == 1 and face.count(other(e)) == 0:
            res += e
    return res

def movesingle(face,e):
    """Move edge e to the front of face face."""
#    print(face,e)
    i = face.index(e)
    return face[i:]+face[:i]

def normalize(face):
    """Find first letter with single occurrence and move it to the front of face face."""
    sings = findsingles(face)
    if sings == "":
        print("Error:",face," - no single occurrence")
        raise
    e = choice(sings)
#    print(sings,e)
    return movesingle(face,e)

def glueorsplit(lf,ss):
    """glue or split lf along the symbols in the list ss"""
    nf = lf
    for s in ss:
#        print(nf,s)
        nf = glorsp(nf,s)
#        print(nf,s)
    return nf

def glorsp(nf,s):
    """glue or split along the symbols s"""
    for i,x in enumerate(nf):
        if s in x or other(s) in x:
            for j,y in enumerate(nf[i+1:]):
                if s in y or other(s) in y:
#                    print("glue")
                    return nf[:i]+nf[i+1:j+i+1]+nf[i+j+2:]+[gluealong(x,y,s)]
#            print("split")
            return nf[:i]+nf[i+1:]+split(x,s)
                    
            

def split(f,e):
    """ split f along an edge e."""
    E = other(e)
    if f.count(e) == 0:
        e = E
    fn = movesingle(f,e)
    if fn.count(e) == 2:
        i = fn[1:].index(e)+1
        f1 = fn[:i]
        f2 = fn[i:]
        return[f1,f2]
    E = other(e)
    i = fn.index(E)
    f1 = fn[:i]
    f2 = E+opposite(fn[i+1:])
    return [f1+f2]

def gluealong(f1,f2,e):
    """Glue faces f1 and f2 along the edge e."""
#    print(e)
    if f1.count(e) == 0:
        f1 = opposite(f1)
    f1 = movesingle(f1,e)
    if f2.count(e) == 0:
        f2 = opposite(f2)
    f2 = movesingle(f2,e)
    return f1+f2

def gluetwo(lf):
    """Glue two faces from the list of faces. One of them is the first face in lf."""
    nlf = [normalize(f) for f in lf]
    shuffle(nlf)
    f1 = nlf[0]
    e = f1[0]
    E = other(e)
    for i,f in enumerate(nlf[1:]):
        if e in f or E in f:
            j = i+1
            f2 = f
            break
    gf = gluealong(f1,f2,e)
#    print(j)
    if j == 1:
        return [gf] + nlf[2:]
    return [gf]+nlf[1:j]+nlf[j+1:]

## FINDS OPTIMIMAL STARTING POSITION
def distr(tt):
    """Computes the maximum distance, the sum of distances and the distribution."""
    res = {}
    ttt = [x.upper() for x in tt]
#   print(ttt)
    st = set(ttt)
    en = 0
    km = -1
    xm = ""
    suma = 0
    wsuma = 0
    for x in st:
        i = ttt.index(x)
        j = ttt.index(x,i+1)
        k = j-i
        suma += k
        wsuma += k*i
        res[k] = res.get(k,0)+1
        if k > en:
            en = k
            xm = x
            km = i
    return (en,suma,res,wsuma)



def maxdist(tt):
    """Computes the maximum distance of two letters in a string."""
    ttt = [x.upper() for x in tt]
#   print(ttt)
    st = set(ttt)
    en = 0
    km = -1
    xm = ""
    for x in st:
        i = ttt.index(x)
        j = ttt.index(x,i+1)
        k = j-i
        if k > en:
            en = k
            xm = x
            km = i
    return (en,xm,km)

def minmaxdist(tt):
    """ Returns a triple:
minmaxdist
shifted tt in an optimal way.
index at whivch the shift is needed."""
    men = len(tt)
    mi = -1
    mtt = ""
    for i,x in enumerate(tt):
        ttt = tt[i:]+ tt[:i]
        (en,xm,km) = maxdist(ttt)
        if en < men:
            men = en
            mtt = ttt
            mi = i
    return(men,mtt,mi)

def minsumdist(tt):
    """ Returns a triple:
minsumdist
shifted tt in an optimal way.
index at whivch the shift is needed
minsum."""
    men = len(tt)*len(tt) +10000 # TODO: Nino was here! quick and dirty fix
    mi = -1
    mtt = ""
#    print(men)
    for i,x in enumerate(tt):
        ttt = tt[i:]+ tt[:i]
        (en,suma,res,wsuma) = distr(ttt)
#        print(wsuma)
#        print(en,xm,km)
#        print(ttt,i,en)
        if wsuma < men:
            men = wsuma
            mtt = ttt
            mi = en
    return(men,mtt,mi)


def mmoptimal(tt):
    """ Return the shifted string with minimal maximal distance."""
    (men,mtt,mi) = minmaxdist(tt)
    return mtt

# The surface associated with the fundamental polygon.
# Applications to TET12
# T. Pisanski, 17.5. 2013

# Examples
tet12 = "acbGAhGcehBe" # Fundamental polygon for TET12
tet12a = "acbadebfeCdF" # Fundamental polygon for TET12_1
tet12b = "acbadeCdfBEf" # Fundamental polygon for TET12_2
tet12c = "adfacefbCdEB" # Fundamental polygon for Tet12_3
triangle = "abcabc"
torus = "abAB"
klein = "abAb"
pplane = "aa"
torus1 = "abcABC"
doubletorus = "abABcdCD"
tripletorus = "abABcdCDefEF"
doubleproj = "aabb"
tripleproj = "aabbcc"

polygons = [tet12,tet12a,tet12b,tet12c,triangle,torus,klein,pplane,torus1,doubletorus,tripletorus,doubleproj,tripleproj]
# polygons = [tet12,tet12a,tet12b,tet12c]


# Functions
def edges(tt):
    """Edges of the fundamental polygon"""
    n = len(tt)
    e=[]
    for i in range(n-1):
        e += [(i,i+1)]
    e += [(n-1,0)]
    return e

def npermut(tt):
    """ Number of permutations for a given sequence and self-paired segments."""
    (anti,par) = noparallel(tt)
    n2 = len(tt)
    n = n2//2
    k = n # number of different types of segments.
    print()
    print("nseg  nself nnons nadmissible nall")
    while k <= n2:
        s = n2 - k # number of self-paired segment types.
        r = n - s # number of non self-paired segment types.
        num = factorial(anti)*factorial(par)*(2**r)
        numall = factorial(n2)//(2**par)
        fract = num/numall
        print(k,s,r,num,numall,fract)
        k += 1
    print()     
        
def factorial(n):
    if n <= 0:
        return 1
    return n*factorial(n-1)


def surf(tt):
    """cn - map to canonical values,
vts - vertices,
edgs - doubled edges
v - number of vertices
e - number of edges
f - number of faces (=1)
chi - Euler characterisitc.
oriented- is the surface oriented (= orientable)?
genus - genus of the surface.
par - number of parallel edges
anti - number of antiparallel edges
"""
    (cn,vts,edgs,gp,prt) = skeleton(tt)
    v = len(vts)
    e = len(edgs)//2
    ts = set(tt)
    oriented = len(tt) == len(ts)
    chi = v - e + 1
    genus = 2 - chi
    if oriented:
        genus = genus//2
    (anti,par) = noparallel(tt)
    return(anti,par,v,e,1,chi,oriented,genus)


def doubletrace(tt):
    """Determine the double trace corresponding to tt"""
    (cn,vts,edgs,gp,prt) = skeleton(tt)
    v = len(vts)
    e = len(edgs)//2
#    print(v,e)
    ts = set(tt)
    oriented = len(tt) == len(ts)
    chi = v - e + 1
    genus = 2 - chi
    if oriented:
        genus = genus//2
    (anti,par) = noparallel(tt)
    trail = [vts.index(edgs[0][0])]
    for (u,w) in edgs:
        trail += [vts.index(w)]
    return(vts,trail)


def surface(tt):
    """cn - map to canonical values,
vts - vertices,
edgs - doubled edges
v - number of vertices
e - number of edges
f - number of faces (=1)
chi - Euler characterisitc.
oriented- is the surface oriented (= orientable)?
genus - genus of the surface.
par - number of parallel edges
anti - number of antiparallel edges
"""
    (cn,vts,edgs,gp,prt) = skeleton(tt)
    dt = doubletrace(tt)
    v = len(vts)
    e = len(edgs)//2
    ts = set(tt)
    oriented = len(tt) == len(ts)
    chi = v - e + 1
    genus = 2 - chi
    if oriented:
        genus = genus//2
    (anti,par) = noparallel(tt)
    trail = [edgs[0][0]]
    for (u,w) in edgs:
        trail += [w]
    print("Eulerian trail: ", trail)
    print()    
    print("fundamental polygon = ",tt)
    print("Double trace = ",dt)
    print("glued vertices = ",gp)
    print("canonical labels = ", cn)
    print("partition = ", prt)
    print("vertices = ",vts)
    print("Eulerian walk = ",edgs)
    print("Number of parallel edges = ",par)
    print("Number of antiparallel edges = ",anti)
#    npermut(tt)
    print( "v = ",v)
    print( "e = ",e)
    print( "f = ",1)
    print( "chi = ",chi)
    print( "oriented = ",oriented)
    print("genus = ", genus)
    print()
    return(anti,par,v,e,dt)

def surfaces(poly):
    start = clock()
    tts = generateall(poly)
    nsol = len(tts)
    stop = clock()
    print("Elapsed time: ","{0:10.4f}".format(stop-start))
    print("Input faces: ",poly)
    print("Number of solutions: ",nsol)
    res = []
    for x in tts:
        # print(x)
        (su,y,t) = minsumdist(x)
        # print(su, y, t)
        res.append((su,doubletrace(y)[1],surf(x)))
    res.sort()
    for r in res:
        if (r[2][1] != 0): continue
        print(r)
    final = clock()
    print("Total time: ","{0:10.4f}".format(final-start))
    
    return res


# Tests    
def tests(p=polygons):
    print("Skeleton and surface from fundamental polygon.") 
    print("==============================================") 
    for pp in p:
        surface(pp)

def experiment(tt,n=1):
    """given a set of faces of a convex polyhedron, find its polypeptide represetnations."""
    res = []
    some = makesingle(tt)
    nn = len(some)
    k = min(5,nn)
    for i in range(n):
        samp = sample(some,k)
#        print("sample = ",samp)
        sx = makesingle(tt,samp)
        found = False
        for x in res:
            if isomorphic(x,sx):
                found = True
#                print(isomorphism(x,sx))
                break
        if not found:
            res.append(mmoptimal(sx))
#            print(i)
    return res

def printall(poly):
    """Print matching indices of polypeptides and their reverses for a given polyhedron."""
    dts = generateall(poly)
    for i,x in enumerate(dts):
            xx = x[::-1]
            for j,y in enumerate(dts):
                    if isomorphic(xx,y):
                            print(i,j)

def newstrings(tests,olds):
    """ From tests select those that are not isomorphic to ones in olds."""
    res = []
    for x in tests:
        for y in olds:
            if isomorphic(x,y):
                break
        else:
            res += [x]
    return res
        
def test(ll=[doubletetra]):
    sp = []
    for x in ll:
        sx = makesingle(x)
        print(sx,noparallel(sx))
        sp.append(sx)
        print(minmaxdist(sx))
        print(surface(sx))
    return sp

# test(polyhedra)
# test([triprism])

# surfaces(triprism)

# surfaces(tetra)

# bipiramida
bipi = ["aeD", "cdF", "bfE", "ahG", "biH", "cgI"]
# surfaces(bipi)    

# 4-pyramid
four_py = ['afE', 'bgF', 'chG', 'deH', 'abcd']
# surfaces(four_py) 

four_py = ['AFe', 'BGf', 'CHg', 'DEh', 'ABCD']
# surfaces(four_py) 


# print(octahedron)
# surfaces(octahedron)
# surfaces(['abc', 'abc'])
# minsumdist('fGClIdhEAfEdlKChGbjIAjKb')

# cube = ["abcd","ijkl","aeIH","bfJE","cgKF","dhLG"]
# surfaces(cube)

bipi_list = [
    ['abcdef', 'gDFhBiEGACHI'],
    ['abcdefBgChDiEHF', 'GAI'],
    ['abcdef', 'gDFh', 'CiEGA', 'IBH'],
    ['abcdefgBDhFiG', 'HCIEA'],
    ['abcdefBgEhACFi', 'HDIG'],
    ['abcdAefg', 'hEiCFHDIGB'],
]

crazy_tetra = ['abcdAeCf', 'DEFB']

# for bp in bipi_list:
 #   surfaces(bp)

bipi3 = [
    ['abcdefghBDFiH', 'GCIEA'],
    ['abcdef', 'ghFB', 'HCEiA', 'GID'],
    ['abcdef', 'gCEhAiGHDIFB'],
    ['abcdAeCfgDhFi', 'GHEIB'],
    ['abcdef', 'gDhAiCEGHIFB'],
    ['abcdef', 'gDhFB', 'HiA', 'GICE'],
    ['abcde', 'fCgEh', 'GiDFA', 'IBH'],
    ['abcdefgBEhGiD', 'CHIFA'],
    ['abcdeAfBgChDi', 'EHFIG'],
    ['abcdAefg', 'hiCFHDIEGB'],
    ['abcdAefg', 'hEGB', 'HDi', 'CFI'],
    ['abcdAefDgEhBFiH', 'CIG'],
    ['abcdefCgEAhGi', 'DFHIB'],
    ['abcdefgBDhEAFiG', 'CIH'],
    ['abcdefg', 'DhiGB', 'IEA', 'CFH'],
    ['abcdeBfDgAhEi', 'GCHIF'],
    ['abcdefgBEhFACiD', 'HIG'],
    ['abcd', 'efgDh', 'GBHiEA', 'ICF'],
    ['abcdAefBghCi', 'EIFGDH'],
    ['abcdAefBgEhFiCH', 'GDI'],
    ['abcdAefDgEhiGCH', 'FIB'],
    ['abcdefCg', 'DFhGiEAHIB'],
    ['abcdefghiCG', 'DIEAFHB'],
    ['abcdefBDghiFGCI', 'HEA'],
    ['abcde', 'fCgEhBiDFAGHI'],
    ['abcdefgD', 'hFBiGHACEI'],
    ['abcdefghCFiH', 'GBDIEA'],
]

pyramid_4 = [['abcdAefg', 'BDEGh', 'CHF'], ['abcdefBgCFhGADHE'], ['abcde', 'BfgEh', 'AGCHDF']]
"""
for x in pyramid_4:
    print("### SEP ### True")
    surfaces(x)
"""

##for x in bipi3:
##    surfaces(x)

# print(gluealong('abc', 'afg', 'a'))

def tast(poly=enpolyhedra):
    print("tast")
    print("poly : ",poly)
    res = []
    for x in poly:
        p = poly[x]
        print("x : ",x)
        print("p : ",p)
        resp = []
        of = makesingle(p)
        print("of : ",of)
        ae = all_embeddings(of)
        print(x,len(ae))
        alls = []
        for rot in ae:
            print(rot)
            faces = facesfromrotation(rot)
            sfac = stringsfromfaces(faces)
            single = makesingle(sfac)
            print("+",sfac,single)
            resp.append(single)
            ga = [canonical(x) for x in generateall([single])]
            freq = {}
            print(len(ga))
            for s in ga:
                if s in alls:
                    print("seen")
                    break
                alls.append(s)
                (anti,par) = noparallel(s)
                freq[anti] = freq.get(anti,0)+1
                if False:
                    print(anti,s)
            print(len(ga))
            ks = list(freq.keys())
            ks.sort()
            for x in ks:
                print(x,freq[x])
        res.append(resp)
    return (poly,res)


def tist(poly=smallpolyhedra):
    """Stil not finished version!!!"""
    res = []
    for x in poly:
        p = poly[x]
        resp = []
        of = makesingle(p)
        ae = all_embeddings(of)
        print(x,len(ae))
        alls = []
        for rot in ae:
            print(rot)
            faces = facesfromrotation(rot)
            sfac = stringsfromfaces(faces)
            single = makesingle(sfac)
            print("+",sfac,single)
            resp.append(single)
            ga = [canonical(x) for x in generateall([single])]
            freq = {}
            print(len(ga))
            for s in ga:
                if s in alls:
                    print("seen")
                    break
                alls.append(s)
                (anti,par) = noparallel(s)
                print(anti,s)
                freq[anti] = freq.get(anti,0)+1
                if False:
                    print(anti,s)
            print(len(ga))
            ks = list(freq.keys())
            ks.sort()
            for x in ks:
                print(x,freq[x])
        res.append(resp)
    return (poly,res)
        
#(p,t) = tast()
#print("Try to run explore(tetra).")
