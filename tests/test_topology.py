import ppmod
import ppmod.utils as u
from ppmod.topology import *
import numpy

def test_permute_segment_left():
    assert permute_segment_left(['a', 'b', 'c']) == \
                            ['b', 'c', 'a']
    

def test_permute_vertices_left():
    assert permute_vertices_left([1,2,3,1], 1) == [2, 3, 1, 2]
    assert permute_vertices_left([1,2,3,1], 2) == [3, 1, 2, 3]  
    
def test_get_permutation_name():
    assert get_permutation_name('1',1) == "1.1"
    assert get_permutation_name('1',2) == "1.2"
    assert get_permutation_name(1,1) == "1.1"

def test_get_segment_distances():
    assert get_segment_distances("a-b-c-a-d-e-f-B-d-g-h-f-c-i-h-E-g-I".split('-')) == \
    [3, 5, 10, 4, 9, 5, 7, 4, 3] 
    assert get_segment_distances(list("cBAaBC")) == \
    [0, 3, 4] 

def test_explore_tetrahedron():
    tops = explore(tetrahedron,verbose=False, return_dataframe=False)
    assert tops == [{'num_AP': 2,
                      'num_P': 4,
                      'num_cross': 0,
                      'reflected': True,
                      'topo': 'ABCADECFEbDf'},
                     {'num_AP': 3,
                      'num_P': 3,
                      'num_cross': 0,
                      'reflected': True,
                      'topo': 'ABCADECFdBef'},
                     {'num_AP': 3,
                      'num_P': 3,
                      'num_cross': 0,
                      'reflected': True,
                      'topo': 'ABCADEbDFceF'}]

def test_name_topologies_and_permutations():
    tops = explore(tetrahedron)
    named_tops = name_topologies_and_permutations(tops) 
    assert len(named_tops)==36
    assert numpy.all(named_tops.index.values == [u'1.1', u'1.2', u'1.3', u'1.4', u'1.5', u'1.6', u'1.7', u'1.8', u'1.9',
       u'1.10', u'1.11', u'1.12', u'2.1', u'2.2', u'2.3', u'2.4', u'2.5',
       u'2.6', u'2.7', u'2.8', u'2.9', u'2.10', u'2.11', u'2.12', u'3.1',
       u'3.2', u'3.3', u'3.4', u'3.5', u'3.6', u'3.7', u'3.8', u'3.9', u'3.10',
       u'3.11', u'3.12'])
    assert numpy.all(named_tops.segments.values == ['ABCADECFEbDf', 'BCADECFEbDfA', 'CADECFEbDfAB', 'ADECFEbDfABC',
       'DECFEbDfABCA', 'ECFEbDfABCAD', 'CFEbDfABCADE', 'FEbDfABCADEC',
       'EbDfABCADECF', 'bDfABCADECFE', 'DfABCADECFEb', 'fABCADECFEbD',
       'ABCADECFdBef', 'BCADECFdBefA', 'CADECFdBefAB', 'ADECFdBefABC',
       'DECFdBefABCA', 'ECFdBefABCAD', 'CFdBefABCADE', 'FdBefABCADEC',
       'dBefABCADECF', 'BefABCADECFd', 'efABCADECFdB', 'fABCADECFdBe',
       'ABCADEbDFceF', 'BCADEbDFceFA', 'CADEbDFceFAB', 'ADEbDFceFABC',
       'DEbDFceFABCA', 'EbDFceFABCAD', 'bDFceFABCADE', 'DFceFABCADEb',
       'FceFABCADEbD', 'ceFABCADEbDF', 'eFABCADEbDFc', 'FABCADEbDFce'])

def test_TCO():
    tops = explore(tetrahedron)
    named_tops = name_topologies_and_permutations(tops) 
    named_tops = calculate_TCO(named_tops)
    assert len(named_tops)==36
    assert numpy.all(named_tops.index.values == [u'1.1', u'1.2', u'1.3', u'1.4', u'1.5', u'1.6', u'1.7', u'1.8', u'1.9',
       u'1.10', u'1.11', u'1.12', u'2.1', u'2.2', u'2.3', u'2.4', u'2.5',
       u'2.6', u'2.7', u'2.8', u'2.9', u'2.10', u'2.11', u'2.12', u'3.1',
       u'3.2', u'3.3', u'3.4', u'3.5', u'3.6', u'3.7', u'3.8', u'3.9', u'3.10',
       u'3.11', u'3.12'])
    assert numpy.all(named_tops.segments.values == ['ABCADECFEbDf', 'BCADECFEbDfA', 'CADECFEbDfAB', 'ADECFEbDfABC',
       'DECFEbDfABCA', 'ECFEbDfABCAD', 'CFEbDfABCADE', 'FEbDfABCADEC',
       'EbDfABCADECF', 'bDfABCADECFE', 'DfABCADECFEb', 'fABCADECFEbD',
       'ABCADECFdBef', 'BCADECFdBefA', 'CADECFdBefAB', 'ADECFdBefABC',
       'DECFdBefABCA', 'ECFdBefABCAD', 'CFdBefABCADE', 'FdBefABCADEC',
       'dBefABCADECF', 'BefABCADECFd', 'efABCADECFdB', 'fABCADECFdBe',
       'ABCADEbDFceF', 'BCADEbDFceFA', 'CADEbDFceFAB', 'ADEbDFceFABC',
       'DEbDFceFABCA', 'EbDFceFABCAD', 'bDFceFABCADE', 'DFceFABCADEb',
       'FceFABCADEbD', 'ceFABCADEbDF', 'eFABCADEbDFc', 'FABCADEbDFce'])    
    assert numpy.all(named_tops.TCO.values == [4.333333333333333, 5.333333333333333, 4.666666666666667,
       5.333333333333333, 4.333333333333333, 4.333333333333333,
       5.333333333333333, 4.666666666666667, 5.333333333333333,
       4.333333333333333, 5.0, 5.0, 4.166666666666667, 5.166666666666667,
       4.5, 5.166666666666667, 4.166666666666667, 4.833333333333333,
       5.166666666666667, 4.5, 5.166666666666667, 4.5, 5.166666666666667,
       4.833333333333333, 3.8333333333333335, 4.833333333333333,
       5.166666666666667, 4.833333333333333, 3.8333333333333335,
       4.833333333333333, 5.166666666666667, 4.833333333333333,
       3.8333333333333335, 4.833333333333333, 5.166666666666667,
       4.833333333333333])

def test_name_of_topology():    
    tops = explore(tetrahedron)
    named_tops = name_topologies_and_permutations(tops) 
    df = calculate_TCO(named_tops)

    pTet12 = "ABCDaEDBFEcF"
    assert name_of_topology(pTet12, df)=='1.10'

    p5 = "A-B-C-A-D-E-c-F-d-B-F-E".replace("-","")
    name_of_topology(p5, df)
    assert name_of_topology(p5, df)=='1.6'

    p33a = "A-B-C-D-A-E-c-F-d-e-B-F".replace("-","")
    assert name_of_topology(p33a, df)=='2.3'
    
def test_convert_vface_to_efaces():
    #tetrahedron    
    vfaces = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
    assert ['Abc', 'CDe', 'Efa', 'BFd']==convert_vface_to_efaces(vfaces)
    
    #pyramid
    vfaces =[[1, 2, 3, 0], [0, 1, 4], [0, 4, 3], [3, 4, 2], [1, 2, 4]]
    assert ['ABcD', 'DEf', 'Fgc', 'GhB', 'AHe']==convert_vface_to_efaces(vfaces)
    
    #prism
    vfaces = [[0, 1, 2], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [0, 2, 5, 3]]
    assert ['ABc', 'DEf', 'AGdh', 'BIeg', 'CIfh']==convert_vface_to_efaces(vfaces)