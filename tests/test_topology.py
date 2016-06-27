import ppmod
import ppmod.utils as u
from ppmod.topology import *

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