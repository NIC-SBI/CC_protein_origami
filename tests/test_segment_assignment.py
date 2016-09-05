import cocopod
import cocopod.utils as u
from cocopod.segment_assignment import *

from collections import OrderedDict

def test_get_complete_pairs_dict_from_topology():
    from collections import OrderedDict
    assert get_complete_pairs_dict_from_topology("ABCDabcd") == \
                    OrderedDict([('A', ['A', 'a']),
                                 ('B', ['B', 'b']),
                                 ('C', ['C', 'c']),
                                 ('D', ['D', 'd'])])

def test_get_pairs_from_topology():
    assert(get_pairs_from_topology("ABCDabcd") == ['A', 'B', 'C', 'D'])
    assert(get_pairs_from_topology("A-B-C-D-a-b-c-d".split("-")) == ['A', 'B', 'C', 'D'])
    
    
def test_segment_assignments_to_dict():
    

    seg_str="""
    A->SEG1:SEG2
    B->SEGA:SEGA    
    """

    d = segment_assignments_to_dict(seg_str)
    assert d == OrderedDict([(u'A', [u'SEG1', u'SEG2']), (u'B', [u'SEGA', u'SEGA'])])

    #test spaces and empty lines
    seg_str="""
    A -> SEG1:SEG2
    
    B->SEGA:SEGA    
    
    C-> SEGB : SEGB    
    """
    d = segment_assignments_to_dict(seg_str)

    d = segment_assignments_to_dict(seg_str)
    assert d == OrderedDict([(u'A', [u'SEG1', u'SEG2']),
                             (u'B', [u'SEGA', u'SEGA']),
                             (u'C', [u'SEGB', u'SEGB'])])


def test_do_assignment_replacements():
    rep_rules = """
    A->1:2
    B->three:four
    C->FIVE : SIX
    D->7:8"""
    assert do_assignment_replacements("ABCDABcd", rep_rules) == \
      ['1', 'three', 'FIVE', '7', '2', 'four', 'SIX', '8']
      
    assert do_assignment_replacements("ABCDABcd", segment_assignments_to_dict(rep_rules)) == \
      ['1', 'three', 'FIVE', '7', '2', 'four', 'SIX', '8']  
      
def test_seq_to_seq_map():
    assert seq_to_seq_map({'1':'ABC', '2':'CDEF'}) == {'1':'ABC', '2':'CDEF'}
    
    xls_file = u.relative_to(__file__, 'data/segments.xlsx')
    assert seq_to_seq_map(xls_file) == \
          {u'P1SN': u'SPED EIRQLEQ ENSQLER ENQRLEQ EIYQLER',
           u'P2SN': u'SPED KIEELKE KNSQLKE KNEELKQ KIYELKE',
           u'P3SN': u'SPED EIQQLEE EISQLEQ KNSELKE KNQELKY',
           u'P4SN': u'SPED KISQLKE KIQQLKQ ENQQLEE ENSQLEY'}

def test_get_annotated_sequence():
    xls_file = u.relative_to(__file__, 'data/segments.xlsx')    
    #TODO increase readibility    
    result="""M
SPED EIRQLEQ ENSQLER ENQRLEQ EIYQLER\t|P1SN
SGPGS
SPED EIQQLEE EISQLEQ KNSELKE KNQELKY\t|P3SN
SGPGS
SPED KISQLKE KIQQLKQ ENQQLEE ENSQLEY\t|P4SN
SGPGS
SPED KIEELKE KNSQLKE KNEELKQ KIYELKE\t|P2SN
LEHHHHHHHH"""

    assert get_annotated_sequence("P1SN-P3SN-P4SN-P2SN".split("-"), 
              xls_file, N_tag="M",C_tag="LEHHHHHHHH") == result
              
               
def test_deannotate_sequence():
    annotated_sequence="""M
SPED EIRQLEQ ENSQLER ENQRLEQ EIYQLER\t|P1SN
SGPGS
SPED EIQQLEE EISQLEQ KNSELKE KNQELKY\t  |P3SN
SGPGS
SPED KISQLKE KIQQLKQ ENQQLEE ENSQLEY\t\t\t\t|  P4SN
SGPGS
SPED KIEELKE KNSQLKE KNEELKQ KIYELKE  |  P2SN
LEHHHHHHHH"""

    result="""M
SPED EIRQLEQ ENSQLER ENQRLEQ EIYQLER
SGPGS
SPED EIQQLEE EISQLEQ KNSELKE KNQELKY
SGPGS
SPED KISQLKE KIQQLKQ ENQQLEE ENSQLEY
SGPGS
SPED KIEELKE KNSQLKE KNEELKQ KIYELKE
LEHHHHHHHH"""

    assert deannotate_sequence(annotated_sequence)==result
    
    assert deannotate_sequence(annotated_sequence, remove_whitespace=True) ==\
    "MSPEDEIRQLEQENSQLERENQRLEQEIYQLERSGPGSSPEDEIQQLEEEISQLEQKNSELKEKNQELKYSGPGSSPEDKISQLKEKIQQLKQENQQLEEENSQLEYSGPGSSPEDKIEELKEKNSQLKEKNEELKQKIYELKELEHHHHHHHH"
    
    


def test_segment_assignment_gui():
    """Test the creation of the GUI"""
    import traitlets    
    try:
        gui = segment_assignment_gui("bDfABCADECFE", auto_display=False)
    except traitlets.traitlets.TraitError:
        pass
    #assert gui
    
def test_sequence_edit_gui():
    """Test the creation of the GUI"""
    import traitlets    
    try:
        gui = sequence_edit_gui("A Text", "A caption", auto_display=False)   
    except traitlets.traitlets.TraitError:
        pass
        