import ppmod
import ppmod.utils as u

def test_relative_to():
  """Test file name handling"""
  import os
  aln_file = u.relative_to(__file__, 'data/test-aln.ali')
  #this file must exist in this relative path (since we put it there)
  assert os.path.isfile(aln_file) 

def test_sequnce_and_knowns():
  aln_file = u.relative_to(__file__, 'data/test-aln.ali')
  seq, knowns = u.sequnce_and_knowns(aln_file)
  assert seq=='tet-p5LD3' 
  #use set so order is not important
  assert set(knowns)==set(['p3_p4', 'p5_p6', 'p7_p8', 'APHshort', 'BCRshort', 'GCN'])

def test_sequnce_and_knowns_spaces():
  aln_file = u.relative_to(__file__, 'data/test-aln-spaces.ali')
  seq, knowns = u.sequnce_and_knowns(aln_file)
  assert seq=='tet-p5LD3' 
  #use set so order is not important
  assert set(knowns)==set(['p3_p4', 'p5_p6', 'p7_p8', 'APHshort', 'BCRshort', 'GCN'])  