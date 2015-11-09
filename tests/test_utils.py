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

def test_align():
  i,ii = u.align('IQQLEEEIAQLEQKNAALKEKNQALKYG','SPEDEIKELEEEIKELEWKNEELKRKNEELKRG',6,6)
  assert (i, ii) == (0, 5) 
  i,ii = u.align('ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK','REKELQKIEEQKAQLQWKAQARKEKLAQLK',6,6)
  assert (i, ii) == (5, 1) 
  i, ii = u.align('ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK','ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLKEKL',6,6)
  assert (i, ii) == (0,0)     
                          
def test_find_pair():
  json = u.relative_to(__file__, 'data/data.json')
  d = u.load_json_data(json)
  p1, p2, pdb = u.find_pair('APHshE', d.segments)
  assert p1 == d.segments[4]['id']-1

def test_score():
  score = u.score('KKLLLVQEI','LEELLEQEK')
  assert score == 4


def test_mdtraj_to_fasta():
    import mdtraj as md
    topology = md.load(u.relative_to(__file__, 'data/p3_p4.pdb')).topology   #read topology
    chain_A = u.mdtraj_to_fasta(topology,0)
    chain_B = u.mdtraj_to_fasta(topology,1)
    assert chain_A=="IQQLEEEIAQLEQKNAALKEKNQALKYG"
    assert chain_B=="IAQLKQKIQALKQENQQLEEENAALEYG"
