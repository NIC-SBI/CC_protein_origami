import cocopod
import cocopod.utils as u

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
  i,ii = u.align('IQQLEEEIAQLEQKNAALKEKNQALKYG', 'SPEDEIKELEEEIKELEWKNEELKRKNEELKRG', 6, 6)
  assert (i, ii) == (0, 5)
  i,ii = u.align('ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK', 'REKELQKIEEQKAQLQWKAQARKEKLAQLK', 6, 6)
  assert (i, ii) == (5, 1)
  i, ii = u.align('ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK', 'ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLKEKL', 6, 6)
  assert (i, ii) == (0, 0)


def test_find_pair():
  json = u.relative_to(__file__, 'data/data.json')
  d = u.load_json_data(json)
  p1, p2, pdb = u.find_pair('APHshE', d.segments)
  assert p1 == d.segments[4]['id']-1

def test_score():
  score = u.score('KKLLLVQEI', 'LEELLEQEK')
  assert score == 4
  
def test_is_str():
    assert u.is_str("A string")
    assert u.is_str(u"A string")
    assert not u.is_str({})
    assert not u.is_str([])
    assert not u.is_str(1)
    

#TODO: fix test
def writepdb():
   import mdtraj as md
   import os

   pdbname = u.relative_to(__file__,'data/APH.pdb')
   pdbname1 = u.relative_to(__file__,'data/APH1.pdb')


   md_obj = md.load(pdbname)
   topology = md_obj.topology
   positions = md_obj.xyz
   seq_target = u.mdtraj_to_fasta(topology,0)[2:17]
   u.writepdb(2, 14, topology, positions, pdbname, 1, 2)
   assert os.path.isfile(pdbname1)
   md_obj = md.load(pdbname1)
   topology = md_obj.topology
   seq_test = u.mdtraj_to_fasta(topology,0)
   assert seq_test == seq_target

def test_mdtraj_to_fasta():
    import mdtraj as md
    topology = md.load(u.relative_to(__file__, 'data/p3_p4.pdb')).topology   
    chain_A = u.mdtraj_to_fasta(topology,0)
    chain_B = u.mdtraj_to_fasta(topology,1)
    assert chain_A=="IQQLEEEIAQLEQKNAALKEKNQALKYG"
    assert chain_B=="IAQLKQKIQALKQENQQLEEENAALEYG"

def test_mdtraj_to_fasta_bcr():
    """Test for reading MSE residue"""
    import mdtraj as md    
    topology = md.load(u.relative_to(__file__, 'data/BCR.pdb')).topology   
    
    chain_A = u.mdtraj_to_fasta(topology,0)
    chain_B = u.mdtraj_to_fasta(topology,1)
    
    assert chain_A=="DIEQELERAKASIRRLEQEVNQERFRMIYLQTLLAK"
    assert chain_B=="DIEQELERAKASIRRLEQEVNQERFRMIYLQTLLAK"


def test_roundrobin():
    assert list(u.roundrobin(list('ABC'), list('D'), list('EF'))) == list("ADEBFC")
    assert list(u.roundrobin(list('ABCD'), list('---'))) == list("A-B-C-D")