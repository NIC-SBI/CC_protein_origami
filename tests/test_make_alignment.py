import cocopod
import cocopod.utils as u 

#TODO: fix test
def test_make_alignemnt():
    import subprocess
    import filecmp
    import Bio.SeqIO
    cmd_line = "python ../cocopod/make_alignment.py --json data/APHsh.json -aln data/APHtest.ali > out.log"
    subprocess.call(cmd_line, shell=True)
    alnf = u.relative_to(__file__, 'data/APHtest.ali')
    with open(alnf, 'r') as f1:
        li = []
        for i in Bio.SeqIO.PirIO.PirIterator(f1):
            li.append(i.seq)
    assert li[0] == "ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK/-----ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK"
    assert li[1] == "ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLKGSGSGSELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK"
      
  
