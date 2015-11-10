import ppmod 

def test_make_alignemnt():
    import subprocess
    import filecmp
    import Bio.SeqIO
    cmd_line = "python ../ppmod/make_alignment.py --json data/APHsh.json -aln data/APHtest.ali > out.log"
    subprocess.call(cmd_line, shell=True)
    with open('data/APHtest.ali', 'r') as f1:
        li = []
        for i in Bio.SeqIO.PirIO.PirIterator(f1):
            li.append(i.seq)
    assert li[0] == "ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK/-----ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK"
    assert li[1] == "ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLKGSGSGSELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK"
