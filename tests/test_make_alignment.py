import ppmod 

def test_make_alignemnt():
    import subprocess
    import filecmp
    cmd_line = "python ../ppmod/make_alignment.py --json data/APHsh.json -aln data/APHtest.ali > out.log"
    subprocess.call(cmd_line, shell=True)
    assert filecmp.cmp('data/APHsh.ali', 'data/APHtest.ali')
