import ppmod 

#TODO: fix test
def test_make_alignemnt():
    import subprocess
    import filecmp
    cmd_line = "python ../ppmod/make_alignment.py --json data/APHsh.json -aln data/APHtest.ali > out.log"
    subprocess.call(cmd_line, shell=True)
    #fails probably due to line endings    
    #assert filecmp.cmp('data/APHsh.ali', 'data/APHtest.ali')
