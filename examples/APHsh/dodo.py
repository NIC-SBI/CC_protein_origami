N_fold = 1
N_homology = 3 


from make_config import model_name
import os
import multiprocessing
import cocopod


MODEL = model_name
SDIR = os.path.dirname(cocopod.__file__)
print(SDIR)
print(os.getcwd())

DOIT_CONFIG = {'default_tasks': ['make_homology_models'], 'num_process': multiprocessing.cpu_count(), 'par_type':'thread',}

from doit.tools import create_folder
from doit import get_var

N_fold = int(get_var('N_fold', N_fold))
N_homology = int(get_var('N_homology', N_homology))


def task_make_config():    
    return { 'file_dep': ['make_config.py'],
             'targets':  [MODEL+".json"],
             'actions':  ["python %(dependencies)s"],
             'clean' : True}

def task_make_helix():
    return {'file_dep': [MODEL+".json"],
            'targets': [MODEL+"-helix.pdb"],
            'actions':  ["python "+SDIR+"/make_helix.py --json %(dependencies)s --out-file %(targets)s"],
            'clean' : True}
            
def task_make_alignment():
    return {'file_dep': [MODEL+".json"],
            'targets': [MODEL+".ali"],
            'actions':  ["python "+SDIR+"/make_alignment.py --json %(dependencies)s -aln %(targets)s"],
            'clean' : True}            

def task_make_initial_models():
    for nf in range(N_fold):
        helix_file = MODEL+"-helix.pdb"
        json_file =  MODEL+".json"       
        out_file = "{m}-{nf:02}/02-final-initial-min.pdb".format(nf=nf, m=MODEL)
        log_file = "{m}-{nf:02}/02-final-initial-min.log".format(nf=nf, m=MODEL)
        name =  "{m}-{nf:02}".format(nf=nf, m=MODEL)          
        cmd_str =  "python {SDIR}/make_initial_model.py --json {json} --helix {helix} \
            --mean 8 --stdev 0.1 \
            --mdsteps 20000 --temp 1200 --shift 2 --out-dir {name} > {log}".format(
                                  m=MODEL, nf=nf, json=json_file, helix=helix_file,
                                  name= name, log=log_file, SDIR=SDIR)       
        yield     {
        'name': name,       
        'file_dep': [helix_file, json_file],
        'targets': [out_file],
        'actions':  [(create_folder, [name]), cmd_str],
        'clean' : True}   


def task_make_homology_models():
    for nf in range(N_fold):
      for nh in range(N_homology):   
        helix_file = MODEL+"-helix.pdb"
        json_file =  MODEL+".json" 
        aln_file =   MODEL+".ali"
        init_file = "{m}-{nf:02}/02-final-initial-min.pdb".format(nf=nf, m=MODEL)        
        out_file = "{m}-{nf:02}/03-homology-model-{nh:02}.pdb".format(nf=nf, nh=nh, m=MODEL) 
        log_file = "{m}-{nf:02}/03-homology-model-{nh:02}.log".format(nf=nf, nh=nh, m=MODEL) 
        name_dir =  "{m}-{nf:02}".format(nf=nf, m=MODEL)          
        name =  "{m}-{nf:02}-{nh:02}".format(nf=nf, nh=nh,  m=MODEL)  
        cmd_str =  "python {SDIR}/make_homology_model.py --json {json} --alnfile {aln} --initial-model {init} \
         --start-index {nh} --end-index {nh}         \
         --repeat 1 --md-level fast --out-dir {name_dir} > {log}".format(
                                  m=MODEL, nf=nf, nh=nh, json=json_file, helix=helix_file, aln=aln_file,
                                  init=init_file, name_dir=name_dir, log=log_file, SDIR=SDIR)       
        yield     {
        'name': name,       
        'file_dep': [init_file, aln_file, json_file],
        'targets': [out_file],
        'actions':  [cmd_str],
        'clean' : True}   
