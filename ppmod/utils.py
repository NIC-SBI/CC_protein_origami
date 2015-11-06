#taken from http://hayd.github.io/2013/dotable-dictionaries/
import os

class Dotable(dict):
    """Generates nested dotable dicts from a json-like object. This makes is possible to write
    a.b[0].c
    """
    __getattr__= dict.__getitem__

    def __init__(self, d):
        self.update(**dict((k, self.parse(v))
                           for k, v in d.items()))

    @classmethod
    def parse(cls, v):
        if isinstance(v, dict):
            return cls(v)
        elif isinstance(v, list):
            return [cls.parse(i) for i in v]
        else:
            return v
        
        
#taken from http://stackoverflow.com/a/13105359/952600
def byteify(input):
    """Takes a jason-like structure and convets unicode str to str"""
    if isinstance(input, dict):
        return {byteify(key):byteify(value) for key,value in input.items()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input        
    
def load_json_data(file_name):
    """Loads data from json and returns a dotable dict with ascii strings"""
    import json
    with open(file_name, 'r') as infile:
        data=json.load(infile)
#    data=byteify(data)
    return Dotable.parse(data)

def pair_ids_from_segments(segments):
    """Return the pair IDs of the segments. The IDs are ZERO based. The matching is based on naming"""
    pair_ids = []
    for n in range(len(segments)):
        for n1 in range(n+1,len(segments)):
            if segments[n]['name']==segments[n1]['pair_name']:
                pair_ids.append((n, n1))
    return pair_ids    
    

import string
import random
import numpy as np
import mdtraj as md
#taken from http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))    


def sequnce_and_knowns(alnfile):
    """Parses the IDs of the known structures and the sequences from the alignment file.
    Returns a tuple:
      (seq, (known1, known2 ...))

    """
    knowns = ()
    with open(alnfile,'r') as f:
        while True:
            line = f.readline()
            if not line: break
            if line[:4] == '>P1;':
                seqname = line[4:].rstrip()
                line = f.readline()            
                if line.split(':')[0] == 'sequence':
                    sequence = seqname
                else:
                     knowns = knowns+(seqname,)
    return (sequence, knowns)                 

def relative_to(file_dir, a_path):
    """Returns the directory of file_dir and appends a_path"""
    #os.path.realpath(
    #os.path.join(os.getcwd(), os.path.dirname(__file__)))
    a_dir = os.path.dirname(file_dir)
    return os.path.join(a_dir,a_path)

def generate_json(name, entire_sequence, segments_str, pairs, out_name=None):
    """generates json file with information about the polypeptide polyhedra.
    Segment parings are calculated, etc...

    Parameters
    ----------
    name : str
        The name of the polyhedral design.
    entire_sequence : str
        The whole sequnce, including his tags etc..
    segments_str: str
        A string of segment sequnces and segment names. The names must be tab separated.
    pairs : dict
        

    Returns
    -------
      Does not return anything. By default generates an out_name file
    """  
    pass    
    
def align(p1f, seq, j, jj):
    """Align the template sequence to the target sequence

    Parameters
    ----------
    p1f : str
        Template sequence string
    seq : str
        Target sequence string
    j : int
        How many different template sequences are checked for optimum alignment (different sequences are obtained by sequentially shortening p1f by one AA)  
    jj : int  
        How many different target sequences are checked for optimum alignment (different sequences are obtained by sequentially shortening seq by one AA)

    Returns
    -------
    i : int
        Start of the alignment on template sequence 
    ii : int
        Start of the alignment on target sequence
    """    
    i = 0
    ii = 0
    scoreold = score(p1f, seq)
    for t in range(j):   
       for tt in range(jj): 
           if t == tt: continue
           scorenew = score(p1f[t:], seq[tt:])  
           if scorenew > scoreold:
               scoreold = scorenew               
               i = t
               ii = tt
    return (i, ii)

def score(p1f, seq):
    """This function scores the quality of the alignment. Higher score means better alignment

    Parameters
    ----------
    p1f : str
        Template sequence string
    seq : str
        Target sequence string
    
    Returns
    -------
    total : int
        Final score of the alignmnet
    """
    total = 0
    for i in range(8):  #check the first ten residues
        if p1f[i] == seq[i]:   #if p1f and seq have a matching residue on position i add a point
            total = total+1  
    return total                          #return final score

def find_pair(pair, segments):
    """Find pair's id, pair_id, and its pdb file

    Parameters
    ----------
    pair : str
        Name of CC segments
    
    Returns
    -------
    p1 : int
        Id number of the CC segment 
    p2 : int
        Pair_id number of the CC segment
    pdbname : str
        Name of the pdb file where the CC segment is written
    """
    for s in segments:
        if s['name'] == pair:
            p1 = s['id']-1         #get segment name
            p2 = s['pair_id']-1
            pdbname = s['pdb_template']   #get template pdb file name
            break     
    return p1, p2, pdbname

def selres(i, t):
    """Select i-th residue from topology

    Parameters
    ----------
    i : int
        Residue index
    t : mdtraj topology
        Topology
    
    Returns
    -------
    res : mdtraj residue
        i-th residue
    """
    res = t.select("resid {:d}".format(i))
    return res

def writepdb(i, l, t, p, path):
    """Write the aligned part of the template sequence to pdb

    Parameters
    ----------
    i : int
        Start of the aligning part on the template sequence
    l : int
        Length of the aligned template sequence
    t : mdtraj topology
        topology of the template 
    p : md traj positions
        atom positions in template    
    """
    cl0 = t.chain(0).n_residues  #cl0 length of the first chain of a CC pair
    topsub1 = t.subset(list(range(selres(i, t)[0],selres(i + l, t)[-1]))) #topology of atoms in the first aligned chain of the CC segment
    topsub2 = t.subset(list(range(selres(cl0 + i, t)[0], selres(cl0 +i + l, t)[-1]))) #topology of atoms in the first aligned chain of the CC segment
    topsubjoin = topsub1.join(topsub2) # join topology subsets
    coord1 = p[0, selres(i, t)[0]:selres(i + l, t)[-1], :]*10 #position of atoms in the first subset
    coord2 = p[0, selres(cl0 + i, t)[0]:selres(cl0 + i + l, t)[-1], :]*10 #position of atoms in the second subset
    coord = np.concatenate((coord1, coord2), axis=0) #join position subsets
    fpdb = md.formats.PDBTrajectoryFile(path, mode='w') 
    fpdb.write(coord, topsubjoin) #write to pdb
    return



_AMINO_ACID_CODES =  {'ACE': None, 'NME':  None, '00C': 'C', '01W':  'X', '02K':
'A', '02L':  'N', '03Y': 'C',  '07O': 'C', '08P':  'C', '0A0': 'D',  '0A1': 'Y',
'0A2': 'K', '0A8':  'C', '0AA': 'V', '0AB': 'V', '0AC':  'G', '0AF': 'W', '0AG':
'L', '0AH':  'S', '0AK': 'D',  '0BN': 'F', '0CS':  'A', '0E5': 'T',  '0EA': 'Y',
'0FL': 'A', '0NC':  'A', '0WZ': 'Y', '0Y8': 'P', '143':  'C', '193': 'X', '1OP':
'Y', '1PA':  'F', '1PI': 'A',  '1TQ': 'W', '1TY':  'Y', '1X6': 'S',  '200': 'F',
'23F': 'F', '23S':  'X', '26B': 'T', '2AD': 'X', '2AG':  'A', '2AO': 'X', '2AS':
'X', '2CO':  'C', '2DO': 'X',  '2FM': 'M', '2HF':  'H', '2KK': 'K',  '2KP': 'K',
'2LU': 'L', '2ML':  'L', '2MR': 'R', '2MT': 'P', '2OR':  'R', '2PI': 'X', '2QZ':
'T', '2R3':  'Y', '2SI': 'X',  '2TL': 'T', '2TY':  'Y', '2VA': 'V',  '2XA': 'C',
'32S': 'X', '32T':  'X', '33X': 'A', '3AH': 'H', '3AR':  'X', '3CF': 'F', '3GA':
'A', '3MD':  'D', '3NF': 'Y',  '3QN': 'K', '3TY':  'X', '3XH': 'G',  '4BF': 'Y',
'4CF': 'F', '4CY':  'M', '4DP': 'W', '4FB': 'P', '4FW':  'W', '4HT': 'W', '4IN':
'W', '4MM':  'X', '4PH': 'F',  '4U7': 'A', '56A':  'H', '5AB': 'A',  '5CS': 'C',
'5CW': 'W', '5HP':  'E', '6CL': 'K', '6CW': 'W', '6GL':  'A', '6HN': 'K', '7JA':
'I', '9NE':  'E', '9NF': 'F',  '9NR': 'R', '9NV':  'V', 'A5N': 'N',  'A66': 'X',
'AA3': 'A', 'AA4':  'A', 'AAR': 'R', 'AB7': 'X', 'ABA':  'A', 'ACB': 'D', 'ACL':
'R', 'ADD':  'X', 'AEA': 'X',  'AEI': 'D', 'AFA':  'N', 'AGM': 'R',  'AGT': 'C',
'AHB': 'N', 'AHH':  'X', 'AHO': 'A', 'AHP': 'A', 'AHS':  'X', 'AHT': 'X', 'AIB':
'A', 'AKL':  'D', 'AKZ': 'D',  'ALA': 'A', 'ALC':  'A', 'ALM': 'A',  'ALN': 'A',
'ALO': 'T', 'ALS':  'A', 'ALT': 'A', 'ALV': 'A', 'ALY':  'K', 'AN8': 'A', 'APE':
'X', 'APH':  'A', 'API': 'K',  'APK': 'K', 'APM':  'X', 'APP': 'X',  'AR2': 'R',
'AR4': 'E', 'AR7':  'R', 'ARG': 'R', 'ARM': 'R', 'ARO':  'R', 'ARV': 'X', 'AS2':
'D', 'AS9':  'X', 'ASA': 'D',  'ASB': 'D', 'ASI':  'D', 'ASK': 'D',  'ASL': 'D',
'ASM': 'X', 'ASN':  'N', 'ASP': 'D', 'ASQ': 'D', 'ASX':  'B', 'AVN': 'X', 'AYA':
'A', 'AZK':  'K', 'AZS': 'S',  'AZY': 'Y', 'B1F':  'F', 'B2A': 'A',  'B2F': 'F',
'B2I': 'I', 'B2V':  'V', 'B3A': 'A', 'B3D': 'D', 'B3E':  'E', 'B3K': 'K', 'B3L':
'X', 'B3M':  'X', 'B3Q': 'X',  'B3S': 'S', 'B3T':  'X', 'B3U': 'H',  'B3X': 'N',
'B3Y': 'Y', 'BB6':  'C', 'BB7': 'C', 'BB8': 'F', 'BB9':  'C', 'BBC': 'C', 'BCS':
'C', 'BE2':  'X', 'BFD': 'D',  'BG1': 'S', 'BH2':  'D', 'BHD': 'D',  'BIF': 'F',
'BIL': 'X', 'BIU':  'I', 'BJH': 'X', 'BL2': 'L', 'BLE':  'L', 'BLY': 'K', 'BMT':
'T', 'BNN':  'F', 'BNO': 'X',  'BOR': 'R', 'BPE':  'C', 'BSE': 'S',  'BTA': 'L',
'BTC': 'C', 'BTR':  'W', 'BUC': 'C', 'BUG': 'V', 'C1X':  'K', 'C22': 'A', 'C3Y':
'C', 'C4R':  'C', 'C5C': 'C',  'C66': 'X', 'C6C':  'C', 'CAF': 'C',  'CAL': 'X',
'CAS': 'C', 'CAV':  'X', 'CAY': 'C', 'CCL': 'K', 'CCS':  'C', 'CDE': 'X', 'CDV':
'X', 'CEA':  'C', 'CGA': 'E',  'CGU': 'E', 'CHF':  'X', 'CHG': 'X',  'CHP': 'G',
'CHS': 'X', 'CIR':  'R', 'CLE': 'L', 'CLG': 'K', 'CLH':  'K', 'CME': 'C', 'CMH':
'C', 'CML':  'C', 'CMT': 'C',  'CPC': 'X', 'CPI':  'X', 'CR5': 'G',  'CS0': 'C',
'CS1': 'C', 'CS3':  'C', 'CS4': 'C', 'CSA': 'C', 'CSB':  'C', 'CSD': 'C', 'CSE':
'C', 'CSJ':  'C', 'CSO': 'C',  'CSP': 'C', 'CSR':  'C', 'CSS': 'C',  'CSU': 'C',
'CSW': 'C', 'CSX':  'C', 'CSZ': 'C', 'CTE': 'W', 'CTH':  'T', 'CUC': 'X', 'CWR':
'S', 'CXM':  'M', 'CY0': 'C',  'CY1': 'C', 'CY3':  'C', 'CY4': 'C',  'CYA': 'C',
'CYD': 'C', 'CYF':  'C', 'CYG': 'C', 'CYJ': 'K', 'CYM':  'C', 'CYQ': 'C', 'CYR':
'C', 'CYS':  'C', 'CZ2': 'C',  'CZZ': 'C', 'D11':  'T', 'D3P': 'G',  'D4P': 'X',
'DA2': 'X', 'DAB':  'A', 'DAH': 'F', 'DAL': 'A', 'DAR':  'R', 'DAS': 'D', 'DBB':
'T', 'DBS':  'S', 'DBU': 'T',  'DBY': 'Y', 'DBZ':  'A', 'DC2': 'C',  'DCL': 'X',
'DCY': 'C', 'DDE':  'H', 'DFI': 'X', 'DFO': 'X', 'DGH':  'G', 'DGL': 'E', 'DGN':
'Q', 'DHA':  'S', 'DHI': 'H',  'DHL': 'X', 'DHN':  'V', 'DHP': 'X',  'DHV': 'V',
'DI7': 'Y', 'DIL':  'I', 'DIR': 'R', 'DIV': 'V', 'DLE':  'L', 'DLS': 'K', 'DLY':
'K', 'DM0':  'K', 'DMH': 'N',  'DMK': 'D', 'DMT':  'X', 'DNE': 'L',  'DNL': 'K',
'DNP': 'A', 'DNS':  'K', 'DOA': 'X', 'DOH': 'D', 'DON':  'L', 'DPL': 'P', 'DPN':
'F', 'DPP':  'A', 'DPQ': 'Y',  'DPR': 'P', 'DSE':  'S', 'DSG': 'N',  'DSN': 'S',
'DSP': 'D', 'DTH':  'T', 'DTR': 'W', 'DTY': 'Y', 'DVA':  'V', 'DYS': 'C', 'ECC':
'Q', 'EFC':  'C', 'EHP': 'F',  'ESB': 'Y', 'ESC':  'M', 'EXY': 'L',  'EYS': 'X',
'F2F': 'F', 'FAK':  'K', 'FB5': 'A', 'FB6': 'A', 'FCL':  'F', 'FGA': 'E', 'FGL':
'G', 'FGP':  'S', 'FH7': 'K',  'FHL': 'K', 'FHO':  'K', 'FLA': 'A',  'FLE': 'L',
'FLT': 'Y', 'FME':  'M', 'FOE': 'C', 'FP9': 'P', 'FRD':  'X', 'FT6': 'W', 'FTR':
'W', 'FTY':  'Y', 'FVA': 'V',  'FZN': 'K', 'GAU':  'E', 'GCM': 'X',  'GFT': 'S',
'GGL': 'E', 'GHG':  'Q', 'GHP': 'G', 'GL3': 'G', 'GLH':  'Q', 'GLJ': 'E', 'GLK':
'E', 'GLM':  'X', 'GLN': 'Q',  'GLQ': 'E', 'GLU':  'E', 'GLX': 'Z',  'GLY': 'G',
'GLZ': 'G', 'GMA':  'E', 'GND': 'X', 'GPL': 'K', 'GSC':  'G', 'GSU': 'E', 'GT9':
'C', 'GVL':  'S', 'H14': 'F',  'H5M': 'P', 'HAC':  'A', 'HAR': 'R',  'HBN': 'H',
'HCS': 'X', 'HFA':  'X', 'HGL': 'X', 'HHI': 'H', 'HIA':  'H', 'HIC': 'H', 'HIP':
'H', 'HIQ':  'H', 'HIS': 'H',  'HL2': 'L', 'HLU':  'L', 'HMR': 'R',  'HPC': 'F',
'HPE': 'F', 'HPH':  'F', 'HPQ': 'F', 'HQA': 'A', 'HRG':  'R', 'HRP': 'W', 'HS8':
'H', 'HS9':  'H', 'HSE': 'S',  'HSL': 'S', 'HSO':  'H', 'HTI': 'C',  'HTN': 'N',
'HTR': 'W', 'HV5':  'A', 'HVA': 'V', 'HY3': 'P', 'HYP':  'P', 'HZP': 'P', 'I2M':
'I', 'I58':  'K', 'IAM': 'A',  'IAR': 'R', 'IAS':  'D', 'IEL': 'K',  'IGL': 'G',
'IIL': 'I', 'ILE':  'I', 'ILG': 'E', 'ILX': 'I', 'IML':  'I', 'IOY': 'F', 'IPG':
'G', 'IT1':  'K', 'IYR': 'Y',  'IYT': 'T', 'IZO':  'M', 'JJJ': 'C',  'JJK': 'C',
'JJL': 'C', 'K1R':  'C', 'KCX': 'K', 'KGC': 'K', 'KNB':  'A', 'KOR': 'M', 'KPI':
'K', 'KST':  'K', 'KYN': 'W',  'KYQ': 'K', 'L2A':  'X', 'LA2': 'K',  'LAA': 'D',
'LAL': 'A', 'LBY':  'K', 'LCK': 'K', 'LCX': 'K', 'LCZ':  'X', 'LDH': 'K', 'LED':
'L', 'LEF':  'L', 'LEH': 'L',  'LEI': 'V', 'LEM':  'L', 'LEN': 'L',  'LET': 'K',
'LEU': 'L', 'LEX':  'L', 'LHC': 'X', 'LLP': 'K', 'LLY':  'K', 'LME': 'E', 'LMF':
'K', 'LMQ':  'Q', 'LP6': 'K',  'LPD': 'P', 'LPG':  'G', 'LPL': 'X',  'LPS': 'S',
'LSO': 'K', 'LTA':  'X', 'LTR': 'W', 'LVG': 'G', 'LVN':  'V', 'LYF': 'K', 'LYK':
'K', 'LYM':  'K', 'LYN': 'K',  'LYR': 'K', 'LYS':  'K', 'LYX': 'K',  'LYZ': 'K',
'M0H': 'C',  'M2L': 'K', 'M2S': 'M',  'M30': 'G', 'M3L': 'K',  'MA': 'A', 'MAA':
'A', 'MAI':  'R', 'MBQ': 'Y',  'MC1': 'S', 'MCG':  'X', 'MCL': 'K',  'MCS': 'C',
'MD3': 'C', 'MD6':  'G', 'MDF': 'Y', 'MDH': 'X', 'MEA':  'F', 'MED': 'M', 'MEG':
'E', 'MEN':  'N', 'MEQ': 'Q',  'MET': 'M', 'MEU':  'G', 'MF3': 'X',  'MGG': 'R',
'MGN': 'Q', 'MGY':  'G', 'MHL': 'L', 'MHO': 'M', 'MHS':  'H', 'MIS': 'S', 'MK8':
'L', 'ML3':  'K', 'MLE': 'L',  'MLL': 'L', 'MLY':  'K', 'MLZ': 'K',  'MME': 'M',
'MMO': 'R', 'MND':  'N', 'MNL': 'L', 'MNV': 'V', 'MOD':  'X', 'MP8': 'P', 'MPH':
'X', 'MPJ':  'X', 'MPQ': 'G',  'MSA': 'G', 'MSE':  'M', 'MSL': 'M',  'MSO': 'M',
'MSP': 'X', 'MT2':  'M', 'MTY': 'Y', 'MVA': 'V', 'N10':  'S', 'N2C': 'X', 'N7P':
'P', 'N80':  'P', 'N8P': 'P',  'NA8': 'A', 'NAL':  'A', 'NAM': 'A',  'NB8': 'N',
'NBQ': 'Y', 'NC1':  'S', 'NCB': 'A', 'NCY': 'X', 'NDF':  'F', 'NEM': 'H', 'NEP':
'H', 'NFA':  'F', 'NHL': 'E',  'NIY': 'Y', 'NLE':  'L', 'NLN': 'L',  'NLO': 'L',
'NLP': 'L', 'NLQ':  'Q', 'NMC': 'G', 'NMM': 'R', 'NNH':  'R', 'NPH': 'C', 'NPI':
'A', 'NSK':  'X', 'NTR': 'Y',  'NTY': 'Y', 'NVA':  'V', 'NYS': 'C',  'NZH': 'H',
'O12': 'X', 'OAR':  'R', 'OAS': 'S', 'OBF': 'X', 'OBS':  'K', 'OCS': 'C', 'OCY':
'C', 'OHI':  'H', 'OHS': 'D',  'OIC': 'X', 'OLE':  'X', 'OLT': 'T',  'OLZ': 'S',
'OMT': 'M', 'ONH':  'A', 'ONL': 'X', 'OPR': 'R', 'ORN':  'A', 'ORQ': 'R', 'OSE':
'S', 'OTB':  'X', 'OTH': 'T',  'OXX': 'D', 'P1L':  'C', 'P2Y': 'P',  'PAQ': 'Y',
'PAS': 'D', 'PAT':  'W', 'PAU': 'A', 'PBB': 'C', 'PBF':  'F', 'PCA': 'E', 'PCC':
'P', 'PCE':  'X', 'PCS': 'F',  'PDL': 'X', 'PEC':  'C', 'PF5': 'F',  'PFF': 'F',
'PFX': 'X', 'PG1':  'S', 'PG9': 'G', 'PGL': 'X', 'PGY':  'G', 'PH6': 'P', 'PHA':
'F', 'PHD':  'D', 'PHE': 'F',  'PHI': 'F', 'PHL':  'F', 'PHM': 'F',  'PIV': 'X',
'PLE': 'L', 'PM3':  'F', 'POM': 'P', 'PPN': 'F', 'PR3':  'C', 'PR9': 'P', 'PRO':
'P', 'PRS':  'P', 'PSA': 'F',  'PSH': 'H', 'PTA':  'X', 'PTH': 'Y',  'PTM': 'Y',
'PTR': 'Y', 'PVH':  'H', 'PVL': 'X', 'PYA': 'A', 'PYL':  'O', 'PYX': 'C', 'QCS':
'C', 'QMM':  'Q', 'QPA': 'C',  'QPH': 'F', 'R1A':  'C', 'R4K': 'W',  'RE0': 'W',
'RE3': 'W', 'RON':  'X', 'RVX': 'S', 'RZ4': 'S', 'S1H':  'S', 'S2C': 'C', 'S2D':
'A', 'S2P':  'A', 'SAC': 'S',  'SAH': 'C', 'SAR':  'G', 'SBL': 'S',  'SCH': 'C',
'SCS': 'C', 'SCY':  'C', 'SD2': 'X', 'SDP': 'S', 'SE7':  'A', 'SEB': 'S', 'SEC':
'U', 'SEG':  'A', 'SEL': 'S',  'SEM': 'S', 'SEN':  'S', 'SEP': 'S',  'SER': 'S',
'SET': 'S', 'SGB':  'S', 'SHC': 'C', 'SHP': 'G', 'SHR':  'K', 'SIB': 'C', 'SLR':
'P', 'SLZ':  'K', 'SMC': 'C',  'SME': 'M', 'SMF':  'F', 'SNC': 'C',  'SNN': 'N',
'SOC': 'C', 'SOY':  'S', 'SRZ': 'S', 'STY': 'Y', 'SUB':  'X', 'SUN': 'S', 'SVA':
'S', 'SVV':  'S', 'SVW': 'S',  'SVX': 'S', 'SVY':  'S', 'SVZ': 'S',  'SYS': 'C',
'T11': 'F', 'T66':  'X', 'TA4': 'X', 'TAV': 'D', 'TBG':  'V', 'TBM': 'T', 'TCQ':
'Y', 'TCR':  'W', 'TDD': 'L',  'TFQ': 'F', 'TH6':  'T', 'THC': 'T',  'THO': 'X',
'THR': 'T', 'THZ':  'R', 'TIH': 'A', 'TMB': 'T', 'TMD':  'T', 'TNB': 'C', 'TNR':
'S', 'TOQ':  'W', 'TPH': 'X',  'TPL': 'W', 'TPO':  'T', 'TPQ': 'Y',  'TQI': 'W',
'TQQ': 'W', 'TRF':  'W', 'TRG': 'K', 'TRN': 'W', 'TRO':  'W', 'TRP': 'W', 'TRQ':
'W', 'TRW':  'W', 'TRX': 'W',  'TRY': 'W', 'TST':  'X', 'TTQ': 'W',  'TTS': 'Y',
'TXY': 'Y', 'TY1':  'Y', 'TY2': 'Y', 'TY3': 'Y', 'TY5':  'Y', 'TYB': 'Y', 'TYI':
'Y', 'TYJ':  'Y', 'TYN': 'Y',  'TYO': 'Y', 'TYQ':  'Y', 'TYR': 'Y',  'TYS': 'Y',
'TYT': 'Y', 'TYW':  'Y', 'TYX': 'X', 'TYY': 'Y', 'TZB':  'X', 'TZO': 'X', 'UMA':
'A', 'UN1':  'X', 'UN2': 'X',  'UNK': 'X', 'VAD':  'V', 'VAF': 'V',  'VAL': 'V',
'VB1': 'K', 'VDL':  'X', 'VLL': 'X', 'VLM': 'X', 'VMS':  'X', 'VOL': 'X', 'WLU':
'L', 'WPA':  'F', 'WRP': 'W',  'WVL': 'V', 'X2W':  'E', 'XCN': 'C',  'XCP': 'X',
'XDT': 'T', 'XPL':  'O', 'XPR': 'P', 'XSN': 'N', 'XX1':  'K', 'YCM': 'C', 'YOF':
'Y', 'YTH':  'T', 'Z01': 'A',  'ZAL': 'A', 'ZCL':  'F', 'ZFB': 'X',  'ZU0': 'T',
'ZZJ': 'A'}        

def mdtraj_to_fasta(topology, chain=None):
    """Convert this topology into FASTA string
    Parameters
    ----------
    chain : Integer, optional, default=None
        If specified, will return the FASTA string for this chain in the
        Topology.
    Returns
    -------
    fasta : String or list of Strings
       A FASTA string for each chain specified.
    """    
    fasta = lambda c: "".join([_AMINO_ACID_CODES[res.name] for res in c.residues
                               if res.is_protein and res.name is not None])
    if chain is not None:
        if not isinstance(chain, int):
            raise ValueError('chain must be an Integer.')
        return fasta(topology._chains[chain])
    else:
        return [fasta(c) for c in topology._chains]
        
