import ppmod.make_json as mj
import ppmod.utils as u
import tempfile
import os
    
def test_load_pairs_from_yaml():
    yaml_str="""
        - pair : 
            - APHshort
            - APHshort
          type : A  
          chains : 
            - A
            - B
          template : APH.pdb      
        - pair : APHshort:APHshort
          type : P  
          chains : A:B
          template : APH.pdb     
          
        - pair : "P7 : P8"
          type : P  
          chains : 'A:  B'
          template : APH.pdb              
    """      
    pairs = mj.load_pairs(yaml_str)  
    assert len(pairs)==3
    assert pairs[0]['pair']==['APHshort','APHshort']
    assert pairs[0]['chains']==['A','B']    
    assert pairs[1]['pair']==['APHshort','APHshort']
    assert pairs[1]['chains']==['A','B']    
    assert pairs[2]['pair']==['P7','P8']
    assert pairs[2]['chains']==['A','B']    
   
   
PAIRS = [
    {"pair" : ["APHshort","APHshort"], "type":"A", "chains":["A", "B"], "template":"APH.pdb"},   
    {"pair" : ["P7mS","P8mS"], "type":"P", "chains":["A", "B"], "template":"p7_p8.pdb"}   

    ]   
def test_find_pair_by_segment_name():

    pair = mj.find_pair_by_segment_name(PAIRS,'APHshort')
    assert pair['pair'] == ["APHshort","APHshort"] 
    
    pair = mj.find_pair_by_segment_name(PAIRS,'P7mS')
    assert pair['pair'] == ["P7mS","P8mS"]
    
    pair = mj.find_pair_by_segment_name(PAIRS,'P8mS')
    assert pair['pair'] == ["P7mS","P8mS"]
        
        
def test_get_other_segment_name():

    other_pair = mj.get_other_segment_name(PAIRS[0],'APHshort')
    assert other_pair == "APHshort" 
    
    other_pair = mj.get_other_segment_name(PAIRS[1],'P7mS')
    assert other_pair == "P8mS"
    
    other_pair = mj.get_other_segment_name(PAIRS[1],'P8mS')
    assert other_pair == "P7mS"



def test_get_segment_chain():
    assert mj.get_segment_chain(PAIRS[0], "APHshort") == "A"
    assert mj.get_segment_chain(PAIRS[1], "P7mS") == "A"
    assert mj.get_segment_chain(PAIRS[1], "P8mS") == "B"
    
def test_generate_json_tet_p5LD3():
    name = "test"    
    entire_sequence = "MMMSPEDEIQSLEEKNSQLKQEISQLEEKNQQLKYGELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLKGKGDGELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLKGKGDGSPEDKISQLKEENQQLEQKIQQLKEENSQLEYGGDGKGLEHHHHH"
    segments_str = """
                SPED EIQSLEE KNSQLKQ EISQLEE KNQQLKY G 			P7mS
                GKGDG
                ELKQLEE ELQAIEE QLAQLQW KAQARKE KLAQLK 		APHshort
                GKGDG
                ELKQLEE ELQAIEE QLAQLQW KAQARKE KLAQLK 		APHshort
                GKGDG
                SPED KISQLKE ENQQLEQ KIQQLKE ENSQLEY G 			P8mS
                GDGKG
                """
    pairs = [
    {"pair" : ["APHshort","APHshort"], "type":"A", "chains":["A", "B"], "template":"APH.pdb"},   
    {"pair" : ["P7mS","P8mS"], "type":"P", "chains":["A", "B"], "template":"p7_p8.pdb"}   

        ]            
    
    
    out_name = tempfile.NamedTemporaryFile().name    
    mj.generate_json(name, entire_sequence, segments_str, pairs, out_name=out_name)
    

    d = u.load_json_data(out_name)    
    assert d.name == name
    assert d.entire_sequence == entire_sequence  
    assert len(d.segments) == 4
    assert d.segments[0].id == 1
    assert d.segments[0].start == 4
    assert d.segments[0].end == 4+33-1
    assert d.segments[0].pair_id == 4
    assert d.segments[0].pair_type == "P"
    assert d.segments[0].pdb_template == "p7_p8.pdb"
    assert d.segments[0].pdb_chain == "A"
    
    assert d.segments[2].id == 3
    assert d.segments[2].start == 76
    assert d.segments[2].end == 109
    assert d.segments[2].pair_id == 2
    assert d.segments[2].pair_type == "A"
    assert d.segments[2].pdb_template == "APH.pdb"
    assert d.segments[2].pdb_chain == "B"   
    
    os.remove(out_name)
    