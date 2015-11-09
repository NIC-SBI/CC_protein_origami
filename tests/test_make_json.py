import ppmod
import ppmod.make_json as mj
    
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
          
        - pair : "APHshort : APHshort"
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
    assert pairs[2]['pair']==['APHshort','APHshort']
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
    
def test_generate_json_tet_p5LD3():
    name = "test"    
    entire_sequence = "MRMKQLEDKVEELERKNYHLENEVSRLKKLVGGDGKGSPEDEIQSLEEKNSQLKQEISQLEEKNQQLKYGGKGDGELERAKQSIRRLEQEVNQERSRMQYLQTGKGDGRMKQLEDKVEELESKNYHLENEVSRLKKLVGGKGDGELKQLEE ELQAIEEQLAQLQWKAQARKEKLAQLKGKGDGSPEDEIQQLEEEISQLEQKNSQLKEKNQQLKYGGKGDGELERAKQSIRRLEQEVNQERSRMQYLQTGKGDGSPEDENSQLEEKISQLKQKNSQLKEEIQQLEYGGKGDGELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLKGKGDGSPEDKISQLKEENQQLEQKIQQLKEENSQLEYGGDGKGSPEDKNSQLKEEIQQLEEENQQLEEKISQLKYGGDGKGSPEDKISQLKQKIQQLKQENQQLEEENSQLEYGLEHHHHHHHH"
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

    mj.generate_json(name, entire_sequence, segments_str, pairs)
    