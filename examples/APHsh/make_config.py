model_name = "APHsh"
annotated_sequence = """
ELKQLEE ELQAIEE QLAQLQW KAQARKE KLAQLK 		APHshort
GKGDG
ELKQLEE ELQAIEE QLAQLQW KAQARKE KLAQLK 		APHshort
"""
pairs_info = """
- pair : APHshort:APHshort
  type : A  
  chains : A:B
  template : APH.pdb  
"""   

if __name__ == "__main__":
    import ppmod.make_json as mj   
    import ppmod.segment_assignment as sa
    entire_sequence = sa.deannotate_sequence(annotated_sequence, remove_whitespace=True)  
    pairs = mj.load_pairs(pairs_info)    
    mj.generate_json(model_name, entire_sequence, annotated_sequence, pairs)
