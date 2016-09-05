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
    import cocopod.make_json as mj   
    import cocopod.segment_assignment as sa
    import cocopod.make_color as mc
    import cocopod.utils as u
    entire_sequence = sa.deannotate_sequence(annotated_sequence, remove_whitespace=True)  
    pairs = mj.load_pairs(pairs_info)    

    #generate json file 
    mj.generate_json(model_name, entire_sequence, annotated_sequence, pairs)
    print("Written: "+model_name+".json")    
    
    #create chimera script for coloring.
    mc.chimera_color(model_name+".json", model_name+".chimera", verbose=False)
    print("Written: "+model_name+".chimera")
    
    #write fasta file
    u.write_fasta_file(model_name+".fasta", model_name, entire_sequence)
    print("Written: "+model_name+".fasta")
