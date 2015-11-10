name = "APHsh"
entire_sequence ="ELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLKGKGDGELKQLEEELQAIEEQLAQLQWKAQARKEKLAQLK"

segments = """
ELKQLEE ELQAIEE QLAQLQW KAQARKE KLAQLK 		APHshort
GKGDG
ELKQLEE ELQAIEE QLAQLQW KAQARKE KLAQLK 		APHshort
"""

pairs = """
- pair : APHshort:APHshort
  type : P  
  chains : A:B
  template : APH.pdb  
"""   


if __name__ == "__main__":
    import ppmod.make_json as mj   
    pairs = mj.load_pairs(pairs)    
    mj.generate_json(name, entire_sequence, segments, pairs)