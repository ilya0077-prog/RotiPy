from RotiPy import RotiPy

if __name__=="__main__":
    
    run = RotiPy()
    run.parse_alphafold_data()
    run.calculate_biopython_residues()