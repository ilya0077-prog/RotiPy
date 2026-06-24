from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
import numpy as np

class RotiPy:
    def __init__(self, 
                 peptide_A='rotifer_data/experimental_group/rotifer_sco_spondin_derived_hexapeptide/fold_rotifer_hexapeptide_model_0.cif', 
                 peptide_B='rotifer_data/control_group/rotifer_poly_alanine_homopolymer/fold_rotifer_hexapeptide_control_0_model_0.cif', 
                 peptide_C='rotifer_data/control_group/rotifer_scrambled_hexapeptide_SSDLDN/fold_2026_06_04_22_47_model_0.cif', 
                 protein_B='prion_data/experimental_group_prion_fibril/7UMQ.cif', 
                 protein_C='prion_data/control_group_human_prion_protein/1QLX.cif', 
                 LIMIT=6):
        self.peptide_A=peptide_A
        self.peptide_B=peptide_B
        self.peptide_C=peptide_C
        self.protein_B=protein_B
        self.protein_C=protein_C
        self.LIMIT=LIMIT
    
    def parse_alphafold_data(self):

        cif_parser=MMCIFParser(QUIET=True)

        structure_A=cif_parser.get_structure('peptide_A',
        self.peptide_A)
        structure_B=cif_parser.get_structure('peptide_B',
        self.peptide_B)
        structure_C=cif_parser.get_structure('peptide_C',
        self.peptide_C)
        structure_7UMQ=cif_parser.get_structure('7UMQ',
        self.protein_B)
        structure_1QLX=cif_parser.get_structure('1QLX',
        self.protein_C)

        self.atom_list_peptide_A=list(structure_A.get_atoms())
        self.atom_list_peptide_B=list(structure_B.get_atoms())
        self.atom_list_peptide_C=list(structure_C.get_atoms())
        self.atom_list_7UMQ=list(structure_7UMQ.get_atoms())
        self.atom_list_1QLX=list(structure_1QLX.get_atoms())

        pdb_parser=PDBParser(QUIET=True)

        mdockpep_res=pdb_parser.get_structure('peptide_A_7UMQ','mdockpep_data/dssndl_7umq.pdb')

        self.atom_list_mdock=list(mdockpep_res.get_atoms())

    def calculate_biopython_residues(self):

        chains=list(set([atom.get_parent().get_parent().id for atom in self.atom_list_mdock]))

        atoms_peptide=[atom for atom in self.atom_list_mdock if atom.get_parent().get_parent().id == ' ']
        atoms_target=[atom for atom in self.atom_list_mdock if atom.get_parent().get_parent().id != ' ']

        self.neighbor_search_target=NeighborSearch(atoms_target)
        cutoff_distance=3.5

        target_atoms=[]
        for atom in atoms_peptide:
            coord=atom.get_coord()
            range=self.neighbor_search_target.search(coord, cutoff_distance, level='A')
            target_atoms.extend(range)

        self.search_target=list(set(target_atoms))

        peptide_search=NeighborSearch(atoms_peptide)
        
        self.search_peptide=[]
        for target_atom in self.search_target:
            close_pep_atoms=peptide_search.search(target_atom.get_coord(), cutoff_distance, level='A')
            self.search_peptide.extend(close_pep_atoms)

        self.search_peptide=list(set(self.search_peptide))

        xyz_peptide=np.array([atom.get_coord() for atom in self.search_peptide])
        xyz_protein=np.array([atom.get_coord() for atom in self.search_target])

        peptide_res_num=list(set([atom.get_parent().get_id()[1] for atom in self.search_peptide]))
        protein_res_num=list(set([atom.get_parent().get_id()[1] for atom in self.search_target]))

        print(f'peptide binding pocket residues -> {peptide_res_num}, target binding pocket residues -> {protein_res_num}')

    @property
    def LIMIT(self):
        return self._LIMIT

    @LIMIT.setter
    def LIMIT(self, LIMIT):
        if LIMIT !=6:
            raise ValueError
        self._LIMIT=LIMIT

        
