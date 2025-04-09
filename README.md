0_get_seq_from_mmcif_s: takes every PDB file provided and extracts the amino acid sequence from the every chain of that PDB file, also checks if all atoms of an amino acid are complete

1_find_crystals_s: uses the edlib library to compare the sequences of all pdb chains obtained by 0_get_seq_from_mmcif_s to all provided gene sequences and gives a similarity score e.g. all human genes

2_structure_data_generator_s: calculates the interfaces in contacting residues and SASA between all chains for every provided PDB .cif data

3_interface_calc_s: takes all genes of interest (here genes detected to participate in co-co assembly) and calculates the total interface for every pdb that matches these genes
