import numpy as np
import pickle
import matplotlib.pyplot as plt
import itertools

####parameter####
seq_similarity = 0.9

#load gene-pdb database
with open("/path/to/humans_all.pick", "rb") as handle:
    pdb_db = pickle.load(handle)

not_in_pdb = 0
not_two_times_in_structure = 0
no_structure_before_onset = 0
bad_seq_fit = 0
genes_c, crystals_c, chains_c, HCc, LCc = [], [], 0, [], []
cc, cp, pc = 0, 0, 0
out_ = open("/path/to/0_total_interfaces.txt", 'w')
out_.write('gene' + '\t' + 'confidence' + '\t' + 'pdb_structure' + '\t' + 'interacting_chains' + '\t' + \
            'contacts_chain_1_to_chain_2' + '\t' + 'contacts_chain_2_to_chain_1' + '\t' + \
            'intramolecular_contacts_of_interface_residues_of_chain_1' + '\t' + \
            'intramolecular_contacts_of_interface_residues_of_chain_2' + '\n')

line_count = 0
with open("/path/to/Bertolini2020_TableS1.txt") as infile:
    for _ in range(2):
        next(infile)
    for line in infile:
        line_count += 1      
        if line_count%100 == 0:
            print(line_count) 
        A = line.split("\t")
        gene = A[0]
        if gene != 'LMNA':
            continue
        onset_RP = int(A[1])+1e7
        onset = onset_RP-30

        if gene not in pdb_db:
            not_in_pdb += 1
            continue
        if A[2] == 'TRUE':
            confidence = 'high'
        else:
            confidence = 'low'
           
        structures = pdb_db[gene]
        for pdb_code, alignments in structures.items():
            if len(alignments) < 2:
                not_two_times_in_structure += 1
                continue
            chains_to_compare = []
            for chain, alignment in alignments.items():
                if alignment[1] > onset:
                    no_structure_before_onset += 1
                    continue
                if alignment[0] <= seq_similarity:
                    bad_seq_fit += 1
                    continue
                chains_to_compare.append(chain)
            if len(chains_to_compare) < 2:
                continue
            try:
                with open("/path/to/assemblies/"+pdb_code[1:3]+"/"+pdb_code+"_contacts_inter.pick", "rb") as handle:
                    pdb_interface = pickle.load(handle)
                with open("/path/to/assemblies/"+pdb_code[1:3]+"/"+pdb_code+"_contacts_intra.pick", "rb") as handle:
                    pdb_contacts = pickle.load(handle)
            except:
                continue
            for i in itertools.permutations(chains_to_compare, 2):
                chain_combination = i[0]+'_'+i[1]
                if chain_combination in pdb_interface:
                    interface = pdb_interface[chain_combination]
                    interface_res_A_ids = np.where(interface[:, 0] <= onset)[0]
                    interface_res_A = interface[interface_res_A_ids, 0]
                    interface_contacts_A = np.sum(interface[interface_res_A_ids, 2])

                    interface_res_B_ids = np.where(interface[:, 1] <= onset)[0]
                    interface_res_B = interface[interface_res_B_ids, 1]
                    interface_contacts_B = np.sum(interface[interface_res_B_ids, 2])

                    intra_res_A = pdb_contacts[i[0]]
                    false_ids_A = np.where(intra_res_A[:, 1]-intra_res_A[:, 0] <= 10)[0]
                    intra_res_A = np.delete(intra_res_A, false_ids_A, axis=0)
                    intra_res_interface_A = np.isin(intra_res_A[:, 1], interface_res_A)
                    intra_res_interface_contacts_A = np.sum(intra_res_A[intra_res_interface_A, 2])

                    intra_res_B = pdb_contacts[i[1]]
                    false_ids_B = np.where(intra_res_B[:, 1]-intra_res_B[:, 0] <= 10)[0]
                    intra_res_B = np.delete(intra_res_B, false_ids_B, axis=0)
                    intra_res_interface_B = np.isin(intra_res_B[:, 1], interface_res_B)
                    intra_res_interface_contacts_B = np.sum(intra_res_B[intra_res_interface_B, 2])
                    if interface_contacts_A == 0 and interface_contacts_B == 0:
                        continue
                    if interface_contacts_A == 0 and intra_res_interface_contacts_A == 0:
                        continue
                    if interface_contacts_B == 0 and intra_res_interface_contacts_B == 0:
                        continue

                    chains_c += 1
                    if gene not in genes_c:
                        genes_c.append(gene)
                    if pdb_code not in crystals_c:
                        crystals_c.append(pdb_code)
                    if A[2] == 'TRUE':
                        if gene not in HCc:
                            HCc.append(gene)
                    else:
                        if gene not in LCc:
                            LCc.append(gene)

                    out_.write(gene + '\t' + confidence + '\t' + pdb_code + '\t' + chain_combination + '\t' + str(interface_contacts_A) + '\t' + \
                                str(interface_contacts_B) + '\t' + str(intra_res_interface_contacts_A) + '\t' + str(intra_res_interface_contacts_B) + '\n')



out_.close()


print(len(genes_c), len(HCc), len(LCc), len(crystals_c), chains_c)
print(not_in_pdb, not_two_times_in_structure, no_structure_before_onset, bad_seq_fit)

