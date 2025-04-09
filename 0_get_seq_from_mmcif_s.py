import re
import os
import numpy as np
import pickle
import time
import collections

#read the resolution index file from the pdb database
out = {}
A = open("/path/to/resolu.idx").read()
B = A.split('\n')
for C in B[1:]:
	D = C.split('\t;\t')
	if len(D) == 2 and len(D[1]) != 0:
		out[D[0]] = float(D[1])
with open("/path/to/resolu.pick", "wb") as handle:
	pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)

#read the entry type file from the pdb database
out = []
A = open("/path/to/pdb_entry_type.txt").read()
B = A.split('\n')
for C in B:
	D = C.split('\t')
	if len(D) != 3:
		continue
	if 'prot' in D[1]:
		out.append(D[0])
with open("/path/to/index.pick", "wb") as handle:
	pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)

aa_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', \
			'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', \
			'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', \
			'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', \
			'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'SEC': 'U'}

with open("/path/to/index.pick", "rb") as handle:
	structures = pickle.load(handle)

bin = 0

out = {}
for p, s in enumerate(structures):
	if p%20 != bin:
		continue
	seq = {}
	A = open("/path/to/"+s[1:-1]+"/"+s+"-assembly1.cif").read()
	id = A.index("_atom_site.pdbx_PDB_model_num")
	B = A[id:].split('\n')
	old_res_id = '-'
	for C in B[1:]:
		if '#' in C:
			break
		D = C.split()
		res = D[5]
		res_id = D[8]
		if res_id == old_res_id:
			continue
		if res_id == '.':
			continue
		chain = D[6]
		if '-' in chain:
			continue
		seq.setdefault(chain, '')
		if res not in aa_code:
			seq[chain] += '-'
		else:
			seq[chain] += aa_code[res]
		old_res_id = res_id

	taken = []
	for i, j in seq.items():
		if j in taken:
			continue
		if len(j) < 25:
			continue
		if j.count('-') > len(j)*0.1:
			continue
		taken.append(j)
		out[s+'_'+i] = j
	if len(out)%100 == 0:
		print(p, len(out))
		with open("/path/to/pdb_seq_"+str(bin)+".pick", 'wb') as handle:
			pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open("/path/to/pdb_seq_"+str(bin)+".pick", 'wb') as handle:
	pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)
