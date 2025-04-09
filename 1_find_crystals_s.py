import pickle
import os
import numpy
import edlib
import resource
import time
import gc
import sys

org = 'humans'

with open("/path/to/pdb_seq_all.pick", "rb") as handle:
	in_ = pickle.load(handle)

out = {}
for p, (i, j) in enumerate(in_.items()):
	hash_tag = abs(hash(j)) % (10 ** 10)
	out.setdefault(hash_tag, [[], ''])
	out[hash_tag][0].append(i)
	out[hash_tag][1] = j

with open("/path/to/pdb_seq_hashed.pick", 'wb') as handle:
	pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)


with open("/path/to/pdb_seq_hashed.pick", "rb") as h:
    seq_crystals = pickle.load(h)
with open("/path/to/aa_seq_"+org+".pick", "rb") as h:
    seq_genes = pickle.load(h)

bin = 19

out = {}
start_time = time.time()
for p, (gene, gene_seq) in enumerate(seq_genes.items()):
	if p%20 != bin:
		continue
	gene_length = len(gene_seq)
	if gene_length < 25:
		continue
	assemblies = {}
	for pp, (hash, cryst_data) in enumerate(seq_crystals.items()):
		cryst_seq = cryst_data[1]
		crystal_length = len(cryst_seq)
		dd = edlib.align(cryst_seq, gene_seq, mode="MW", task="locations")
		test_score = 1-(dd['editDistance']/gene_length)
		if test_score < 0.3:
			continue
		ddd = edlib.align(cryst_seq, gene_seq, mode="HW", task="locations")
		score = 1-(ddd['editDistance']/gene_length)
		align_start = ddd['locations'][0][0]
		align_stop = ddd['locations'][0][1]
		align_length = align_stop-align_start
		coverage = align_length/gene_length
		for A in cryst_data[0]:
			B = A.split('-')
			if len(B) == 3:
				crystal, assembly, chain = B
			else:
				crystal, assembly, chain_a, chain_b = B
				chain = chain_a+'-'+chain_b
			assembly_key = crystal+'-'+assembly
			assemblies.setdefault(assembly_key, {})
			assemblies[assembly_key][chain] = [score, align_start, align_stop, coverage]
	if len(assemblies) > 0:
		out[gene] = assemblies
		print(p, gene, len(assemblies))
	if len(out)%100 == 0:
		with open("/path/to/"+org+"_"+str(bin)+".pick", 'wb') as handle:
			pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open("/path/to/"+org+"_"+str(bin)+".pick", 'wb') as handle:
	pickle.dump(out, handle, protocol=pickle.HIGHEST_PROTOCOL)
