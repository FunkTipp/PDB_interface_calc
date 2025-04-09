import re
import os
import numpy as np
import pickle
import time
import collections
from scipy.spatial import distance
import matplotlib.pyplot as plt
import itertools
import freesasa
radii_class = freesasa.Classifier.getStandardClassifier('naccess')

def get_contacts(coord_dict_in, dist, gap):
	all_dists = distance.squareform(distance.pdist(coord_dict_in[:,1:], metric='euclidean'))
	all_gaps = coord_dict_in[:,0] - coord_dict_in[:,0][np.newaxis].T
	cond_ids = np.where((all_dists <= dist) & (all_gaps >= gap))
	res_array = np.array([coord_dict_in[cond_ids[0],0], coord_dict_in[cond_ids[1],0]]).T
	set_res_array, count = np.unique(res_array, return_counts=True, axis=0)
	#array row = [resA, resB, strengh]
	return np.concatenate((set_res_array, count[np.newaxis].T), axis=1).astype(np.int_)

def get_surface(coords):
	surA = freesasa.calcCoord(np.ndarray.flatten(coords[:, 1:4]), coords[:, 4].reshape((len(coords[:, 4]), 1)))
	surB = []
	for surface_id, res_id in enumerate(coords[:, 0]):
		surB.append([res_id, surA.atomArea(surface_id)])
	surC = np.array(surB)
	res_ids, res_ids_inv = np.unique(surC[:, 0], return_inverse=True)
	surD = np.zeros((len(res_ids), 2))
	surD[:, 0] = res_ids
	np.add.at(surD[:, 1], res_ids_inv, surC[:, 1])
	return surD

aa_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', \
			'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', \
			'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', \
			'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', \
			'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'SEC': 'U'}

with open("path/to/list/of/pdbs.pick", "rb") as handle:
	assemblies = pickle.load(handle)

path = '/path/to'
print(len(assemblies))
bin = 3
for p, assembly in enumerate(assemblies):

	if p%4 != bin:
		continue
	print(p, assembly)

	coordinate_dict_p = {}
	A = open("/path/tp/pdbs/"+assembly+".cif").read()
	id = A.index("_atom_site.pdbx_PDB_model_num")
	B = A[id:].split('\n')
	old_res_id = '-'
	for C in B[1:]:
		if '#' in C:
			break
		D = C.split()
		res = D[5]
		if res not in aa_code:
			continue
		chain = D[6]
		res_id_pre = D[8]
		if res_id_pre == '.':
			continue
		res_id = int(res_id_pre)
		X = float(D[10])
		Y = float(D[11])
		Z = float(D[12])
		if D[2] == 'S':
			atom = 'O'
		else:
			atom = D[2]
		radius = radii_class.radius(D[-1], atom)
		coordinate_dict_p.setdefault(chain, [])
		coordinate_dict_p[chain].append([res_id, X, Y, Z, radius])
	coordinate_dict = {}
	for i, j in coordinate_dict_p.items():
		coordinate_dict[i] = np.array(j)

	#intra chain features
	intra_contacts_dict = {}
	single_chain_surface_dict = {}
	for chain, coords in coordinate_dict.items():
		#contacts
		coords_in = coords[:, :3]
		intra_contacts_dict[chain] =  get_contacts(coords_in, 5, 6)
		#surface
		single_chain_surface_dict[chain] = get_surface(coords)

	if os.path.isdir(path + assembly[1:3]) == False:
		os.mkdir(path + assembly[1:3])

	with open(path + assembly[1:3]+"/"+assembly+"_contacts_intra.pick", "wb") as handle:
		pickle.dump(intra_contacts_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
	with open(path + assembly[1:3]+"/"+assembly+"_surface_intra.pick", "wb") as handle:
		pickle.dump(single_chain_surface_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

	if len(coordinate_dict) < 2:
		continue

	#inter chain features
	inter_contacts_dict = {}
	paired_chain_surface_dict = {}
	comb = itertools.combinations(list(coordinate_dict.keys()), 2)
	for co in comb:
		arA = coordinate_dict[co[0]]
		arB = coordinate_dict[co[1]]
		#contacts
		cA = arA[:, 1:4]
		cB = arB[:, 1:4]
		dists = distance.cdist(cA, cB)
		interface_stepA = np.where(dists <= 5)
		if len(interface_stepA[0]) == 0:
			continue
		interface_stepB = np.array((arA[interface_stepA[0], 0], arB[interface_stepA[1], 0])).T
		int_res, int_count = np.unique(interface_stepB, return_counts=True, axis=0)
		#out columns
		## 0: id chain A in crystal
		## 1: id chain B in crystal
		## 2: sum contacts <= 5 angstrom
		inter_contacts_dict[co[0]+'_'+co[1]] = np.c_[int_res, int_count].astype(int)

		#surface
		id_arA = arA.copy()
		id_arB = arB.copy()
		id_arA[:, 0] += 10000000
		id_arB[:, 0] += 20000000
		dual_coords = np.r_[id_arA, id_arB]
		paired_chain_surface_dict[co[0]+'_'+co[1]] = get_surface(dual_coords)

	if len(inter_contacts_dict) != 0:
		with open(path + assembly[1:3]+"/"+assembly+"_contacts_inter.pick", "wb") as handle:
			pickle.dump(inter_contacts_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

	inter_chain_surface_dict = {}
	for chain_ids, interfaces in paired_chain_surface_dict.items():
		chainA, chainB = chain_ids.split('_')
		ref_A = single_chain_surface_dict[chainA].copy()
		ref_A[:, 0] += 10000000
		ref_B = single_chain_surface_dict[chainB].copy()
		ref_B[:, 0] += 20000000
		combAB = np.r_[ref_A, ref_B]
		surface_dist_A = combAB[:, 1]-interfaces[:, 1]
		surface_dist_B = np.c_[combAB[:, 0], surface_dist_A]
		chainA_id = np.where((surface_dist_B[:, 0] < 20000000) & (surface_dist_B[:, 1] != 0))[0]
		chainB_id = np.where((surface_dist_B[:, 0] > 19999999) & (surface_dist_B[:, 1] != 0))[0]
		if len(chainA_id) == 0 or len(chainA_id) == 0:
			continue
		outA = surface_dist_B[chainA_id, :]
		outA[:, 0] -= 10000000
		inter_chain_surface_dict[chainA+'_'+chainB] = outA
		outB = surface_dist_B[chainB_id, :]
		outB[:, 0] -= 20000000
		inter_chain_surface_dict[chainB+'_'+chainA] = outB
	if len(inter_chain_surface_dict) != 0:
		with open(path + assembly[1:3]+"/"+assembly+"_surface_inter.pick", "wb") as handle:
			pickle.dump(inter_chain_surface_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

	all_coords = np.empty((0, 5))
	all_chain_ids = {}
	for p, (c, x) in enumerate(coordinate_dict.items()):
		all_chain_ids[c] = (p+1)*10000000
		id_x = x.copy()
		id_x[:, 0] += (p+1)*10000000
		all_coords = np.r_[all_coords, id_x]
	all_chains_surface = get_surface(all_coords)

	all_surfaces = {}
	for chain, interface in single_chain_surface_dict.items():
		all_A = interface.copy()
		all_A[:, 0] += all_chain_ids[chain]
		all_B = np.where((all_chain_ids[chain]-1 < all_chains_surface[:, 0]) & (all_chains_surface[:, 0] < all_chain_ids[chain]+10000000))[0]
		all_C = all_A[:, 1]-all_chains_surface[all_B, 1]
		all_D = np.where(all_C != 0)[0]
		all_surfaces[chain] = np.c_[all_A[all_D, 0]-all_chain_ids[chain], all_C[all_D]]

	with open(path + assembly[1:3]+"/"+assembly+"_surface_all.pick", "wb") as handle:
		pickle.dump(all_surfaces, handle, protocol=pickle.HIGHEST_PROTOCOL)
