import os
import sys
import math
import pprint
import random
import logging
import itertools
import operator
import numpy as np
from datetime import datetime
from predictPrPlus.predictor.CSP_io import *
cimport cython
cimport numpy as np
from libc.stdlib cimport rand, RAND_MAX

# ctypes definitions
DTYPE_double = np.double
DTYPE_intp = np.intp

ctypedef np.double_t DTYPE_double_t
ctypedef np.intp_t DTYPE_intp_t

class peptide(object):
	def __init__(self, siteid, resid, resname, Type, x, y, z, GBI):
		values = []
		for i in x,y,z,GBI:
			if i is not None:
				values.append(float(i))
			else:
				values.append(None)
		self.siteid = siteid
		self.resid = resid
		self.resname = resname
		self.Type = Type
		self.x = values[0]
		self.y = values[1]
		self.z = values[2]
		self.GBI = values[3]
		return

	def __repr__(self):
		return '%s, %s, %s, %s, %s, %s, %s, %s' %(self.siteid, self.resid, self.resname, self.Type, self.x, self.y, self.z, self.GBI)

def build_peptide(options):
	'''
	Construct a list of peptide objects for a FASTA type input sequence

	Notes:
		1) Fixed charges are defined by an 'X' placed at the required position in the sequence
		2) Fixed charges are offset along the y-axis by 5 A - can make this a user variable

	charge_state_range = [3]#range(2,8) # NB: charge state is TOTAL charge state - ie charge state = fixed charges plus protons
					 this is defined in the randomise module -
					 -- number of protons assigned = charge state - fised charges

	'''

	cdef DTYPE_intp_t index, basic_sites

	sequence = options.sequence

	res_spacing = options.resSpacing
	siteid, ionisable_sites = 0, 0
	res_data, ionisable = resdata(options)
	peplist, fixedlist = [], []


	# parse sequence string and create list of ionisation nodes
	for index, element in enumerate(sequence.replace('X','')):
		x_coord = res_spacing * (index+1)

		'''
		When coding negative mode
			- need to remove the N-term entry that is added by default for all peptides
			- add a new C-term entry at the end by default
		'''
		if index == 0:
			peplist.append(
							peptide(siteid,
									index+1,
									'N-term',
									'amide',
									0,
									0,
									None,
									res_data['amide']['intrinsic']
									)
							)
			siteid += 1

		if element in ionisable:
			peplist.append(
							peptide(siteid,
									index+1,
									element,
									'amide',
									x_coord,
									0,
									None,
									res_data['amide']['intrinsic']
									)
							)
			siteid += 1
			peplist.append(
							peptide(siteid,
									index+1,
									element,
									'SC',
									x_coord - res_spacing/2,
									0,
									None,
									res_data[element]['intrinsic']
									)
							)
			siteid += 1
			ionisable_sites += 1

		else:
			peplist.append(
							peptide(siteid,
									index+1,
									element,
									'amide',
									x_coord,
									0,
									None,
									res_data['amide']['intrinsic']
									)
							)
			siteid += 1


	# add entries for fixed charge sites
	for index, element in enumerate(sequence):
		x_coord = res_spacing * (index) + 1.9
		if element == 'X':
			fixedlist.append(
							peptide(siteid,
									None,
									'X',
									'fixed',
									x_coord,
									5,
									None,
									0)
							)
			siteid += 1

	# create state array
	state = create_state_array(fixedlist, peplist)

	return peplist, fixedlist, ionisable_sites, state

def build_peptide_3D(options):
	'''
	Construct a peptide object from a PDB input file

	Notes:

		1) PDB atom type definitions are used to identify the important atoms in a structures
			---> atom types must match PDB format

		2) Important atoms are:
			- Positive ion:
				- backbone carbonyl O for all residues
				- N-terminal amine
				- ARG = CZ
				- HIS = ND1
				- LYS = NZ
				- GLN = OE1
				- TRP = NE1
				- PRO = N

			- Negative ion:
				- TODO

		3) Fixed charges are not directly extracted from PDB FILTERS
			- this would be difficult to achieve accurately in the general case given that metal ions present may/may not be neutralised by coordinating ligands
			- TODO ---> add functionality allowing the user to specify additional charges and coordinates thereof


	-----> rewrite this function
	'''

	PDB  = []
	res_data, ionisable = resdata(options)
	atoms, numResidues, sequence = parse_pdb(options)

	for atom in atoms:

		SLresname = translate(atom.resName)

		# handle special cases of N- and C- termini
		if atom.resIndex == 1 and options.ionType == 'Pos':
			if atom.atomName == 'N': PDB.append([atom , 'N-term'])

		if atom.resIndex == numResidues and options.ionType == 'Neg':
			if atom_type == 'OXT': PDB.append([atom, 'C-term'])

		# parse file body
		if SLresname in ionisable:
		 	# add nodes for SC and amide
			if res_data[SLresname]['atom'] == atom.atomName:
				PDB.append([atom, 'SC'])
			elif res_data['amide']['atom'] == atom.atomName:
				PDB.append([atom, 'amide'])

		# add amide node for non-SC ionisable residues
		else:
			if res_data['amide']['atom'] == atom.atomName: PDB.append([atom, 'amide'])


	# for data in atomLines:

	# 	resname = data[3]

	# 	SLresname = translate(resname)
	# 	atom_type = data[2]

	# 	# handle special cases of N- and C- termini
	# 	if data[5] == '1' and options.ionType == 'Pos':
	# 		if atom_type == 'N': PDB.append(data + ['N-term'])
	# 		# Should the C-Term COOH be ionisable in +'ve mode???	TODO

	# 	if data[5] == str(numResidues) and options.ionType == 'Neg':
	# 		if atom_type == 'OXT': PDB.append(data + ['C-term'])

	# 	# parse body of file
	# 	if SLresname in ionisable:
	# 	 	# add nodes for SC and amide
	# 		if res_data[SLresname]['atom'] == atom_type:
	# 			PDB.append(data + ['SC'])
	# 		elif res_data['amide']['atom'] == atom_type:
	# 			PDB.append(data + ['amide'])

	# 	# add amide node for non-SC ionisable residues
	# 	else:
	# 		if res_data['amide']['atom'] == atom_type: PDB.append(data + ['amide'])


	siteid, ionisable_sites = 0, 0
	peplist, fixedlist = [], []

	for i, element in enumerate(PDB):
		atom, Type = element

		atom.resName = translate(atom.resName)
		# siteid = i
		# resid = int(element[5])
		# resName = translate(element[3])
		# x = element[6]
		# y = element[7]
		# z = element[8]

		if atom.resName in ionisable and Type == 'SC':
			ionisable_sites += 1
		#x_coord = res_specing * count2 # NB residue coordinates specified in angstroms here then converted to metres in charge_repulsion function

		peplist.append(
			peptide(
					siteid,	# siteID
					atom.resIndex, #resid, # resID
					atom.resName, #resName, # resName
					Type, # type
					atom.x, # x
					atom.y, # y
					atom.z, # z
					#GI, # GBI
					res_data[atom.resName]['intrinsic'] if Type == 'SC' else res_data[Type]['intrinsic'] # GBI
					)
		)
		siteid += 1

	# add in fixed sites if needed
	for i, element in enumerate(options.pdbFixed):
		x,y,z = element

		fixedlist.append(
						peptide(siteid,
								None,
								'X',
								'fixed',
								x,
								y,
								z,
								0)
							)
		siteid += 1
		sequence += 'X'

	amides = [x for x in peplist if x.Type == 'amide']
	SCs = [x for x in peplist if x.Type == 'SC']

	# sequence = sequence.split(' ')
	# sequence = [translate(x.strip()) for x in sequence if x != '']
	# sequence = ''.join(sequence)

	options.sequence = sequence
	state = create_state_array(fixedlist, peplist)
	print 'shape of state array is'
	print state.shape
	return peplist, fixedlist, ionisable_sites, state

def create_state_array(fixedlist, peplist):
	'''
	Create state array from list of peptide objects

	State array specs
	Row #
	0 	siteIDs - integer
	1 	fixed - binary
	2	charged - binary
	3 	amide - binary
	4 	SC - binary
	5 	amideN - binary
	6 	SCN - binary
	7 	len parameters for amideN and SCN
		--> [7,0] = peptide length (not including fixed sites)
		--> [7,1] = fixedLength
		--> [7,2] = amidesN length
		--> [7,3] = SCsN length
		--> [7,4] = number of mobile charges: i.e. not fixed

	Rules:

	Fixed, Amide and SC are mutually exclusive
	AmidesN, SCN and charged are mutually exclusive
	charged and neutral are mutually exclusive
	'''

	cdef DTYPE_intp_t fixedLength, pepLength, arrayLength, i
	cdef int[:,::1] state

	fixedLength = len(fixedlist)
	pepLength = len(peplist)

	# determine the combined length of the peptide and fixed lists - this becomes the number of columns in state array
	arrayLength = len(peplist) + len(fixedlist)

	# create state array
	state = np.zeros( 8 * arrayLength, dtype = np.int32).reshape(8,arrayLength)

	# add peptide and fixed siteIDs to state array - not necessary but might be useful for debugging
	for i in range(arrayLength):
		for pep in peplist:
			if pep.siteid == i:
				state[0,i] = i

		for f in fixedlist:
			if f.siteid ==i:
				state[0,i] = i

	# add other parameters to state array
	for i in range(arrayLength):
		for pep in peplist:
			if pep.siteid == i:
				if pep.Type == 'amide':
					state[3,i] = 1
				elif pep.Type == 'SC':
					state[4,i] = 1

		for f in fixedlist:
			if f.siteid ==i:
				state[1,i] = 1

	state[7,0] = pepLength
	state[7,1] = fixedLength

	return state

def get_energies(peptide, fixed, options):
	'''
	Determine energies of all possible pairwise interactions and store as np array

	Inputs:
		1) list of peptide objects
		2) list of fixed charge peptide objects
		3) options parameters

	Output:
		- n*n array of doubles where n is the number of ionisable sites (including fixed charges)
		- array indices directly correspond to the siteid values of peptide objects

	Notes: This function is only called once for each peptide simulated
	'''
	cdef DTYPE_double_t q, vp, rp, constant, e_energy, k
	cdef DTYPE_intp_t key1, key2, j, array_len
	cdef double[:,::1] energies
	cdef double[:] GBIs

	# constants
	q = float(1.60 * math.pow(10,-19))
	vp = float(8.85 * math.pow(10,-12))
	rp = float(options.RVP)
	constant = 4*math.pi*vp*rp

	# create list of all possible charge sites
	charges = peptide + fixed

	print len(peptide), len(fixed)
	# get all possible combinations of 2 sites
	combinations = list(itertools.combinations(charges,2))

	# create array of charge sites
	array_len = len(charges)
	energies = np.zeros((array_len, array_len), dtype = np.double)

	# determine energies of each pairwise interaction and add to array
	for i in combinations:
		key1, x1, y1, z1 = i[0].siteid, i[0].x, i[0].y, i[0].z
		key2, x2, y2, z2 = i[1].siteid, i[1].x, i[1].y, i[1].z

		Dx = x2 - x1
		Dy = y2 - y1

		if options.PDB:
			Dz = z2 - z1
			distance = math.sqrt(math.pow(Dx,2) + math.pow(Dy,2) + math.pow(Dz,2)) * math.pow(10,-10)
		else:
			distance = math.sqrt(math.pow(Dx,2) + math.pow(Dy,2)) * math.pow(10,-10)

			# 10E-10 to convert angstroms to metres

		# formulation of this equation depends on ion type
		e_energy = np.divide(float(q**2), constant*float(math.fabs(distance))) * 6.022 * math.pow(10,20)

		# need to account for attractive force between negative charges and fixed positive charges if present
		if options.ionType == 'Neg':
			if i[0].Type == 'fixed' and i[1].Type == 'fixed':
				pass
			elif i[0].Type == 'fixed' or i[1].Type == 'fixed':
				e_energy = e_energy * -1

		energies[key1, key2] = e_energy
		energies[key2, key1] = e_energy



	# create array for storing GBI values
	GBIs = np.zeros(len(peptide), dtype = np.double)
	for j, res in enumerate(peptide):
		k = res.GBI
		if options.ionType == 'Pos':
			GBIs[j] = k
		elif options.ionType == 'Neg':
			'''
			Note for negative mode
			------------------------------------------
			Need to minimise sum(CR0) + sum(GA)
				-for positive mode, minimise sum(CR0) - sum(GG)

			This is calculated in the charge_repulsion function but don't want to have a ion type check there as the function is called many times

			Negating GA values will produce the same effect
				sum(CR) - sum(-GA)

			Remember to account for this when calculating GAapp
			'''
			GBIs[j] = k * -1

	return energies, GBIs

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double charge_repulsion(double[:,::1] energies,
	double[:] intrinsic,
	int[:] all_charges_array,
	int[:,::1] combs,
	int[:,::1] test_state):

	'''
	Determine the energy of a specified charge configuration

	Inputs:
		1) pairwise electrostatic energy array
		2) site GBI array
		3) *empty* array for storing charge arrays
		4) combinations array defining pairwise interactions of elements in all_charges_array
		5) test charge configuration array

	Outputs:
		1) double representing energy of test charge configuration
	'''

	# type declarations
	cdef double GB, CR, free_energy, gb, cr
	cdef DTYPE_intp_t index, i, row_val, col_val, site, site1, site2

	# init CR and GBI values to 0
	GB, CR, index = 0,0,0

	# populate all_charges_array with charged sites from state array
	for i in range(test_state.shape[1]):
		if test_state[2,i] == 1 or test_state[1,i] == 1:
			all_charges_array[index] = i
			index += 1

	# determine total intrinsic GB
	for i in range(all_charges_array.shape[0]):
		site = all_charges_array[i]
		gb = intrinsic[site]
		GB += gb

	# determine total electrostatic repulsion
	for i in range(combs.shape[1]):

		row_val = combs[0,i]
		col_val = combs[1,i]

		site1 = all_charges_array[row_val]
		site2 = all_charges_array[col_val]

		cr = energies[site1, site2]
		CR += cr

	free_energy = CR - GB * 4.184
	return free_energy

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:,::1] randomise(int charge_state,
	int[:,::1] state):
	'''
	Randomly assign initial charge sites

	Inputs:
		1) charge state
		2) peptide state array

	Outputs:
		1) State array with charge assignments
	'''

	cdef DTYPE_intp_t i, x, num, index

	# neutral = peptide and charge = [] are effectively resetting the peptide charges to 0
	# NB: random proton assignment can be to SCs - not just amides - does this need to be changed?

	# reset all charge values to 0
	for i in range(state.shape[1]):
		state[2,i] = 0

	x = 0
	# get number of mobile charges to be added - i.e. non-fixed charges
	num = len(range(charge_state - state[7,1]))
	state[7,4] = num

	charged_indices = []

	while (x < num):
		index = rand() % state[7,0]
		if index not in charged_indices:
			charged_indices.append(index)
			state[2,index] = 1
			x += 1

	return state

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:,::1] walk(int[:,::1] state,
	double pr_amide,
	double[:,::1] energies,
	double[:] GBIs,
	int[:] all_charges_array,
	int[:,::1] combs):

	'''
	TODO - need to check that the random number generators do not select fixed charge sites for protonation
	'''

	'''
	Move one charge to a test position

	Inputs:
		1) peptide state array
		2) probability of moving charge to an amide site
		3) pairwise electrostatic energies array
		4) GBI array
		5) *empty* array for storing charge arrays
		6) charge site combinations array

	Outpust:
		1) updated peptide state array
	'''
	cdef int[:,::1] new_state
	cdef DTYPE_intp_t index_of_charge_to_replace, charge_site_to_append_to_neutral, count, j, prob1, prob2, len_amide, len_SCsN, index_of_neutral_site_to_protonate, E_min

	new_state = state

	index_of_charge_to_replace = rand() % new_state[7,4]

	charge_site_to_append_to_neutral = 0

	# find the nth protonated site
	count = 0
	for j in range(state.shape[1]):
		if new_state[2,j] == 1:
			if count == index_of_charge_to_replace:
	#			print 'neutralising site: %s' %j

				new_state[2,j] = 0
				charge_site_to_append_to_neutral = j
			count += 1

	prob1 = rand() % 100
	prob2 = rand() % 100

	len_amide = new_state[7,2]
	len_SCsN = new_state[7,3]

	# select site to protonate
	if prob1 >= pr_amide: # >= pr_amide = probability that proton moved to sc i.e. pr_sc
		if len_SCsN > 0:
			if len_SCsN > 0:
				if prob2 >= 5:
					# protonate SC

					index_of_neutral_site_to_protonate = rand() % len_SCsN

					count = 0
					for j in range(state.shape[1]):
						if new_state[6,j] == 1:
							if count == index_of_neutral_site_to_protonate:
								new_state[2,j] = 1
								break
							count += 1

					return new_state

				else:
					prob2 = 1

			if prob2 <5:
				# protonate SC

				index_of_neutral_site_to_protonate = rand() % len_SCsN
				count = 0
				for j in range(state.shape[1]):
					if new_state[6,j] == 1:
						if count == index_of_neutral_site_to_protonate:
							new_state[2,j] = 1
							break
						count += 1

				return new_state

		else:
			# protonate amide - catch all option -> don't think this is actuall ever called

			index_of_neutral_site_to_protonate = rand() % len_amide
			count = 0
			for j in range(state.shape[1]):
				if new_state[5,j] == 1:
					if count == index_of_neutral_site_to_protonate:
						new_state[2,j] = 1
						break
					count += 1
			return new_state

	else:
		if len_amide > 0:
			# protonate minimum amide
			E_min = minimum_amide(energies, all_charges_array, combs, new_state, index_of_charge_to_replace)
			new_state[2,E_min] = 1
			return new_state

		else:
			print('skipping walk function')
			pass

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int minimum_amide(double[:,::1] energies,
	int[:] all_charges_array,
	int[:,::1] combs,
	int[:,::1] new_state,
	DTYPE_intp_t index_of_charge_to_replace):
	'''
	Find the amide site that will result in the lowest energy configuration

	Inputs:
		1) pairwise electrostatic energies array
		2) *empty* array for storing charge arrays
		3) charge site combinations array
		4) peptide state array
		5) index of element in peptide state array to be replaced

	Outputs:
		1) intiger site ID of minimum energy amide site
	'''

	cdef DTYPE_intp_t last_all_charges_index, i, index, min_amide_site, j, site1, site2, row_val, col_val
	cdef DTYPE_double_t min_energy, CR, cr

	last_all_charges_index = all_charges_array.shape[0] - 1
	for i in range(all_charges_array.shape[0]):
		all_charges_array[i] = 0

	# calculate all charges list for fixed sites
	index = 0
	for i in range(new_state.shape[1]):
		if new_state[2,i] == 1 or new_state[1,i] == 1:
			all_charges_array[index] = i
			index += 1

	# iterate through amidesN and calculate pairwise interactions
	# initialise to a big number
	min_energy = 9999999999999999999
	min_amide_site = 0

	for i in range(new_state.shape[1]):
		# test if site i a neutral amide
		if new_state[5,i] == 1:
			#amide = new_state[5,i]
			all_charges_array[last_all_charges_index] = i
			CR = 0

			for j in range(combs.shape[1]):

				row_val = combs[0,j]
				col_val = combs[1,j]

				site1 = all_charges_array[row_val]
				site2 = all_charges_array[col_val]

				cr = energies[site1, site2]
				CR += cr
			# since all sites of interest are amides (that have the same GBI value) subtraction of GBI can be neglected and energy is defined only by sum(CRs)
			if CR < min_energy:
				min_energy = CR
				min_amide_site = i

	return min_amide_site

def probabilities(charge_state, basic_sites):
	''' Define probabilities of proton movements to either SC or backbone as a function of different charge states '''

	if charge_state >= basic_sites:
		pr_sc = 50
		pr_amide = 50
	elif charge_state <= float(basic_sites)/2:
		pr_sc = 95
		pr_amide = 5
	else:
		# Equation giving probability that proton transferred to backbone amide
		(x1,y1) = (float(basic_sites)/2, 5)
		(x2,y2) = (float(basic_sites), 50)
		m = float((y2-y1) / (x2-x1))
		c = 50 - m*x2
		pr_amide = float(charge_state*m + c)
		pr_sc = 100 - pr_amide

	return pr_sc, pr_amide

def energy_exception(final,initial):
	# defines probabilities for accepting changes of higher energy
	dE = (float(final)-float(initial))
	prob = rand() % 100
	if prob >= 10:
		N = 16 * 4.184
	else:
		N = 50 * 4.184
	probability = 1 - dE/N
	if probability < 0:
		probability = 0
	return probability*1000

class siteResults(object):
	def __init__(self, siteid, GBapp, GBI, CR):
		# NB: CR in this instance represents the total electrostatic repulsion on site [siteid] rather that the total colomb repulsion of all sites combined
		# GBapp, GBI and CR represent the average of all configs within cutoff that have this site protonated
		self.siteid = siteid
		self.GBapp = GBapp
		self.GBI = GBI
		self.CR = CR
	def __repr__(self):
		return (4* '%s, ' %(self.siteid, self.GBapp, self.GBI, self.CR))

def write_results(min_config, low_energy_configs, minimum, charge_state, charge_state_range, energies, intrinsic, peptide, fixedlist, charge_freqs, of1, of2, options, q):

	import logging

	'''
	Calculate most reactive site from final configuration and write results to output files

	Inputs:
		1) min_config - single charge config list containing siteids of ionised residues
		2) charge_state - integer - charge state optimised
		3) minimum - double - total energy of minimum configuration
		4) pairwise electrostatic energies array
		5) peptide - list of peptide objects defining fixed charge sites
		6) fixedlist - list of peptide objects defining fixed charge sites
		7) charge_state_range - list of integers representing charge states to be optimised
	'''

	name = options.outFile
	cutoff = options.cutoff
	residue, charge, histogram, histogram_relative = [], [], [], []

	print('Optimisation of +%s charge state complete' %str(charge_state))
	print('Minimum configuration is: %s' %(min_config))
	print('Minimum energy is : %s kJ/mol' %minimum)

	# remove elements from low energy configs list that are not within energy cutoff

	#print(len(low_energy_configs))

	low_energy_configs = [x for x in low_energy_configs if calc_energy(x, fixedlist, energies, intrinsic) < (minimum+(cutoff*4.184))]
	num_low_energy_structures = len(low_energy_configs)

	#print(len(low_energy_configs))

	# write all low energy confs to log file
	for element in low_energy_configs:
		logging.info([sorted(element)])

	# GBapp = apparent GB of most reactive site
	# proton_intrinsic = GBI of this site
	# proton_colomb = total colomb repulsion experienced by this site
	GBapp, proton_intrinsic, proton_colomb = get_results(low_energy_configs, charge_state, energies, intrinsic, fixedlist, options)

	# add GBapp result to MP queue if necessary
	if q is not None:
		q.put(['GBapp',charge_state, GBapp])


	# get total energies for this state
	colomb, GBint, free_energies = thermochem(low_energy_configs, energies, intrinsic, fixedlist)

	# if this is the first set of result to be written, write file headers etc...
	if charge_state == min(charge_state_range):
		of3 = open('%s_peptide.dat' %(name), 'wt', 0) # prints peptide objects - record of params for each site

		# GBapp file headers
		of1.write('#Sequence: %s\n' %options.sequence)
		of1.write('#Ion Polarity: %s\n' %(options.ionType))
		of1.write('#Charge state range: %s-%s (including %s fixed charge(s))\n' %(min(charge_state_range), max(charge_state_range)+len(fixedlist), len(fixedlist)))

		if options.ionType == 'Pos':
			of1.write('#Charge,    GBapp,   Charge(Colomb), Charge(GBint),   Sum(Colomb), Sum(GBint), Structures\n')
		else:
			of1.write('#Charge,    GAapp,   Charge(Colomb), Charge(GAint),   Sum(Colomb), Sum(GAint), Structures\n')

		# protonation frequency file headers
		of2.write('#Sequence: %s\n' %options.sequence)
		of2.write('#Charge state range: %s-%s (including %s fixed charge(s))\n' %(min(charge_state_range), max(charge_state_range)+len(fixedlist), len(fixedlist)))
		of2.write("#Charge, Residue, Resname, Type, Frequency, Relative frequency\n")

		# write data to ptpdie file
		of3.write('siteid, resid, resname, Type, x, y, z, proton, %s\n' %('GBI' if options.ionType == 'Pos' else 'GAI'))

		nodes = peptide + fixedlist

		for p in nodes:
			of3.write('%s, %s, %s, %s, %s, %s, %s, %s\n' %(p.siteid, p.resid, p.resname, p.Type, p.x, p.y, p.z, p.GBI))
		of3.close()

	# write data to GBapp file
	print  >> of1, '{0:^8} {1:^12} {2:^14} {3:^14} {4:^13} {5:^12} {6:^9}'.format(int(charge_state), '{0:.2f}'.format(GBapp), '{0:.2f}'.format(proton_colomb), '{0:.2f}'.format(proton_intrinsic), '{0:.2f}'.format(colomb), '{0:.2f}'.format(GBint), num_low_energy_structures)

	# determine maximum frequency value in charge_freqs dictionary

	# generate histogram data
	hist = frequency(low_energy_configs, peptide)

	max_val_key = max(hist.iteritems(), key=operator.itemgetter(1))[0]
	max_val = hist[max_val_key]

	# write frequency data to file
	print >> of2, '{0:^8} {1:^8} {2:^8} {3:^7} {4:^7} {5:^20}'.format(int(charge_state), -1, 0, 0, 0, 0)

	for i, pep in enumerate(peptide):
		print >> of2, '{0:^8} {1:^8} {2:^8} {3:^7} {4:^7} {5:^20}'.format(charge_state, pep.resid, pep.resname, pep.Type, '{0:.3f}'.format(hist[pep.siteid]), '{0:.3f}'.format(float(hist[pep.siteid])/float(max_val)))

	print >> of2, '{0:^8} {1:^8} {2:^8} {3:^7} {4:^7} {5:^20}'.format(int(charge_state), len(options.sequence)+1, 0, 0, 0, 0)

	return

def get_results(low_energy_configs, charge_state, energies, intrinsic, fixedlist, options):

	'''
	Determine the most unstable charged site - i.e most acidic site for + ions

	Inputs:
		1) min_config - list containing siteids of charged residues
		2) charge state being optimised
		3) pairwise electrostatic energies array
		4) list of peptide objects defining fixed sites

	Outputs:
		1) GBapp value of most unstable site
		2) GBI of this site
		3) total colomb repulsion experiencecd by this site
	'''

	GBapp_dictionary, results = {}, [],

	# create a list of all sites that are protonated in at least one config
	#charge_list = list(itertools.chain.from_iterable(low_energy_configs))
	charge_list = []
	for config in low_energy_configs:
		for site in config:
			if site not in charge_list:
				charge_list.append(site)

	# iterate through sites
	for site in sorted(charge_list):

		# get low energy configs that have this site protonated
		site_charged_configs = []
		for config in low_energy_configs:
			if site in config:
				site_charged_configs.append(config)

		# create holding lists for GBapp, GBI and CR data
		GBapp, GBIs, CRs = [],[],[]

		# calculate GBapp, GBI and CR for each of these configs
		for config in site_charged_configs:
			siteGBapp, totalGBI, totalCR = get_GBapp(config, site, energies, intrinsic, fixedlist, options)
			GBapp.append(siteGBapp)
			GBIs.append(totalGBI)
			CRs.append(totalCR)

		# calculate average values for the configs with this site charged
		av_GBapp = sum(GBapp)/len(GBapp)
		av_GBI = sum(GBIs)/len(GBIs)
		av_CRs = sum(CRs)/len(CRs)

		# add GBapp data to site dictionary
		GBapp_dictionary[site] = av_GBapp

		# append results object to results lists
		data = siteResults(site, av_GBapp, av_GBI, av_CRs)
		results.append(data)

	# get config with minimum GBapp value
	results.sort(key=lambda x: x.GBapp)

	# for i, x in enumerate(results):
	# 	if i < 10:
	# 		print(x)


	min_site = results[0]

	return min_site.GBapp, min_site.GBI, min_site.CR

def get_GBapp(min_config, target, energies, intrinsic, fixedlist, options):
	'''
	Calculate energy (GBapp) of a given charge site

	Inputs:
		1) min_config - list contining siteids of all charged sites - NB, includes fixed sites
		2) target - integer - siteid of charged site of interest
		3) pairwise electrostatic energies array
		4) list of peptide objects defining fixed sites

	Outputs:
		1) apparent GB of target site
		2) intrinsic GB of target site
		3) total CR experienced by target site
	'''

	all_charges = min_config + [x.siteid for x in fixedlist]
	combinations = [[x,target] for x in all_charges if x != target]
	CRs = []

	for element in combinations:
		CRs.append(energies[element[0], element[1]])

	GBI = intrinsic[target] * 4.184
	CR = sum(CRs)

	if options.ionType == 'Pos':
		GBapp = GBI - CR
	elif options.ionType == 'Neg':
		'''
		* -1 to create positive GA values
		These were negated in get energies function to save time -> see note there
		'''
		GBapp = GBI * -1 + CR

	return GBapp, GBI, CR

def thermochem(low_energy_configs, energies, intrinsic, fixedlist):
	'''
	Calculate global thermochem parameters
	'''
	colomb, GBintrinsic, free_energies = [],[],[]

	# calculate total CR and GBI for each charge state
	for config in low_energy_configs:
		charges = config + [x.siteid for x in fixedlist]
		combinations = list(itertools.combinations(charges, 2))
		GBIc = [intrinsic[x] for x in config]
		CRs = []

		for element in combinations:
			CRs.append(energies[element[0], element[1]])

		GBI = sum(GBIc) * 4.184
		CR = sum(CRs)
		free_energy = CR - GBI

		# append data to top level lists
		colomb.append(CR)
		GBintrinsic.append(GBI)
		free_energies.append(free_energy)

	return sum(colomb)/len(colomb), sum(GBintrinsic)/len(GBintrinsic), free_energies

def frequency(low_energy_configs, peptide):
	'''
	Create protonation frequency dictionary from min_energy confs lists
	'''
	hist = {}

	# init all values to 0
	for element in peptide:
		hist[element.siteid] = 0

	for index, element in enumerate(low_energy_configs):
		for x in element:
			hist[x] = hist[x] + 1

	return hist

def get_labels(peptide):
	return [element.resname for element in peptide]

def calc_energy(config, fixedlist, energies, intrinsic):
	'''
	Convenience function - calculate energy directly from a python list of charge site ids
	'''
	all_charges = config + [x.siteid for x in fixedlist]
	combinations = list(itertools.combinations(all_charges, 2))

	CRs = []
	for element in combinations:
		CRs.append(energies[element[0], element[1]])

	GBIs = [intrinsic[x] for x in config]

	return sum(CRs) - (sum(GBIs) * 4.184)

def data_chunks(charges, size):
	'''
	Divides charge states between cores
		- Implemented in MPI runs

	Inputs:
		1) list of integers representing charge states to be optimised
		2) number of cores to be used in the simulation

	Outputs:
		1) nested lists containing charge states to be optimised by each core
	'''

	size = size - 1 # deduct rank 0 that will be master node
	scatterLength = int(len(charges)/size)
	data = [[] for x in range(size)]
	assign = 0
	for item in charges:
		data[assign].append(item)
		assign += 1
		if (assign+1) > len(range(size)):
			assign = 0
	return data

@cython.boundscheck(False)
@cython.wraparound(False)
cdef intersect(int[:,::1] state):

	'''
	Identify the intersection of amides or SCs and neutrals
	Implemented in cythonised version - np.intersect1D adds too much overhead
	'''


	cdef DTYPE_intp_t amideNlen, SCsNlen, i

	amideNlen = 0
	SCsNlen = 0

	for i in range(state.shape[1]):

		# get amide intersection
		# check that site i is amide and is not charged
		if state[3,i] == 1 and state[2,i] == 1:
			# reset amideN entry to 0 if site is charged
			state[5,i] = 0
		elif state[3,i] == 1 and state[2,i] == 0:
			state[5,i] = 1
			amideNlen += 1

		# get SC intersection
		if state[4,i] == 1 and state[2,i] == 1:
			# reset SCN entry to 0 if site is charged
			state[6,i] = 0
		elif state[4,i] == 1 and state[2,i] == 0:
			state[6,i] = 1
			SCsNlen += 1

	state[7,2] = amideNlen
	state[7,3] = SCsNlen

	#return

def get_num_charges(state):
	num_charges = 0
	for i in range(state.shape[1]):
		if state[2,i] == 1 or state[1,i] == 1:
			num_charges += 1
	return num_charges

def get_charged_sites(test_state):
	'''
	Get siteids for charged sites from state array
		- NB: this does not include fixed charges
	'''
	config = []
	for k in range(test_state.shape[1]):
		if test_state[2,k] == 1:
			config.append(k)
	return sorted(config)

def main(options, q = None):

	'''
	State array specs

	Row #

	0 	siteIDs - integer
	1 	fixed - binary
	2	charged - binary
	3 	amide - binary
	4 	SC - binary
	5 	amideN - binary
	6 	SCN - binary
	7 	len parameters for amideN and SCN
		--> [7,0] = peptide length (not including fixed sites)
		--> [7,1] = fixedLength
		--> [7,2] = amidesN length
		--> [7,3] = SCsN length
		--> [7,4] = number of mobile charges: i.e. not fixed

	Rules:

	Fixed, Amide and SC are mutually exclusive
	AmidesN, SCN and charged are mutually exclusive
	charged and neutral are mutually exclusive
	'''

	import time
	import logging
	startTime = datetime.now()

	# run sanity check and cleanup input
	options = sanity_check(options)

	try:
		from mpi4py import MPI
		if MPI.COMM_WORLD.Get_rank() == 0:
			print '\nMPI import successful'
		if MPI.COMM_WORLD.size == 1:
			MPI = None
	except:
		MPI = None

	cdef double[:,::1] energies
	cdef double[:] GBIs
	cdef int[:,::1] combs, prev_state, state, test_state
	cdef int[:] all_charges_array
	cdef DTYPE_intp_t index_array_length, v, iv
	cdef double minimum, current, test
	cdef int end, loop, cycle, min_count, low_E_config_counter, accept, reject, iteration_counter, k, prob3, exception_probability

	optCycle = options.optSteps
	mcCycle = options.mcSteps
	cutoff = options.cutoff
	name = options.outFile
	CSR = options.CSR

	# create output files
	of1 = open('%s_GBapp.dat' %(name),'wt', 0) # stores calculated GB values for each proton
	of2 = open('%s_configs.dat' %(name),'wt', 0) # per residue protonation frequency data

	# prepare charge state list
	charge_state_range = range(CSR[0], CSR[1], CSR[2])

	# create state array for input peptide
	if options.PDB:
		peptide, fixedlist, basic_sites, state = build_peptide_3D(options)
		#print('Need to adapt 3d peptide build function for use with numpy arrays ---> exiting')
		#sys.exit()

	else:
		peptide, fixedlist, basic_sites, state = build_peptide(options)

	# calculate pairwise interaction energies
	energies, GBIs = get_energies(peptide, fixedlist, options)

	if MPI is not None:
		comm = MPI.COMM_WORLD
		rank = comm.Get_rank()
		size = comm.size
		if rank == 0:
			data = None
			chunk = data_chunks(charge_state_range, size)
			for rank in range(1,size):
				comm.send(chunk[rank-1], dest = rank)
		else:
			data = comm.recv(source = 0)
	else:
		rank = 0
		size = 1
		data = charge_state_range
		print('mpi exception')
		global count
		count = 0

	if MPI is not None:
		comm = MPI.COMM_WORLD
		rank = comm.Get_rank()
		size = comm.size

	if options.logfile:
		logging.basicConfig(filename = options.logfile + '-' + str(rank),
							filemode = 'w',
							level = logging.DEBUG,
							format='%(levelname)s %(asctime)s %(message)s',
							datefmt='%m/%d/%Y %H:%M:%S')

		logging.info('Program started')
	logging.info('command line: {}'.format(' '.join(sys.argv)))

	if rank == 0:
		print("\nSequence: %s\n" %str(options.sequence))
		print("Peptide length is: %s" %(len(str(options.sequence).replace('X',''))))
		print('The number of possible protonation sites is: %s' %(state.shape[1]))
		print('The number of fixed charge sites is: %s' %(state[7,1]))
		print('Charge state range (including %s fixed charge(s)) is: %s\n' %(state[7,1], charge_state_range))
		print ('%s processors will be used\n' % size)

		logging.info("Sequence: %s" %str(options.sequence))
		logging.info("Peptide length is: %s" %(len(str(options.sequence).replace('X',''))))
		logging.info('The number of possible protonation sites is: %s' %(state.shape[0]))
		logging.info('The number of fixed charge sites is: %s' %(state[7,1]))
		logging.info('Charge state range (including %s fixed charge(s)) is: %s\n' %(state[7,1], charge_state_range))
		logging.info('%s processors will be used\n' % size)

	if data is not None:
		for charge_state in data:
			print('Rank %s is processing charge state +%s' %(rank, charge_state))
			logging.info('Processing charge state +%s' %(charge_state))
			if charge_state >= len(peptide):
				print('+%s charge state is > #possible protonation sites - passing' %(charge_state))
				pass

			pr_sc, pr_amide = probabilities(charge_state, basic_sites)

			# determine inital proton and neutal arrays
			randomise(charge_state, state)

			# get intersection of neutral amide and SC arrays
			intersect(state)

			# get total number of charges
			num_charges = get_num_charges(state)

			# calculate charge list index combinations for charge repulsion function
			pairwise = list(itertools.combinations(range(num_charges), 2))
			entries = []
			for i in pairwise:
				entries.append([i[0],i[1]])

			combinations_t = zip(*entries)
			index_array_length = int(len(entries))
			rows = np.zeros(index_array_length, dtype = np.int32)
			cols = np.zeros(index_array_length, dtype = np.int32)
			combs = np.zeros([2,index_array_length], dtype = np.int32)

			for iv in range(index_array_length):
				v = combinations_t[0][iv]
				combs[0,iv] = v

			for iv in range(index_array_length):
				v = combinations_t[1][iv]
				combs[1,iv] = v


			# set up arrays for holding GBI and CR data - do this once then overwrite with each cycle
			# --> np array creation methods seem to be somewhat slow
			all_charges_array = np.zeros(charge_state, dtype = np.int32)

			# calculate energy of starting structure
			current = charge_repulsion(energies, GBIs, all_charges_array, combs, state)

			# track energy of minimum energy config identified
			minimum = current

			# holding list for siteids of minimum energy structure
			min_config = []

			# holding list for siteid lists of all minimum energy structures identified - nested lists
			low_energy_configs = []

			# init prev_state to a 2D array
			prev_state = np.zeros([state.shape[0], state.shape[1]], dtype=np.int32)

			if q is not None:
				''' append charge state to MP queue '''
				q.put(['charge_state',charge_state])

			end, loop, cycle, min_count = 0, 0, 0, 0
			iteration_counter = 1
			accept = 0
			reject = 0
			low_E_config_counter = 0

			charge_freqs = {}
			for element in range(len(peptide)):
				charge_freqs[element] = 0

			# Start optimisation loop
			while True:

				# No new mimina found >> begin Monte Carlo randomisation
				if loop > optCycle:
					logging.info('Rank %s, charge %s, end loop number = %s' %(rank, int(charge_state), end))
					randomise(charge_state, state)
					current = charge_repulsion(energies, GBIs, all_charges_array, combs, state)
					intersect(state)
					loop = 0
					end += 1
					if q is not None:
						q.put(['mc', end])

				# mx number of randomisations reached - write results
				if end > mcCycle:
					if q is not None:
						q.put(['done'])

					if rank != 0:
						comm.send([min_config, low_energy_configs, minimum, charge_state, charge_state_range, energies, GBIs, peptide, fixedlist, charge_freqs, of1, of2, options], dest = 0)
						print('sending data from rank %s to 0' %rank)
					else:
						#print(min_config)
						write_results(min_config, low_energy_configs, minimum, charge_state, charge_state_range, energies, GBIs, peptide, fixedlist, charge_freqs, of1, of2, options, q)
						count += 1
					break

				# copy state information to holding array
				prev_state[...] = state

				# move a charge along the peptide
				test_state = walk(state, pr_amide, energies, GBIs, all_charges_array, combs)

				# get energy of test config
				test = charge_repulsion(energies, GBIs, all_charges_array, combs, test_state)

				difference = test-current

				# test if config is within cutoff energy of minimum
				if math.fabs(minimum-test) < (cutoff*4.184): # NB save structures w/in 10 kcal/mol
					config = get_charged_sites(test_state)
					if config not in low_energy_configs:
						logging.info('New Minima: Rank %s, charge %s, energy %s, sites %s' %(rank, charge_state, current, sorted(config)))
						low_energy_configs.append(config)
					min_count += 1

				if test <= current:
					current = test

					if current <= minimum:
						if current < minimum:

							# get charged sites
							config = get_charged_sites(test_state)

							# update min_config with new charges
							min_config = config

							# reset loop counter
							loop = 0
							#time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

							# write data to log file and append to low_energy_confs list
							if config not in low_energy_configs:
								logging.info('New Minima: Rank %s, charge %s, energy %s, sites %s' %(rank, charge_state, current, sorted(config)))
								low_energy_configs.append(config)

							min_count += 1

							# append data to multiprocess queue if running from GUI mode
							if q is not None:
								# create reporting class
								q.put(['opt',iteration_counter, current])
								q.put(['config',low_E_config_counter, config])
								low_E_config_counter += 1

						minimum = current

					state[...] = test_state
					intersect(state)
					accept += 1


				else:
					prob3 = rand() % 1000
					exception_probability = energy_exception(test,current)

					# accept higher energy move if exception condition met
					if prob3 < exception_probability:
						current = test
						state[...] = test_state
						intersect(state)
						accept += 1

					# reject move and revert to previous state
					else:
						loop += 1
						state[...] = prev_state
						reject += 1

						# increment charge site frequency dictionary
						for k in range(test_state.shape[1]):
							if test_state[2,k] == 1:
								charge_freqs[k] = charge_freqs[k] + 1

				iteration_counter += 1

			print('Accept: %s, Reject: %s, Total: %s\n\n' %(accept, reject, iteration_counter))
			#print('Accept pc: %s, Reject pc: %s' %(accept/iteration_counter*100, reject/iteration_counter*100))

	if MPI is not None:
		if rank == 0:
			count2 = 0
			output = []
			charges = [x for x in charge_state_range]
			while int(count2) < int(len(charges)):
					output_data = comm.recv(source = MPI.ANY_SOURCE)
					output.append(output_data)
					for i in sorted(charges):
						index = []
						if count2 == len(charges):
							break
						charge = charges[count2]
						for element in output:
							if element[3] == charge:
								min_config, low_energy_configs, minimum, charge_state, charge_state_range, energies, GBIs, peptide, fixedlist, charge_freqs, of1, of2, options = element
								write_results(min_config, low_energy_configs, minimum, charge_state, charge_state_range, energies, GBIs, peptide, fixedlist, charge_freqs, of1, of2, options, q)
								count2 += 1
								index.append(output.index(element))
						if len(index) > 0:
							for item in index:
								output.pop(item)
					time.sleep(0.5)

	print('Execution time is: ')
	print datetime.now() - startTime
	of1.close()
	of2.close()
	return
