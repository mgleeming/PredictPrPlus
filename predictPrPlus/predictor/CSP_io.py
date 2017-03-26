import os, sys
import itertools
import pprint
import operator
from prody import parsePDB


class Atom (object):
	def __init__(self):
		return

def alias_resname(resname):
	alias = ['HSP', 'HSE', 'HSD']
	if resname in alias:
		return 'HIS'
	else:
		return resname

def parse_pdb(options):
	'''
	Retrieve 'ATOM' entries from PDB file and determine number of residues
	'''
	# parse atom lines
	atoms, sequence = [], ''

	alias = ['HSP', 'HSE', 'HSD']

	pdb = parsePDB(options.PDB)

	numResidues = len(pdb.getResnames())
	numatoms = pdb.numAtoms()

	# maintain separate index counter in case atom.isprotein is False
	index = 0
	for i in xrange(numatoms):
		atom = pdb[i]
		if atom.isprotein:

			# alias atom name if necessary
			if atom.getResname() in alias:
				atom.setResname('HIS')

			a = Atom()
			a.atomIndex = index
			a.atomName = atom.getName()
			a.resName = atom.getResname()
			a.resIndex = atom.getResindex() + 1 # prody res indices start at 0
			a.x, a.y, a.z = atom.getCoords()

			if a.atomName == 'CA':
				sequence += translate(a.resName)

			atoms.append(a)
			index += 1

	return atoms, numResidues, sequence

def resdata(options):
	'''
	Compile dict of relevant GB/GA values
	'''
	res_data = {}
	if options.ionType == 'Pos':
		if options.gbDict:
			for k,v in options.gbDict.iteritems():
				if v['GBI'] is not None:
					res_data[k] = {'intrinsic': v['GBI'], 'atom': v['atom']}
		else:
			# should create a defaults file and read this in at statup
			res_data['R'] = {'atom': 'CZ', 'intrinsic': 251.3}
			res_data['H'] = {'atom': 'ND1', 'intrinsic': 244.8}
			res_data['K'] = {'atom': 'NZ', 'intrinsic': 241.0}
			res_data['Q'] = {'atom': 'OE1', 'intrinsic': 237.4}
			res_data['W'] = {'atom': 'NE1', 'intrinsic': 234.3}
			res_data['P'] = {'atom': 'N', 'intrinsic': 234.3}
			res_data['amide'] = {'atom': 'O', 'intrinsic': 221.6}
			res_data['N-term'] = {'atom': 'N', 'intrinsic': 221.6}


	elif options.ionType == 'Neg':
		if options.gbDict:
			for k,v in options.gbDict.iteritems():
				if v['GAI'] is not None:
					res_data[k] = {'intrinsic': v['GAI'], 'atom': v['atom']}
		else:
			res_data['D'] = {'atom': 'OD1', 'intrinsic': 326.5}
			res_data['E'] = {'atom': 'OE1', 'intrinsic': 328.7}
			res_data['Y'] = {'atom': 'OH', 'intrinsic': 336.4}
			res_data['amide'] = {'atom': 'N', 'intrinsic': 337.6}
			res_data['C-term'] = {'atom': 'OXT', 'intrinsic': 328.7}

			# update dict with GBdict entries if any
			keys = res_data.keys()
	else:
		raise Exception ('Ion type specification error')

	ionisable = []
	for key in res_data.keys():
			ionisable.append(key)
	return res_data, ionisable

def translate(resname):

	rescodes = {}
	rescodes['ALA'] = 'A'
	rescodes['ARG'] = 'R'
	rescodes['ASN'] = 'N'
	rescodes['ASP'] = 'D'
	rescodes['CYS'] = 'C'
	rescodes['GLN'] = 'Q'
	rescodes['GLU'] = 'E'
	rescodes['GLY'] = 'G'
	rescodes['HIS'] = 'H'
	rescodes['ILE'] = 'I'
	rescodes['LEU'] = 'L'
	rescodes['LYS'] = 'K'
	rescodes['MET'] = 'M'
	rescodes['PHE'] = 'F'
	rescodes['PRO'] = 'P'
	rescodes['SER'] = 'S'
	rescodes['THR'] = 'T'
	rescodes['TRP'] = 'W'
	rescodes['TYR'] = 'Y'
	rescodes['VAL'] = 'V'
	return rescodes[resname]


def sanity_check(options):
	'''
	Sanity check all inputs and prepare input data for simulation

	Tasks:
		1) check output files defined do not already exist
		2) check for a valid input sequence
		3) check CSR input
			- check min > 1 and max < number of ionisable sites
		4) check RVP > 1
	'''

	# strip whitespace from sequence and check validity
	if options.sequence:
		options.sequence = ''.join([i for i in options.sequence.strip().upper() if i.isalpha()])

		# check items in sequence are valid
		residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X']

		for i, res in enumerate(options.sequence):
			if res not in residues:
				raise Exception ('Invalid sequence\n\tResidue %s not found' %i)

	# Check CSR
	chargeError = 'Error in charge state input'

	assert len(options.CSR) == 3, chargeError

	try:
		options.CSR = [int(x) for x in options.CSR]
	except:
		raise Exception (chargeError)

	assert options.CSR[0] < options.CSR[1] and (options.CSR[1] - options.CSR[0]) >= options.CSR[2], chargeError

	# Check RVP
	RVPError = 'Error in RVP input'
	try:
		options.RVP = float(options.RVP)
	except:
		raise Exception (RVPError)

	assert options.RVP >= 1, RVPError

	# Check resSpacing
	resSpacingError = 'Error in resSpacing input'
	try:
		options.resSpacing = float(options.resSpacing)
	except:
		raise Exception (RVPError)

	assert options.resSpacing > 0, resSpacingError

	# Check cutoff input
	cutoffError = 'Error in energy cutoff input'
	try:
		options.cutoff = float(options.cutoff)
	except:
		raise Exception (cutoffError)

	# Check mcSteps and optSteps inputs
	optStepsError = 'Error in optSteps input'
	mcStepsError = 'Error in mcSteps input'
	try:
		options.mcSteps = int(options.mcSteps)
	except:
		raise Exception (optStepsError)
	try:
		options.optSteps = int(options.optSteps)
	except:
		raise Exception (mcStepsError)

	# Check if output files already exists
	'''
	Remove for development
	options.outFile = options.outFile.split('.')[0]
	if os.path.isfile('%s_GBapp.dat' %options.outFile): raise Exception ('Output file "%s_GBapp.dat" already exists' %options.outFile)
	if os.path.isfile('%s_configs.dat' %options.outFile): raise Exception ('Output file "%s_configs.dat" already exists' %options.outFile)
	if os.path.isfile('%s_peptide.dat' %options.outFile): raise Exception ('Output file "%s_peptide.dat" already exists' %options.outFile)
	'''
	return options
