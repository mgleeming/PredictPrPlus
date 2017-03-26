import sys
import argparse
from predictPrPlus.predictor import CSP
# default parameter values
defaultOPT = 1000
defaultMC = 500
defaultRVP = 2
defaultCutoff = 3
defaultRS = 3.8
defaultWS = 10
defaultName = 'outdata'
defaultIonType = 'Pos'

parser = argparse.ArgumentParser(description = 'Optimise protein charge distributions')

parser.add_argument('--CSR',
			help = 'range of charge states to be optimised. Enter lower_limit upper limit',
			nargs = '+',
			required = True)
parser.add_argument('--outFile',
			help = 'name of output file',
			default = defaultName,
			type = str)
parser.add_argument('--optSteps',
			help = 'number of optimisation steps for each MonteCarlo cycle',
			default = defaultOPT,
			type = int)
parser.add_argument('--mcSteps',
			help = 'number of MonteCarlo cycles',
			default = defaultMC,
			type = int)
parser.add_argument('--RVP',
			help = 'relative vacuum permeability value used in charge repulsion calculation',
			default = defaultRVP,
			type = float)
parser.add_argument('--cutoff',
			help = 'energy cutoff for saving minimum energy conformations',
			default = defaultCutoff,
			type = float)
parser.add_argument('--resSpacing',
			help = 'distance between backbone carbonyl groups of adjacent residues',
			default = defaultRS,
			type = float)
parser.add_argument('--walkStep',
			help = 'maximum standard proton movement in optimisation cycles',
			default = defaultWS,
			type = int)
parser.add_argument('--logfile',
			help = 'file name for logfile',
			type = str)
parser.add_argument('--PDB',
			help = 'Local PDB filename for 3d optimisation',
			type = str)
parser.add_argument('--sequence',
			help = 'peptide sequence to optimise',
			type = str)
parser.add_argument('--ionType',
			help = 'Ion polarity to model. Either "Pos" or "Neg"',
			default = defaultIonType,
			type = str)


def gui_init(q, args):
	options = gui_options(args)
	CSP.main(options, q)
	q.put(['Done'])
	return

class gui_options(object):
	def __init__(self, args):
		for k, v in args.iteritems():
			setattr(self, k, v)

#	def __repr__(self):
#		return ( 11 * '%s, ' %(self.CSR, self.outFile, self.optSteps, self.mcSteps, self.RVP, self.cutoff, self.resSpacing, self.walkStep, self.logfile, self.sequence, self.PDB))

# def run():
# 	#options = gui_options([2,5,1], 'outdata', 1000, 25, 2, 3, 3.8, 10, 'log', sequence = 'GRGRGRGRG')
# 	options = gui_options([2,30,1], 'outdata', 1000, 500, 2, 3, 3.8, 10, None, sequence = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG')
# 	main(options)

if __name__ == '__main__':
	options = parser.parse_args()
	CSP.main(options)
