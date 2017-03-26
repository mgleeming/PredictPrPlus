import sip
sip.setapi('QString', 2)

from PyQt4 import QtCore, QtGui
import multiprocessing
from multiprocessing import Process, Queue

import os, sys, math, random, time
import pyqtgraph as pg

from gui import CSP_gui
from gui import intrinsicEditor as intE
from gui import pdbFixedEditor as pdbFE

from predictor import run_CSP_sim

class guiprocs (QtGui.QDialog, CSP_gui.Ui_Dialog):
	'''
	guiprocs inherits from CSP_GUI.UI class
	'''
	def __init__(self, parent = None):

		''' set plot figure options '''
		pg.setConfigOption('background','w')

		''' init GUI '''
		super(guiprocs, self).__init__(parent)
		self.setupUi(self)
		self.setWindowTitle('predictPrPlus (v1.0.0)')

		''' make connections '''
		self.form_connections()

		''' multiprocess job holding list '''
		self.jobs = []

		''' disable entry fields until an input option is selected '''
		self.PDB_file_line.setEnabled(False)
		self.selectPDBButton.setEnabled(False)
		self.seq_entry.setEnabled(False)
		self.addFixedButton.setEnabled(False)

		''' set ion type default to pos '''
		self.ionPositive.setChecked(True)

		''' set FASTA as default input '''
		self.inputFASTA.setChecked(True)

		''' init multiprocessing queue'''
		self.q = Queue(maxsize = 0)

		''' set timer '''
		self.timer = QtCore.QTimer()
		self.timer.timeout.connect(self.update_GUI_data)

		''' set default values in parameter fields '''
		sequence = self.seq_entry.setText('MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG')
		RVP = self.RVP.setText('2')
		res_spacing = self.res_spacing.setText('3.8')
		cutoff = self.cutoff.setText('3')
		MC_steps = self.MC_steps.setText('500')
		opt_step = self.opt_step.setText('1000')
		#walk_step = self.walk_step.setText('10')
		OFN = self.OFN.setText('test')
		start = self.CSR_start.setText('2')
		stop = self.CSR_stop.setText('5')
		step = self.CSR_step.setText('1')
		write_log = self.write_log.isChecked()

		''' data storage for plot updates '''
		self.opt_x, self.opt_y = [], []
		self.freq_x, self.freq_y = [], []
		self.charge_states, self.GBapps = [], []
		self.MC_steps_progress_bar = int(self.MC_steps.text())

		''' Create dict for custom GB/GA values '''
		self.GBdict = {}
		self.pdbFixed = []

		''' set plot config options '''
		self.opt_plot.showGrid(x = True, y = True)
		self.opt_plot.showLabel('bottom', show = True)
		#self.opt_plot.showLabel('top', show = True)
		self.opt_plot.showLabel('left', show = True)
		self.opt_plot.showLabel('right', show = True)
		self.opt_plot.setLabel(axis = 'bottom', text = 'Iteration')
		self.opt_plot.setLabel(axis = 'left', text = 'Total Energy (kJ/mol)')
		self.opt_plot.setLabel(axis = 'top', text = '')
		self.opt_plot.setLabel(axis = 'right', text = '')
		self.opt_plot.setLogMode(x = True)

		self.freq_plot.showGrid(x = True, y = True)
		self.freq_plot.showLabel('bottom', show = True)
		#self.freq_plot.showLabel('top', show = True)
		self.freq_plot.showLabel('left', show = True)
		self.freq_plot.showLabel('right', show = True)
		self.freq_plot.setLabel(axis = 'bottom', text = 'Charge State')
		self.freq_plot.setLabel(axis = 'left', text = 'GBapp (kJ/mol)')
		self.freq_plot.setLabel(axis = 'top', text = '')
		self.freq_plot.setLabel(axis = 'right', text = '')

	def main(self):
		self.show()

	def form_connections(self):
		self.run.clicked.connect(self.get_input_data)
		self.inputFASTA.toggled.connect(self.enableInput)
		self.inputPDB.toggled.connect(self.enableInput)
		self.selectPDBButton.clicked.connect(self.selectPDBFile)
		self.editIntrinsicButton.clicked.connect(self.getGBdict)
		self.addFixedButton.clicked.connect(self.fixedPDBs)
		return

	def getGBdict(self):
		dialog = intrinsicEditor(parent = self, residueDict = self.GBdict)
		self.GBdict = dialog.editIntrinsic(self.GBdict)
		return

	def enableInput(self):
		if self.inputFASTA.isChecked() == True:
			self.seq_entry.setEnabled(True)
			self.PDB_file_line.setEnabled(False)
			self.selectPDBButton.setEnabled(False)
			self.addFixedButton.setEnabled(False)
			self.res_spacing.setEnabled(True)

		elif self.inputPDB.isChecked() == True:
			self.PDB_file_line.setEnabled(True)
			self.selectPDBButton.setEnabled(True)
			self.addFixedButton.setEnabled(True)
			self.seq_entry.setEnabled(False)
			self.res_spacing.setEnabled(False)
		return

	def selectPDBFile(self):
		PDBFile = QtGui.QFileDialog.getOpenFileName()
		self.PDB_file_line.setText(str(PDBFile))
		return

	def fixedPDBs(self):
		dialog = pdbFixedEditor(parent = self, fixedData = self.pdbFixed)
		self.pdbFixed = dialog.fixedEditor(self.pdbFixed)
		return

	def run_job(self, function):#, input):
		''' Call this to start a new job '''
		runner.job_function = function
		#runner.job_input = input
		runner_thread.start()

	def get_input_data(self):
		''' Pull input parameters from GUI and init optimisation '''

		# prompt user for output directory
		directory = str(QtGui.QFileDialog.getExistingDirectory(self,"Select output directory")) + '/'

		sequence = self.seq_entry.toPlainText()
		RVP = self.RVP.text()
		res_spacing = self.res_spacing.text()
		cutoff = self.cutoff.text()
		MC_steps = self.MC_steps.text()
		opt_step = self.opt_step.text()
		#walk_step = self.walk_step.text()
		OFN = directory + str(self.OFN.text())
		start = self.CSR_start.text()
		stop = self.CSR_stop.text()
		step = self.CSR_step.text()
		CSR = [str(start),str(stop),str(step)]
		PDB = str(self.PDB_file_line.text())
		ionPolarity = 'Pos' if self.ionPositive.isChecked() == True else 'Neg'

		if self.write_log.isChecked():
			logfile = directory + 'log_' + OFN
		else:
			logfile = False

		''' clear plot data '''
		del self.opt_x[:]
		del self.opt_y[:]
		del self.freq_x[:]
		del self.freq_y[:]

		self.opt_plot.clear()
		self.freq_plot.clear()
		self.clear_plots()

		args = {
			'CSR' : CSR,
			'outFile': OFN,
			'optSteps': int(opt_step),
			'mcSteps': int(MC_steps),
			'RVP' : float(RVP),
			'cutoff' : float(cutoff),
			'resSpacing' : float(res_spacing),
			#'walkStep' : float(walk_step),
			'logfile' : logfile,
			'ionType' : ionPolarity,
			'sequence' : sequence if str(sequence) != '' and self.inputFASTA.isChecked() == True else None,
			'PDB' : PDB if str(PDB) != '' and self.inputPDB.isChecked() == True else None,
			'pdbFixed' : self.pdbFixed,
			'gbDict' : self.GBdict
		}

		''' run simulation '''
		self.p1 = Process(target = run_simulation, args = (self.q, args,))
		self.p1.start()

		self.timer.start(50)
		return

	def update_GUI_data(self):
		'''
		pyqtgraph plot reference
		http://www.pyqtgraph.org/documentation/graphicsItems/plotitem.html
		'''

		opt_data = []
		config_data = []
		updates = []

		# get data from queue
		while not self.q.empty():
			update = self.q.get()
			updates.append(update)

		for update in updates:
			if update[0] == 'done': self.clear_plots()

		for update in updates:
			if update[0] == 'opt': opt_data.append(update)
			elif update[0] == 'config': config_data.append(update)
			elif update[0] == 'charge_state': self.CS_processing_value.setText(str(update[1]))
			elif update[0] == 'mc':
				pc = '%s%%' %int(float(update[1]) / (float(self.MC_steps_progress_bar)+1) * 100)
				self.progress_processing_value.setText(str(pc))
			elif update[0] == 'GBapp':
				self.charge_states.append(int(update[1]))
				self.GBapps.append(float(update[2]))
			elif update[0] == 'Done':
				time.sleep(1)
				self.p1.terminate()


		for i in reversed(opt_data):
			self.opt_x.append(i[1]) # iteration number
			self.opt_y.append(i[2]) # current energy
			self.energy_processing_value.setText(str(i[2]))

		for i in reversed(config_data):
			for site in i[2]:
				self.freq_x.append(int(site))
				self.freq_y.append(int(i[1]))


		self.opt_plot.plot(self.opt_x, self.opt_y, pen = None, symbol = 'o')
		#self.freq_plot.plot(self.freq_x, self.freq_y, pen = None, symbol = 'o')
		self.freq_plot.plot(self.charge_states, self.GBapps, pen = None, symbol = 'o')
		return

	def clear_plots(self):
		del self.opt_x[:]
		del self.opt_y[:]
		del self.freq_x[:]
		del self.freq_y[:]
		self.opt_plot.clear()
		self.freq_plot.clear()

class pdbFixedEditor(QtGui.QDialog, pdbFE.Ui_Dialog):
	'''
	Specification of fixed charge sites for PDB input runs
	'''
	def __init__(self, parent = None, fixedData = None):
		super(pdbFixedEditor, self).__init__(parent)
		self.setupUi(self)
		self.setWindowTitle('Specify fixed charges')

		self.make_dialog_connections()

		# set tablewidget behaviour properties
		self.coordList.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
		self.coordList.resizeRowsToContents()
		self.coordList.resizeColumnsToContents()
		self.coordList.horizontalHeader().setStretchLastSection(True)
		self.fixedData = fixedData

		if self.fixedData is not None:
			for c, i in enumerate(self.fixedData):
				rowcount = self.coordList.rowCount()
				self.coordList.insertRow(rowcount)

				self.coordList.setItem(rowcount, 0, QtGui.QTableWidgetItem(i[0]))
				self.coordList.setItem(rowcount, 1, QtGui.QTableWidgetItem(i[1]))
				self.coordList.setItem(rowcount, 2, QtGui.QTableWidgetItem(i[2]))



	def make_dialog_connections(self):
		self.done.clicked.connect(self.accept)
		self.cancel.clicked.connect(self.reject)
		self.removeCharge.clicked.connect(self.remove_selected_target)
		self.removeAll.clicked.connect(self.remove_all)
		self.addCharge.clicked.connect(self.add_new_node)
		return

	def remove_selected_target(self):
		# find highlighted row and remove entry
		highlighted_row = self.coordList.selectionModel().selectedRows()[0].row()
		self.coordList.removeRow(highlighted_row)
		return

	def remove_all(self):
		numRows = self.coordList.rowCount() + 1
		for i in reversed(range(numRows)):
			self.coordList.removeRow(i)
		return

	def add_new_node(self):
		xc = str(self.xEntry.text())
		yc = str(self.yEntry.text())
		zc = str(self.zEntry.text())

		try:
			x = float(xc)
			y = float(yc)
			z = float(zc)
		except:
			return

		nextRow = self.coordList.rowCount()
		self.coordList.insertRow(nextRow)
		self.coordList.setItem(nextRow, 0, QtGui.QTableWidgetItem(str(xc)))
		self.coordList.setItem(nextRow, 1, QtGui.QTableWidgetItem(str(yc)))
		self.coordList.setItem(nextRow, 2, QtGui.QTableWidgetItem(str(zc)))

		# reset point specification fields
		self.xEntry.setText('')
		self.yEntry.setText('')
		self.zEntry.setText('')
		return

	def returnValues(self):
		rows = self.coordList.rowCount()

		fixedList = []
		for i in xrange(rows):
			x = str(self.coordList.item(i,0).text())
			y = str(self.coordList.item(i,1).text())
			z = str(self.coordList.item(i,2).text())
			fixedList.append([x,y,z])
		return fixedList

	@staticmethod
	def fixedEditor(pdbFixed, parent = None):
		dialog = pdbFixedEditor(parent = parent, fixedData = pdbFixed)
		result = dialog.exec_()
		pdbFixed = dialog.returnValues()
		return pdbFixed

class intrinsicEditor(QtGui.QDialog, intE.Ui_Dialog):
	'''
	Panel allowing user to define custom intrinsic GA/GB values
	'''
	def __init__(self, residueDict = {}, parent = None ):
		super(intrinsicEditor, self).__init__(parent)
		self.setupUi(self)
		self.setWindowTitle('Intrinsic GB/GA values')

		self.make_dialog_connections()

		# set table behaviour options
		self.intrinsicTable.resizeRowsToContents()
		#self.intrinsicTable.resizeColumnsToContents()
		#self.intrinsicTable.horizontalHeader().setStretchLastSection(True)

		self.residueDict = residueDict

		if self.residueDict == {}:
			# define defaults
			self.residueDict = {

			'G': {'GBI': None, 'GAI': None, 'atom': 'CA'},
			'A': {'GBI': None, 'GAI': None, 'atom': 'CB'},
			'S': {'GBI': None, 'GAI': None, 'atom': 'OG'},
			'T': {'GBI': None, 'GAI': None, 'atom': 'OG1'},
			'C': {'GBI': None, 'GAI': None, 'atom': 'SG'},
			'V': {'GBI': None, 'GAI': None, 'atom': 'CB'},
			'L': {'GBI': None, 'GAI': None, 'atom': 'CG'},
			'I': {'GBI': None, 'GAI': None, 'atom': 'CB'},
			'M': {'GBI': None, 'GAI': None, 'atom': 'CG'},
			'P': {'GBI': 234.3, 'GAI': None, 'atom': 'N'},
			'F': {'GBI': None, 'GAI': None, 'atom': 'CB'},
			'Y': {'GBI': None, 'GAI': 336.4, 'atom': 'OH'},
			'W': {'GBI': 234.3, 'GAI': None, 'atom': 'NE1'},
			'D': {'GBI': None, 'GAI': 326.5, 'atom': 'OD1'},
			'E': {'GBI': None, 'GAI': 328.7, 'atom': 'OE1'},
			'N': {'GBI': None, 'GAI': None, 'atom': 'OD1'},
			'Q': {'GBI': 237.4, 'GAI': None, 'atom': 'OE1'},
			'H': {'GBI': 244.8, 'GAI': None, 'atom': 'ND1'},
			'K': {'GBI': 241.0, 'GAI': None, 'atom': 'NZ'},
			'R': {'GBI': 251.3, 'GAI': None, 'atom': 'CZ'},
			'amide': {'GBI': 221.6, 'GAI': 337.6, 'atom': 'O'},
			'N-term': {'GBI': 221.6, 'GAI': None, 'atom': 'N'},
			'C-term': {'GBI': None, 'GAI': 328.7, 'atom': 'OXT'}

			}

		self.dict_keys = sorted(self.residueDict.keys())

		# add rows to list
		self.intrinsicTable.setRowCount(len(self.dict_keys))

		# add defaults to table
		for i, residue in enumerate(self.dict_keys):
			GB = self.residueDict[residue]['GBI']
			GA = self.residueDict[residue]['GAI']

			GB = GB if GB != None else ''
			GA = GA if GA != None else ''

			self.intrinsicTable.setItem(i, 0, QtGui.QTableWidgetItem(str(residue)))
			self.intrinsicTable.setItem(i, 1, QtGui.QTableWidgetItem(str(GB)))
			self.intrinsicTable.setItem(i, 2, QtGui.QTableWidgetItem(str(GA)))

			# center text
			self.intrinsicTable.item(i,0).setTextAlignment(QtCore.Qt.AlignCenter)
			self.intrinsicTable.item(i,1).setTextAlignment(QtCore.Qt.AlignCenter)
			self.intrinsicTable.item(i,2).setTextAlignment(QtCore.Qt.AlignCenter)

		# handle column widths
		self.intrinsicTable.setColumnWidth(1,100)
		self.intrinsicTable.setColumnWidth(2,100)
		self.intrinsicTable.horizontalHeader().setResizeMode(0, QtGui.QHeaderView.Stretch)


	def make_dialog_connections(self):
		self.intrinsicDone.clicked.connect(self.accept)
		self.intrinsicCancel.clicked.connect(self.reject)
		return

	def returnValues(self):
		rows = self.intrinsicTable.rowCount()

		# get entries from table
		for i in xrange(rows):
			res = str(self.intrinsicTable.item(i,0).text())
			GB = str(self.intrinsicTable.item(i,1).text())
			GA = str(self.intrinsicTable.item(i,2).text())

			GB = float(GB) if GB.strip() != '' else None
			GA = float(GA) if GA.strip() != '' else None

			self.residueDict[res]['GBI'] = GB
			self.residueDict[res]['GAI'] = GA

		return self.residueDict

	@staticmethod
	def editIntrinsic(gbdict, parent = None):
		dialog = intrinsicEditor(parent = parent, residueDict = gbdict)
		result = dialog.exec_()
		residueDict = dialog.returnValues()
		return residueDict

''' multiprocess target function '''
def run_simulation(q, args):
	run_CSP_sim.gui_init(q, args)
	return


def main():
	app = QtGui.QApplication(sys.argv)
	gui = guiprocs()
	gui.main()
	sys.exit(app.exec_())
	return


if __name__ == '__main__':
	multiprocessing.freeze_support()
	main()
