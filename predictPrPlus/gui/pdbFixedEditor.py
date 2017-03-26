# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'pdbFixedEditor.ui'
#
# Created: Thu Sep  1 16:33:38 2016
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(690, 395)
        self.gridLayout_2 = QtGui.QGridLayout(Dialog)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_6 = QtGui.QLabel(Dialog)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.verticalLayout.addWidget(self.label_6)
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.coordList = QtGui.QTableWidget(Dialog)
        self.coordList.setObjectName(_fromUtf8("coordList"))
        self.coordList.setColumnCount(3)
        self.coordList.setRowCount(0)
        item = QtGui.QTableWidgetItem()
        self.coordList.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        self.coordList.setHorizontalHeaderItem(1, item)
        item = QtGui.QTableWidgetItem()
        self.coordList.setHorizontalHeaderItem(2, item)
        self.gridLayout.addWidget(self.coordList, 0, 0, 3, 1)
        self.removeCharge = QtGui.QPushButton(Dialog)
        self.removeCharge.setObjectName(_fromUtf8("removeCharge"))
        self.gridLayout.addWidget(self.removeCharge, 0, 1, 1, 1)
        self.removeAll = QtGui.QPushButton(Dialog)
        self.removeAll.setObjectName(_fromUtf8("removeAll"))
        self.gridLayout.addWidget(self.removeAll, 1, 1, 1, 1)
        spacerItem = QtGui.QSpacerItem(20, 108, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 2, 1, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.label = QtGui.QLabel(Dialog)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout.addWidget(self.label)
        self.xEntry = QtGui.QLineEdit(Dialog)
        self.xEntry.setText(_fromUtf8(""))
        self.xEntry.setObjectName(_fromUtf8("xEntry"))
        self.horizontalLayout.addWidget(self.xEntry)
        self.label_4 = QtGui.QLabel(Dialog)
        self.label_4.setText(_fromUtf8(""))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.horizontalLayout.addWidget(self.label_4)
        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout.addWidget(self.label_2)
        self.yEntry = QtGui.QLineEdit(Dialog)
        self.yEntry.setObjectName(_fromUtf8("yEntry"))
        self.horizontalLayout.addWidget(self.yEntry)
        self.label_5 = QtGui.QLabel(Dialog)
        self.label_5.setText(_fromUtf8(""))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.horizontalLayout.addWidget(self.label_5)
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.horizontalLayout.addWidget(self.label_3)
        self.zEntry = QtGui.QLineEdit(Dialog)
        self.zEntry.setObjectName(_fromUtf8("zEntry"))
        self.horizontalLayout.addWidget(self.zEntry)
        self.gridLayout.addLayout(self.horizontalLayout, 3, 0, 1, 1)
        self.addCharge = QtGui.QPushButton(Dialog)
        self.addCharge.setObjectName(_fromUtf8("addCharge"))
        self.gridLayout.addWidget(self.addCharge, 3, 1, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.done = QtGui.QPushButton(Dialog)
        self.done.setObjectName(_fromUtf8("done"))
        self.horizontalLayout_2.addWidget(self.done)
        self.cancel = QtGui.QPushButton(Dialog)
        self.cancel.setObjectName(_fromUtf8("cancel"))
        self.horizontalLayout_2.addWidget(self.cancel)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.gridLayout_2.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label_6.setText(_translate("Dialog", "Enter coordinates of fixed charge sites", None))
        item = self.coordList.horizontalHeaderItem(0)
        item.setText(_translate("Dialog", "X coordinate", None))
        item = self.coordList.horizontalHeaderItem(1)
        item.setText(_translate("Dialog", "Y coordinate", None))
        item = self.coordList.horizontalHeaderItem(2)
        item.setText(_translate("Dialog", "Z coordinate", None))
        self.removeCharge.setText(_translate("Dialog", "Remove Charge", None))
        self.removeAll.setText(_translate("Dialog", "Remove All", None))
        self.label.setText(_translate("Dialog", "X  ", None))
        self.label_2.setText(_translate("Dialog", "Y ", None))
        self.label_3.setText(_translate("Dialog", "Z ", None))
        self.addCharge.setText(_translate("Dialog", "Add Charge", None))
        self.done.setText(_translate("Dialog", "Done", None))
        self.cancel.setText(_translate("Dialog", "Cancel", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

