# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'intrinsicEditor.ui'
#
# Created: Sat Mar 18 17:11:07 2017
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
        Dialog.resize(389, 486)
        self.gridLayout = QtGui.QGridLayout(Dialog)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(Dialog)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout.addWidget(self.label_2)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.intrinsicTable = QtGui.QTableWidget(Dialog)
        self.intrinsicTable.setObjectName(_fromUtf8("intrinsicTable"))
        self.intrinsicTable.setColumnCount(3)
        self.intrinsicTable.setRowCount(0)
        item = QtGui.QTableWidgetItem()
        self.intrinsicTable.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        self.intrinsicTable.setHorizontalHeaderItem(1, item)
        item = QtGui.QTableWidgetItem()
        self.intrinsicTable.setHorizontalHeaderItem(2, item)
        self.verticalLayout_2.addWidget(self.intrinsicTable)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.intrinsicDone = QtGui.QPushButton(Dialog)
        self.intrinsicDone.setObjectName(_fromUtf8("intrinsicDone"))
        self.horizontalLayout.addWidget(self.intrinsicDone)
        self.intrinsicCancel = QtGui.QPushButton(Dialog)
        self.intrinsicCancel.setObjectName(_fromUtf8("intrinsicCancel"))
        self.horizontalLayout.addWidget(self.intrinsicCancel)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.gridLayout.addLayout(self.verticalLayout_2, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label.setText(_translate("Dialog", "Enter gas-phase basicity or acidity values (kcal/mol) for ", None))
        self.label_2.setText(_translate("Dialog", "residue side chains. Otherwise leave fields blank", None))
        item = self.intrinsicTable.horizontalHeaderItem(0)
        item.setText(_translate("Dialog", "Residue", None))
        item = self.intrinsicTable.horizontalHeaderItem(1)
        item.setText(_translate("Dialog", "GBint", None))
        item = self.intrinsicTable.horizontalHeaderItem(2)
        item.setText(_translate("Dialog", "GAint", None))
        self.intrinsicDone.setText(_translate("Dialog", "Done", None))
        self.intrinsicCancel.setText(_translate("Dialog", "Cancel", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

