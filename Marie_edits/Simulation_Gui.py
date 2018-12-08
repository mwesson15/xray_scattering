from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from Simulation import Simulation

class FDTDWidget(QMainWindow):
    def __init__(self, parent=None):
        super(FDTDWidget, self).__init__(parent)
        # button1 = QPushButton('Button 1')
        self.fontsize=14
        button1 = QPushButton('Run Simulation')
        button1.setFont(QFont('SansSerif',self.fontsize))
        self.grid = QGridLayout()
        self.params = []

        self.h = self.make_dialog('Height (nm)')
        self.w = self.make_dialog('Width (nm)')
        self.sm = self.make_dropdown('Substrate Material', 'silicon', 'gold')
        self.nm = self.make_dropdown('Nanostructure Material', 'none', 'gold', 'silicon')
        self.ns = self.make_dropdown('Nanostructure Shape', 'none', 'rectangle', 'sphere')
        self.se = self.make_dialog('Source Energy (eV)')
        self.se.setText("7350")

        self.ia = self.make_dialog('Incident Angle (deg)')



        grid = self.grid
        grid.addWidget(button1, len(self.params)+1,2)
        grid.setSpacing(15)
        main_frame = QWidget()
        main_frame.setLayout(grid)
        # self.connect(botton1, SIGNAL("clicked()"),self.on_button)
        button1.clicked.connect(self.on_button)

        self.setWindowTitle('FDTD Simulation')   

        self.setCentralWidget(main_frame)

    def on_button(self, n):
        sim = Simulation()
        # print(self.h.text(), self.w.text(), self.sm.currentText(), self.nm.currentText(), self.se.text(), self.ia.text())
        sim.run_simulation(self.h.text(), self.w.text(), self.sm.currentText(), self.nm.currentText(), self.ns.currentText(),self.se.text(), self.ia.text())

    def make_dialog(self, parameter):
        line = QLineEdit()
        line.setFont(QFont('SansSerif',self.fontsize))
        self.params.append(parameter)
        line_label = QLabel(parameter)
        line_label.setFont(QFont('SansSerif',self.fontsize))
        self.grid.addWidget(line_label, len(self.params), 1, Qt.Alignment(2))
        self.grid.addWidget(line, len(self.params), 2)
        return line
    def make_dropdown(self, parameter, *args):
        line = QComboBox()
        line.setFont(QFont('SansSerif',self.fontsize))
        self.params.append(parameter)
        for i in args:
            line.addItem(i)
        line_label = QLabel(parameter)
        line_label.setFont(QFont('SansSerif',self.fontsize))
        self.grid.addWidget(line_label, len(self.params), 1, Qt.Alignment(2))
        self.grid.addWidget(line, len(self.params), 2)
        return line

    # def button_click(self):
    #     # shost is a QString object
    #     shost = self.le.text()
    #     print shost






if __name__ == "__main__":
    import sys
    # if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)

# if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
#     PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication(sys.argv)
    form = FDTDWidget()
    form.show()
    app.exec_()