from ConvertData import convertData, binnedConvertData
from driver import main00
from Box_single_pion_studies import CalculateBoxAffinity

import ROOT as root
import numpy as np
from array import array
'''
USER SET FLAGS AND PATHS
'''

Binned = False #set true to use pre-binned kinematcs; false to use tree_MC and bin after with Box.py
xlsxFileName = "xlsx/May23_test_500k"
inRootFileName = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/May_23/for_box_500k.root"
outRootFileName = "root_files/May23_test.root"

plotFileName = "driver_test_500k.pdf"

'''
STEP 1: Create xlsx file with kinematics from MC
'''
if(Binned):
    binnedConvertData(xlsxFileName + ".xlsx", inRootFileName)
else:
    convertData(xlsxFileName + ".xlsx", inRootFileName)
    
'''
STEP 2: Calculate ratios using driver, save them to a root file
'''

tabs,ratios = main00(xlsxFileName)
tab = tabs[0]

print(f"creating root file with kinematics needed to run Box.py")

file = root.TFile.Open(outRootFileName, "RECREATE")

tree = root.TTree("tree_driver","tree_driver")

R0_t = array('d',[0])
R1_t = array('d',[0])
R2_t = array('d',[0])
Q2_t = array('d',[0])
z_t = array('d',[0])
x_t = array('d',[0])
qTQ_hadron_t = array('d',[0])
pT_t = array('d',[0])
tmdaff_t = array('d',[0])

tree.Branch('R0', R0_t,'R0_t/D')
tree.Branch('R1', R1_t,'R1_t/D')
tree.Branch('R2', R2_t,'R2_t/D')
tree.Branch('x', x_t,'x_t/D')
tree.Branch('Q2', Q2_t,'Q2_t/D')
tree.Branch('pT', pT_t,'pT_t/D')
tree.Branch('z', z_t,'z_t/D')
tree.Branch('qTQ_hadron', qTQ_hadron_t,'qTQ_hadron_t/D')
tree.Branch('tmdaff', tmdaff_t,'tmdaff_t/D')
for i in range(tab.shape[0]): #iterate over each event
    R0_t[0] = tab['R0'][i]
    R1_t[0] = tab['R1'][i]
    R2_t[0] = tab['R2'][i]
    Q2_t[0] = tab['Q2'][i]
    pT_t[0] = tab['pT'][i]
    z_t[0] = tab['z'][i]
    x_t[0] = tab['x'][i]
    qTQ_hadron_t[0] = tab['qTQ'][i]
    tmdaff_t[0] = tab['tmdaff'][i]
    tree.Fill()

tree.Write()
file.Close()

CalculateBoxAffinity("Analysis/BoxAffinity" + outRootFileName, True, plotFileName)