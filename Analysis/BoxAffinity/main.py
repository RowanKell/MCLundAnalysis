# EXAMPLE USAGE
# python3 main.py --no-useArgs


from ConvertData import convertData, binnedConvertData
from driver import main00, maxmin
from Box_single_pion_studies import CalculateBoxAffinity

import ROOT as root
import numpy as np
from array import array
import os
def checkdir(path):
    if not os.path.exists(path): 
        os.makedirs(path)
import datetime

x = datetime.datetime.now()
today = x.strftime("%B_%d")
import argparse
parser = argparse.ArgumentParser(description='Affinity Calculation Program')

parser.add_argument('--outRootName', type=str, default="no_root_name",
                        help='Name for root file produced by driver.py')
parser.add_argument('--fileFromLundAnalysis', type=str, default="no_root_name",
                        help='Name for root file produced by LundAnalysis.C that contains kinematics')
parser.add_argument('--plotName', type=str, default="no_plot_name",
                        help='Name for plot created by Box.py')
parser.add_argument('--xlsxFileName', type=str, default="no_xlsx_name",
                        help='Name for excel file created by ConvertData.py')
parser.add_argument('--useDriver', action=argparse.BooleanOptionalAction,
                        help='If True, uses driver.py to calculate ratios instead of LundAnalysis.C')
parser.add_argument('--createMaxMin', action=argparse.BooleanOptionalAction,
                        help='If True, creates maxmin tree for plotting histogram')
parser.add_argument('--multipleFiles', action=argparse.BooleanOptionalAction,
                        help='If True, runs with whole directory')
parser.add_argument('--useArgs', action=argparse.BooleanOptionalAction,
                        help='If True, uses command line arguments to set parameters')


args = parser.parse_args()
useArgs = args.useArgs


'''
USER SET FLAGS AND PATHS
'''
if not useArgs:
    useDriver = True
#     fileFromLundAnalysis = "Files_Spring_24/August_15/"
    fileFromLundAnalysis = "Slurm_Spring_24/August_15/Run_1_single_pion/" #for multipleFiles
    createMaxMin = False
    multipleFiles = True
    if(useDriver):
        Binned = False #set true to use pre-binned kinematcs; false to use tree_MC and bin after with Box.py
        xlsxFileName = "xlsx/July_8_100_driver_original_R2_three_files"
        inRootFileName = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/OutputFiles/" + fileFromLundAnalysis
        outDayDir ="root_files/" + today + "/"
        checkdir(outDayDir)
#         outRootFileName = "/" + outDayDir + "old_R2_driver_MCNP_all_low.root"
        outRootFileName = "/" + outDayDir + "old_R2_driver_low.root"
        plotFileName = "driver_july_8_old_R2_6_files.pdf"
        plot_title = "Driver Affinity old R2 high affinity bin"
        calcAff = False
        highAff = False
#         useMCNP = False
        useMCNP = False
    else:
        inRootFileName = "/OutputFiles/" + fileFromLundAnalysis
        plotFileName = "driver_july_8_all_files_low.pdf"
        plot_title = "driver Affinity all files low TMD aff binning"
else:
    useDriver = args.useDriver
    fileFromLundAnalysis = args.fileFromLundAnalysis
    createMaxMin = args.createMaxMin
    multipleFiles = args.multipleFiles
    if(useDriver):
        Binned = False #set true to use pre-binned kinematcs; false to use tree_MC and bin after with Box.py
        xlsxFileName = args.xlsxFileName
        inRootFileName = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/OutputFiles/" + fileFromLundAnalysis
        outRootFileName = args.outRootName
        plotFileName = args.plotName
    else:
        inRootFileName = "/OutputFiles/" + fileFromLundAnalysis
        plotFileName = args.plotName

'''
DRIVER METHOD
'''
if(useDriver):
    '''
    UPDATE JULY 8: removing xlsx step bc they can only have 1m rows which is too few. Now just use root files
    STEP 1: Create xlsx file with kinematics from MC
    '''
#     if(Binned):
#         binnedConvertData(xlsxFileName + ".xlsx", inRootFileName)
#     else:
#         convertData(xlsxFileName + ".xlsx", inRootFileName,multipleFiles)

    '''
    STEP 2: Calculate ratios using driver, save them to a root file
    '''

    tabs,ratios = main00(inRootFileName,highAff,useMCNP)
    tab = tabs[0]

    print(f"creating root file with kinematics needed to run Box.py")

    file = root.TFile.Open("/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Analysis/BoxAffinity" + outRootFileName, "RECREATE")

    tree = root.TTree("tree_driver","tree_driver")

    R0_t = array('d',[0])
    R1_t = array('d',[0])
    R2_t = array('d',[0])
    Q2_t = array('d',[0])
    z_t = array('d',[0])
    x_t = array('d',[0])
    zeta_t = array('d',[0])
    xi_t = array('d',[0])
    qTQ_hadron_t = array('d',[0])
    pT_BF_t = array('d',[0])
    tmdaff_t = array('d',[0])
    
    M_ki_t  = array('d',[0])
    M_kf_t  = array('d',[0])
    delta_k_T_t  = array('d',[0])
    ki_T_t  = array('d',[0])
    
    theta_deltak_t  = array('d',[0])
    theta_H_t  = array('d',[0])
    theta_ki_t  = array('d',[0])
    
    tree.Branch('R0', R0_t,'R0_t/D')
    tree.Branch('R1', R1_t,'R1_t/D')
    tree.Branch('R2', R2_t,'R2_t/D')
    tree.Branch('x', x_t,'x_t/D')
    tree.Branch('Q2', Q2_t,'Q2_t/D')
    tree.Branch('pT_BF', pT_BF_t,'pT_BF_t/D')
    tree.Branch('z', z_t,'z_t/D')
    tree.Branch('qTQ_hadron', qTQ_hadron_t,'qTQ_hadron_t/D')
    tree.Branch('tmdaff', tmdaff_t,'tmdaff_t/D')
    
    #["M_ki","M_kf","delta_k_T","ki_T"]
    tree.Branch('M_ki', M_ki_t,'M_kf_t/D')
    tree.Branch('M_kf', M_kf_t,'M_kf_t/D')
    tree.Branch('delta_k_T', delta_k_T_t,'delta_k_T_t/D')
    tree.Branch('ki_T', ki_T_t,'ki_T_t/D')
    
    tree.Branch('xi', xi_t,'xi_t/D')
    tree.Branch('zeta', zeta_t,'zeta_t/D')
    
    tree.Branch('theta_ki', theta_ki_t,'theta_ki/D')
    tree.Branch('theta_H', theta_H_t,'theta_H/D')
    tree.Branch('theta_deltak', theta_deltak_t,'theta_deltak/D')
    
    for i in range(tab.shape[0]): #iterate over each event
        R0_t[0] = tab['R0'][i]
        R1_t[0] = tab['R1'][i]
        R2_t[0] = tab['R2'][i]
        Q2_t[0] = tab['Q2'][i]
        pT_BF_t[0] = tab['pT_BF'][i]
        z_t[0] = tab['z'][i]
        x_t[0] = tab['x'][i]
        zeta_t[0] = tab['zeta'][i]
        xi_t[0] = tab['xi'][i]
        qTQ_hadron_t[0] = tab['qTQ'][i]
        tmdaff_t[0] = tab['tmdaff'][i]
        
        M_ki_t[0] = tab['M_ki'][i]
        M_kf_t[0] = tab['M_kf'][i]
        delta_k_T_t[0] = tab['delta_k_T'][i]
        ki_T_t[0] = tab['ki_T'][i]
        
        theta_ki_t[0] = tab['theta_ki'][i]
        theta_H_t[0] = tab['theta_H'][i]
        theta_deltak_t[0] = tab['theta_deltak'][i]
        
        tree.Fill()
        
    tree.Write()
    file.Close()
    if(createMaxMin):#create maxmin tree and write to same root file
        maxmin("/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Analysis/BoxAffinity" + outRootFileName, ratios)

    if(calcAff):
        CalculateBoxAffinity("Analysis/BoxAffinity" + outRootFileName, True, plotFileName,plot_title)
    
'''
BOX Method
'''

if not useDriver:
    CalculateBoxAffinity(inRootFileName, False, plotFileName,multipleFiles,plot_title)