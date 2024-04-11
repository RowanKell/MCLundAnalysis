from ROOT import TFile,TTree,TCanvas,TH1F,gStyle,TLatex,gPad,TLegend,TLorentzVector,TH2F,TLine,TF1,TBox,RDataFrame,TPad,TF2
import ROOT
import numpy as np
import awkward as awk
import uproot
import pandas as pd
import matplotlib.pyplot as plot
from pandas import read_excel 
from copy import deepcopy
from ipywidgets import *
import logging, os 
logging.disable(logging.WARNING) 
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import tensorflow as tf
print("tf.__version__", tf.__version__)

# simple version for working with CWD
file_dir = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Slurm_Spring_24/April_8/Run_1/"
num_files = len([name for name in os.listdir(file_dir) if not os.path.isdir(name)])
file_names = [name for name in os.listdir(file_dir) if not os.path.isdir(name)]

# file_dir = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Slurm_Spring_24/Feb29/Run_1/"
# num_files = 1
# file_names = ["file_1.root"]

tree_MC_list = []
tree_x_list = []
tree_z_h_list = []
tree_Mh_list = []
tree_qTdivQ_list = []


for i in range(num_files):
    try:
        tree_MC_list.append(uproot.open(file_dir + file_names[i]+ ":tree_MC"))
        tree_x_list.append(uproot.open(file_dir + file_names[i]+ ":tree_x_bins"))    
        tree_z_h_list.append(uproot.open(file_dir + file_names[i]+ ":tree_z_h_bins"))    
        tree_Mh_list.append(uproot.open(file_dir + file_names[i]+ ":tree_Mh_bins"))
        tree_qTdivQ_list.append(uproot.open(file_dir + file_names[i]+ ":tree_qTQ_bins"))
    except uproot.exceptions.KeyInFileError as e:
        print(f"exception: {e}\nexception for file {file_names[i]}; continuing")
        continue
        
    
#0 is z, 1 is x, 2 is mh
zarray = np.array([np.array([np.zeros(7)] * 6)] * num_files)
xarray = np.array([np.array([np.zeros(7)] * 6)] * num_files)
Mharray = np.array([np.array([np.zeros(7)] * 7)] * num_files)
qTdivQarray = np.array([np.array([np.zeros(9)] * 7)] * num_files)


xkinematics = np.array(["z_h", "Q2", "pT", "R0", "R1_p", "R2"])
zkinematics = np.array(["x", "Q2", "pT", "R0", "R1_p", "R2"])
Mhkinematics = np.array(["x", "z_h", "Q2", "pT", "R0", "R1_p", "R2"])
qTdivQkinematics = np.array(["x", "z_h", "Q2", "pT", "R0", "R1_p", "R2"])

#These arrays each hold an array for each variable, meaning that the first bin of variables is in the first index of every kinematics array

#piplus
#z
for i in range(num_files):
    z_iter = 0
    for var in zkinematics:
        zarray[i][z_iter] = tree_z_h_list[i][var].array(library='np')
        z_iter += 1
    x_iter = 0
    for var in xkinematics:
        xarray[i][x_iter] = tree_x_list[i][var].array(library='np')
        x_iter += 1
    Mh_iter = 0
    for var in Mhkinematics:
        Mharray[i][Mh_iter] = tree_Mh_list[i][var].array(library='np')
        Mh_iter += 1
    qTdivQ_iter = 0
    for var in qTdivQkinematics:
        qTdivQarray[i][qTdivQ_iter] = tree_qTdivQ_list[i][var].array(library='np')
        qTdivQ_iter += 1

zarray_t = np.array([np.array([np.zeros(6)] * 7)] * num_files)
xarray_t = np.array([np.array([np.zeros(6)] * 7)] * num_files)
Mharray_t = np.array([np.array([np.zeros(7)] * 7)] * num_files)
qTdivQarray_t = np.array([np.array([np.zeros(7)] * 9)] * num_files)
for i in range(num_files):
    xarray_t[i] = np.transpose(xarray[i])
    zarray_t[i] = np.transpose(zarray[i])
    Mharray_t[i] = np.transpose(Mharray[i])
    qTdivQarray_t[i] = np.transpose(qTdivQarray[i])
    
collinear_region_name = 'collinear'
current_region_name = 'current'
target_region_name = 'target'
TMD_region_name = 'TMD'
soft_region_name = 'soft'
collinear_lable_name = 'collinearaff'
target_lable_name = 'targetaff'
current_lable_name = 'currentaff'
TMD_lable_name = 'tmdaff'
soft_lable_name = 'softaff'

tmd_model_name = '../../SIDIS-Affinity/models/final_%s' % TMD_region_name
tmd_model = tf.keras.models.load_model(tmd_model_name)
target_model_name = '../../SIDIS-Affinity/models/final_%s' % target_region_name
target_model = tf.keras.models.load_model(target_model_name)
collinear_model_name = '../../SIDIS-Affinity/models/final_%s' % collinear_region_name
collinear_model = tf.keras.models.load_model(collinear_model_name)
current_model_name = '../../SIDIS-Affinity/models/final_%s' % current_region_name
current_model = tf.keras.models.load_model(current_model_name)
soft_model_name = '../../SIDIS-Affinity/models/final_%s' % soft_region_name
soft_model = tf.keras.models.load_model(soft_model_name)

Mhbins = np.linspace(0.3,1.3,7)
xbins = np.array([0.1,0.13,0.16,0.19,0.235,0.3,0.5])
zbins = np.array([0.35,0.43,0.49,0.55,0.62,0.7,0.83])
qTdivQbins = np.array([0.1,0.3,0.5,0.8,1.5,2,2.5,3,4])
# qTdivQbins = np.linspace(0.1,0.7,7)

def calculator(array, region, binType, binnedVariable = 0):
    R0max = 0.3
    R1max = 0.3
    R2max = 0.3
    if binType == "x":
        z = array[0]
        Q2 = array[1]
        pT = array[2]
        x = binnedVariable
    elif binType == "z":
        x = array[0]
        Q2 = array[1]
        pT = array[2]
        z = binnedVariable
    elif (binType == "Mh") or (binType == "qTdivQ"):
        x = array[0]
        z = array[1]
        Q2 = array[2]
        pT = array[3]
        
    test_features = pd.DataFrame({'pT':pT,'Q2':Q2,'x':x,'z':z,'R0max':R0max,'R1max':R1max,'R2max':R2max},index=[0])

    if region == 'tmd':
        prediction = tmd_model.predict(test_features).flatten()
        
    elif region == 'target':
        prediction = target_model.predict(test_features).flatten()
        
    elif region == 'collinear':
        prediction = collinear_model.predict(test_features).flatten()

    elif region == 'soft':
        prediction = soft_model.predict(test_features).flatten()

    else:
        prediction = current_model.predict(test_features).flatten()

    return prediction[0] #returns affinity value

colxaffinity = np.zeros(7)
colzaffinity = np.zeros(7)
colMhaffinity = np.zeros(7)
colqTdivQaffinity = np.zeros(9)
TMDxaffinity = np.zeros(7)
TMDzaffinity = np.zeros(7)
TMDMhaffinity = np.zeros(7)
TMDqTdivQaffinity = np.zeros(9)
Currentxaffinity = np.zeros(7)
Currentzaffinity = np.zeros(7)
CurrentMhaffinity = np.zeros(7)
CurrentqTdivQaffinity = np.zeros(9)

region = "collinear"
region2 = "tmd"
region3 = "current"
for file in range(num_files):
    for i in range(7):
        colzaffinity[i] += calculator(zarray_t[file][i], region, "z", zbins[i])
        colxaffinity[i] += calculator(xarray_t[file][i], region, "x", xbins[i])
        colMhaffinity[i] += calculator(Mharray_t[file][i], region, "Mh")
        
        TMDzaffinity[i] += calculator(zarray_t[file][i], region2, "z", zbins[i])
        TMDxaffinity[i] += calculator(xarray_t[file][i], region2, "x", xbins[i])
        TMDMhaffinity[i] += calculator(Mharray_t[file][i], region2, "Mh")
        
        Currentzaffinity[i] += calculator(zarray_t[file][i], region3, "z", zbins[i])
        Currentxaffinity[i] += calculator(xarray_t[file][i], region3, "x", xbins[i])
        CurrentMhaffinity[i] += calculator(Mharray_t[file][i], region3, "Mh")
    for i in range(9):
        CurrentqTdivQaffinity[i] += calculator(qTdivQarray_t[file][i], region3, "qTdivQ")
        TMDqTdivQaffinity[i] += calculator(qTdivQarray_t[file][i], region2, "qTdivQ")
        colqTdivQaffinity[i] += calculator(qTdivQarray_t[file][i], region, "qTdivQ")
for i in range(9):
    print(f"TMDqTdivQaffinity[{i}] total: {TMDqTdivQaffinity[i]}")
#Now to average the affinity across all files (prob do this earlier)
for i in range(7):
    colzaffinity[i] = colzaffinity[i] / num_files
    colxaffinity[i] = colxaffinity[i] / num_files
    colMhaffinity[i] = colMhaffinity[i] / num_files
    
    TMDzaffinity[i] = TMDzaffinity[i] / num_files
    TMDxaffinity[i] = TMDxaffinity[i] / num_files
    TMDMhaffinity[i] = TMDMhaffinity[i] / num_files 
    
    Currentzaffinity[i] = Currentzaffinity[i] / num_files 
    Currentxaffinity[i] = Currentxaffinity[i] / num_files 
    CurrentMhaffinity[i] = CurrentMhaffinity[i] / num_files
for i in range(9):
    CurrentqTdivQaffinity[i] = CurrentqTdivQaffinity[i] / num_files
    print(f"TMDqTdivQaffinity[{i}] total: {TMDqTdivQaffinity[i]}; num_files: {num_files}\n div: {TMDqTdivQaffinity[i] / num_files}\n\n")
    TMDqTdivQaffinity[i] = TMDqTdivQaffinity[i] / num_files 
    colqTdivQaffinity[i] = colqTdivQaffinity[i] / num_files
    
    
# fig, ((ax1, ax2),(ax3, ax4)) = plot.subplots(2, 2, figsize = (12,12),dpi=60)
# fig.suptitle("Dihadron Affinity in the Collinear region")
# ax1.set(ylabel = "Affinity")
# ax1.scatter(Mhbins, colMhaffinity)
# ax1.axhline(y=0, color="gray", lw = 1)
# ax1.set_title("Mh binning")
# ax1.set(xlabel = "Mh (GeV)")
# ax2.scatter(xbins, colxaffinity)
# ax2.axhline(y=0, color="gray", lw = 1)
# ax2.set_title("x binning")
# ax2.set(xlabel = "x")
# ax3.scatter(zbins, colzaffinity)
# ax3.axhline(y=0, color="gray", lw = 1)
# ax3.set_title("z_h binning")
# ax3.set(xlabel = "z_h")
# fig.savefig("/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Analysis/Affinity_plots/col.jpeg")

fig2, ((ax12, ax22),(ax32,ax42)) = plot.subplots(2, 2, figsize = (12,12),dpi=60)
fig2.suptitle("Dihadron Affinity in the TMD region")
ax12.set(ylabel = "Affinity")
ax12.scatter(Mhbins, TMDMhaffinity)
ax12.axhline(y=0, color="gray", lw = 1)
ax12.set_title("Mh binning")
ax12.set(xlabel = "Mh (GeV)")
ax22.scatter(xbins, TMDxaffinity)
ax22.axhline(y=0, color="gray", lw = 1)
ax22.set_title("x binning")
ax22.set(xlabel = "x")
ax32.scatter(zbins, TMDzaffinity)
ax32.axhline(y=0, color="gray", lw = 1)
ax32.set_title("z_h binning")
ax32.set(xlabel = "z_h")
# print(f"x: {len(qTdivQbins)}; y: {len(TMDqTdivQaffinity)}")
for i in range(9):
    print(f"x: {qTdivQbins[i]}; y = {TMDqTdivQaffinity[i]}")
ax42.scatter(qTdivQbins, TMDqTdivQaffinity)
ax42.axhline(y=0, color="gray", lw = 1)
ax42.set_title("q_T/Q binning")
ax42.set(xlabel = "q_T/Q")
fig2.savefig("/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Analysis/PlotAffinityCalc/April_8/TMD.svg")

# fig3, (ax13, ax23, ax33) = plot.subplots(1, 3, figsize = (20, 5))
# fig3.suptitle("Dihadron Affinity in the Current region")
# ax13.set(ylabel = "Affinity")
# ax13.scatter(Mhbins, CurrentMhaffinity)
# ax13.axhline(y=0, color="gray", lw = 1)
# ax13.set_title("Mh binning")
# ax13.set(xlabel = "Mh (GeV)")
# ax23.scatter(xbins, Currentxaffinity)
# ax23.axhline(y=0, color="gray", lw = 1)
# ax23.set_title("x binning")
# ax23.set(xlabel = "x")
# ax33.scatter(zbins, Currentzaffinity)
# ax33.axhline(y=0, color="gray", lw = 1)
# ax33.set_title("z_h binning")
# ax33.set(xlabel = "z_h")
# fig3.savefig("/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Analysis/Affinity_plots/Cur.jpeg")