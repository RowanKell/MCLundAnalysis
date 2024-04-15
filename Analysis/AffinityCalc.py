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

# simple version for working with CWD
file_dir = "/work/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/April_14/Run_1_single_pion/"
num_files = len([name for name in os.listdir(file_dir) if not os.path.isdir(name)])
file_names = [name for name in os.listdir(file_dir) if not os.path.isdir(name)]
# file_dir = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/April_11/Run_1/"
# num_files = 1
# file_names = ["file_1.root"]

product = False
single_pion = True
tree_MC_list = []
tree_x_list = []
tree_z_h_list = [] #Notes - for product, the bin z_h is the dihadron z_h, while the z we need for affinity is z_1
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
xarray = np.array([np.array([np.zeros(7)] * 3)] * num_files)
Mharray = np.array([np.array([np.zeros(7)] * 4)] * num_files)
qTdivQarray = np.array([np.array([np.zeros(9)] * 4)] * num_files)
if(product):
    zarray = np.array([np.array([np.zeros(7)] * 4)] * num_files)#for product, need to include z for pion in array
    zarray_2 = np.array([np.array([np.zeros(7)] * 4)] * num_files)
    xarray_2 = np.array([np.array([np.zeros(7)] * 3)] * num_files)
    Mharray_2 = np.array([np.array([np.zeros(7)] * 4)] * num_files)
    qTdivQarray_2 = np.array([np.array([np.zeros(9)] * 4)] * num_files)
else:
    zarray = np.array([np.array([np.zeros(7)] * 3)] * num_files)



if(product):
    xkinematics = np.array(["z_h_1", "Q2", "pT_1"])
    zkinematics = np.array(["x", "Q2", "pT_1","z_h_1"])
    Mhkinematics = np.array(["x", "z_h_1", "Q2", "pT_1"])
    qTdivQkinematics = np.array(["x", "z_h_1", "Q2", "pT_1"])
    
    xkinematics_2 = np.array(["z_h_2", "Q2", "pT_2"])
    zkinematics_2 = np.array(["x", "Q2", "pT_2","z_h_2"])
    Mhkinematics_2 = np.array(["x", "z_h_2", "Q2", "pT_2"])
    qTdivQkinematics_2 = np.array(["x", "z_h_2", "Q2", "pT_2"])
elif(single_pion):
    xkinematics = np.array(["z_h_1", "Q2", "pT_1"])
    zkinematics = np.array(["x", "Q2", "pT_1"])
    Mhkinematics = np.array(["x", "z_h_1", "Q2", "pT_1"])
    qTdivQkinematics = np.array(["x", "z_h_1", "Q2", "pT_1"])
else:#dihadron case
    xkinematics = np.array(["z_h", "Q2", "pT"])
    zkinematics = np.array(["x", "Q2", "pT"])
    Mhkinematics = np.array(["x", "z_h", "Q2", "pT"])
    qTdivQkinematics = np.array(["x", "z_h", "Q2", "pT"])
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
    if(product):
        z_iter = 0
        for var in zkinematics_2:
            zarray_2[i][z_iter] = tree_z_h_list[i][var].array(library='np')
            z_iter += 1
        x_iter = 0
        for var in xkinematics_2:
            xarray_2[i][x_iter] = tree_x_list[i][var].array(library='np')
            x_iter += 1
        Mh_iter = 0
        for var in Mhkinematics_2:
            Mharray_2[i][Mh_iter] = tree_Mh_list[i][var].array(library='np')
            Mh_iter += 1
        qTdivQ_iter = 0
        for var in qTdivQkinematics_2:
            qTdivQarray_2[i][qTdivQ_iter] = tree_qTdivQ_list[i][var].array(library='np')
            qTdivQ_iter += 1

xarray_t = np.array([np.array([np.zeros(3)] * 7)] * num_files)
Mharray_t = np.array([np.array([np.zeros(4)] * 7)] * num_files)
qTdivQarray_t = np.array([np.array([np.zeros(4)] * 9)] * num_files)
if(product):
    zarray_t = np.array([np.array([np.zeros(4)] * 7)] * num_files)
    zarray_t_2 = np.array([np.array([np.zeros(4)] * 7)] * num_files)
    xarray_t_2 = np.array([np.array([np.zeros(3)] * 7)] * num_files)
    Mharray_t_2 = np.array([np.array([np.zeros(4)] * 7)] * num_files)
    qTdivQarray_t_2 = np.array([np.array([np.zeros(4)] * 9)] * num_files)
else:
    zarray_t = np.array([np.array([np.zeros(3)] * 7)] * num_files)
for i in range(num_files):
    xarray_t[i] = np.transpose(xarray[i])
    zarray_t[i] = np.transpose(zarray[i])
    Mharray_t[i] = np.transpose(Mharray[i])
    qTdivQarray_t[i] = np.transpose(qTdivQarray[i])
    if(product):
        xarray_t_2[i] = np.transpose(xarray[i])
        zarray_t_2[i] = np.transpose(zarray[i])
        Mharray_t_2[i] = np.transpose(Mharray[i])
        qTdivQarray_t_2[i] = np.transpose(qTdivQarray[i])
TMD_region_name = 'TMD'
TMD_lable_name = 'tmdaff'

tmd_model_name = '../../SIDIS-Affinity/models/final_%s' % TMD_region_name
tmd_model = tf.keras.models.load_model(tmd_model_name)

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
        if(product):
            z = array[3]
        else:
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
TMDxaffinity = np.zeros(7)
TMDzaffinity = np.zeros(7)
TMDMhaffinity = np.zeros(7)
TMDqTdivQaffinity = np.zeros(9)

region = "collinear"
region2 = "tmd"
region3 = "current"
#For product affinity - multiply output of calculator for both pions
if(product):
    print("entering product")
    for file in range(num_files):
        for i in range(7):
            TMDzaffinity[i] += calculator(zarray_t[file][i], region2, "z", zbins[i]) * calculator(zarray_t_2[file][i], region2, "z", zbins[i])
            TMDxaffinity[i] += calculator(xarray_t[file][i], region2, "x", xbins[i]) * calculator(xarray_t_2[file][i], region2, "x", xbins[i])
            TMDMhaffinity[i] += calculator(Mharray_t[file][i], region2, "Mh") * calculator(Mharray_t_2[file][i], region2, "Mh")
        for i in range(9):
            TMDqTdivQaffinity[i] += calculator(qTdivQarray_t[file][i], region2, "qTdivQ") * calculator(qTdivQarray_t_2[file][i], region2, "qTdivQ")
else:
    for file in range(num_files):
        for i in range(7):
            TMDzaffinity[i] += calculator(zarray_t[file][i], region2, "z", zbins[i])
            TMDxaffinity[i] += calculator(xarray_t[file][i], region2, "x", xbins[i])
            TMDMhaffinity[i] += calculator(Mharray_t[file][i], region2, "Mh")
        for i in range(9):
            TMDqTdivQaffinity[i] += calculator(qTdivQarray_t[file][i], region2, "qTdivQ")

#Now to average the affinity across all files (prob do this earlier)
for i in range(7):
    TMDzaffinity[i] = TMDzaffinity[i] / num_files
    TMDxaffinity[i] = TMDxaffinity[i] / num_files
    TMDMhaffinity[i] = TMDMhaffinity[i] / num_files 
for i in range(9):
    TMDqTdivQaffinity[i] = TMDqTdivQaffinity[i] / num_files 

fig2, ((ax12, ax22),(ax32,ax42)) = plot.subplots(2, 2, figsize = (10,10),dpi=60)
fig2.suptitle("Single Pion Affinity in the TMD region")
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
ax42.scatter(qTdivQbins, TMDqTdivQaffinity)
ax42.axhline(y=0, color="gray", lw = 1)
ax42.set_title("q_T/Q binning")
ax42.set(xlabel = "q_T/Q")
fig2.savefig("/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/Analysis/PlotAffinityCalc/April_14/TMD_one_file_single_pion.svg")