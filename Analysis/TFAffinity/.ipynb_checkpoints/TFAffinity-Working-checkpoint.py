#As of March 20th, 2023 this program works -Rowan
#This works for calculating affinity using the tensor flow model
#However, the root files must have 3 ttrees, one for each bin variable
#Additionally, although the trees don't need to contain the ratios for use here, for some reason everything
#breaks a bit when LundAnalysis.C is run and the bin variables don't save the ratios...
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
import time
logging.disable(logging.WARNING) 
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import tensorflow as tf
# print("tf.__version__", tf.__version__)

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

tmd_model_name = './models_copy/final_%s' % TMD_region_name
tmd_model = tf.keras.models.load_model(tmd_model_name)
target_model_name = './models_copy/final_%s' % target_region_name
target_model = tf.keras.models.load_model(target_model_name)
collinear_model_name = './models_copy/final_%s' % collinear_region_name
collinear_model = tf.keras.models.load_model(collinear_model_name)
current_model_name = './models_copy/final_%s' % current_region_name
current_model = tf.keras.models.load_model(current_model_name)
soft_model_name = './models_copy/final_%s' % soft_region_name
soft_model = tf.keras.models.load_model(soft_model_name)

fileDirectory = "../../OutputFiles/Slurm/May_31/Run_1/"
# fileDirectory = "../../OutputFiles/Slurm/May_25/qTQ_hadron/"
fileCount = 0
#counting number of files in target directory
for path in os.listdir(fileDirectory):
    # check if current path is a file
    if os.path.isfile(os.path.join(fileDirectory, path)):
        fileCount += 1
count = 0
#0 is z, 1 is x, 2 is mh
zarray = [[0 for i in range(3)] for j in range(7)]
xarray = [[0 for i in range(3)] for j in range(7)]
Mharray = [[0 for i in range(4)] for j in range(7)]
Q2array = [[0 for i in range(3)] for j in range(8)]
qTQarray = [[0 for i in range(4)] for j in range(7)]


xkinematics = np.array(["z_h", "Q2", "pT"])
zkinematics = np.array(["x", "Q2", "pT"])
Mhkinematics = np.array(["x", "z_h", "Q2", "pT"])
Q2kinematics = np.array(["x", "z_h", "pT"])
qTQkinematics = np.array(["x", "z_h", "Q2", "pT"])
#loop over files
for path in os.listdir(fileDirectory):
    if os.path.isfile(os.path.join(fileDirectory, path)):
        count += 1;
        inFileName = os.path.join(fileDirectory, path)
        inFile = ROOT.TFile.Open(inFileName,"READ")
        
        #grab trees
        tree_z_h_bins = inFile.Get("tree_z_h_bins")
        tree_x_bins = inFile.Get("tree_x_bins")
        tree_Mh_bins = inFile.Get("tree_Mh_bins")
        tree_Q2_bins = inFile.Get("tree_Q2_bins")
        tree_qTQ_bins = inFile.Get("tree_qTQ_bins")
#         print("On file #%d" % (count))
        try:
            #Use binnum as an interator over the number of bins in tree
            for binnum in range(0, tree_z_h_bins.GetEntries()):
                #put tree on current entry
                tree_z_h_bins.GetEntry(binnum)
                #loop over all variables (and use varnum as the iterator)
                for varnum in range(0, len(zkinematics)):
                    #sum up the values from each file
                    zarray[binnum][varnum] += getattr(tree_z_h_bins, zkinematics[varnum])

            for binnum in range(0, tree_x_bins.GetEntries()):
                tree_x_bins.GetEntry(binnum)
                for varnum in range(0, len(xkinematics)):
                    xarray[binnum][varnum] += getattr(tree_x_bins, xkinematics[varnum])

            for binnum in range(0, tree_Mh_bins.GetEntries()):
                tree_Mh_bins.GetEntry(binnum)
                for varnum in range(0, len(Mhkinematics)):
                    Mharray[binnum][varnum] += getattr(tree_Mh_bins, Mhkinematics[varnum])
            #Q2
            for binnum in range(0, tree_Q2_bins.GetEntries()):
                tree_Q2_bins.GetEntry(binnum)
                for varnum in range(0, len(Q2kinematics)):
                    Q2array[binnum][varnum] += getattr(tree_Q2_bins, Q2kinematics[varnum])
            #qTQ
            for binnum in range(0, tree_qTQ_bins.GetEntries()):
                tree_qTQ_bins.GetEntry(binnum)
                for varnum in range(0, len(qTQkinematics)):
                    qTQarray[binnum][varnum] += getattr(tree_qTQ_bins, qTQkinematics[varnum])
        except AttributeError as e:
            print("Skipping file #%d - Encountered Attribute Error" % count)
            print(e)
            continue
            
for binnum in range(0, len(zarray)):
    for varnum in range(0, len(zkinematics)):
        #calculate mean - already have sum so need to divide by number of files
        zarray[binnum][varnum] = zarray[binnum][varnum] / fileCount
for binnum in range(0, len(xarray)):
    for varnum in range(0, len(xkinematics)):
        xarray[binnum][varnum] = xarray[binnum][varnum] / fileCount
for binnum in range(0, len(Mharray)):
    for varnum in range(0, len(Mhkinematics)):
        Mharray[binnum][varnum] = Mharray[binnum][varnum] / fileCount
for binnum in range(0, len(Q2array)):
    for varnum in range(0, len(Q2kinematics)):
        Q2array[binnum][varnum] = Q2array[binnum][varnum] / fileCount
for binnum in range(0, len(qTQarray)):
    for varnum in range(0, len(qTQkinematics)):
        qTQarray[binnum][varnum] = qTQarray[binnum][varnum] / fileCount

Mhbins = np.linspace(0.3,1.3,7)
xbins = np.array([0.1,0.13,0.16,0.19,0.235,0.3,0.5])
zbins = np.array([0.35,0.43,0.49,0.55,0.62,0.7,0.83])
qTQbins = np.linspace(0.1,0.7,7)
Q2bins = np.array([1,1.4,2,2.8,4,5.6,7.9,11.1])

varName = np.array(["x", "z", "Q2", "pT", "R0max", "R1max", "R2max"])

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
    elif binType == "Mh":
        x = array[0]
        z = array[1]
        Q2 = array[2]
        pT = array[3]
    elif binType == "qTQ":
        x = array[0]
        z = array[1]
        Q2 = array[2]
        pT = array[3]
    elif binType == "Q2":
        x = array[0]
        z = array[1]
        pT = array[2]
        Q2 = binnedVariable
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
    
colxaffinityplus = np.zeros(7)
colzaffinityplus = np.zeros(7)
colMhaffinityplus = np.zeros(7)
colqTQaffinityplus = np.zeros(7)
colQ2affinityplus = np.zeros(8)
TMDxaffinityplus = np.zeros(7)
TMDzaffinityplus = np.zeros(7)
TMDMhaffinityplus = np.zeros(7)
TMDqTQaffinityplus = np.zeros(7)
TMDQ2affinityplus = np.zeros(8)
Currentxaffinityplus = np.zeros(7)
Currentzaffinityplus = np.zeros(7)
CurrentMhaffinityplus = np.zeros(7)
CurrentqTQaffinityplus = np.zeros(7)
CurrentQ2affinityplus = np.zeros(8)

region = "collinear"
region2 = "tmd"
region3 = "current"
for i in range(7):
    colzaffinityplus[i] = calculator(zarray[i], region, "z", zbins[i])
    colxaffinityplus[i] = calculator(xarray[i], region, "x", xbins[i])
    colMhaffinityplus[i] = calculator(Mharray[i], region, "Mh")
    colqTQaffinityplus[i] = calculator(qTQarray[i], region, "qTQ")
    TMDzaffinityplus[i] = calculator(zarray[i], region2, "z", zbins[i])
    TMDxaffinityplus[i] = calculator(xarray[i], region2, "x", xbins[i])
    TMDMhaffinityplus[i] = calculator(Mharray[i], region2, "Mh")
    TMDqTQaffinityplus[i] = calculator(qTQarray[i], region2, "qTQ")
    Currentzaffinityplus[i] = calculator(zarray[i], region3, "z", zbins[i])
    Currentxaffinityplus[i] = calculator(xarray[i], region3, "x", xbins[i])
    CurrentMhaffinityplus[i] = calculator(Mharray[i], region3, "Mh")
    CurrentqTQaffinityplus[i] = calculator(qTQarray[i], region3, "qTQ")
for i in range(8):
    CurrentQ2affinityplus[i] = calculator(Q2array[i], region3, "Q2", Q2bins[i])
    TMDQ2affinityplus[i] = calculator(Q2array[i], region2, "Q2", Q2bins[i])
#     print(Q2array[i])
    colQ2affinityplus[i] = calculator(Q2array[i], region, "Q2", Q2bins[i])
#     print(colzaffinityplus[i])
#     print(colxaffinityplus[i])
#     print(colMhaffinityplus[i])
#     print(TMDzaffinityplus[i])
#     print(TMDxaffinityplus[i])
#     print(TMDMhaffinityplus[i])
#     print(Currentzaffinityplus[i]) 
#     print(Currentxaffinityplus[i]) 
#     print(CurrentMhaffinityplus[i])


# fig2, ax42 = plot.subplots(1, 1, figsize = (10, 10), dpi=60)
fig2, ((ax12, ax22),(ax42, ax32)) = plot.subplots(2, 2, figsize = (12, 12), dpi=60)
fig2.suptitle("MODEL Pi+ Affinity in the TMD region")
ax12.set(ylabel = "Affinity")
ax12.scatter(Mhbins, TMDMhaffinityplus)
ax12.axhline(y=0, color="gray", lw = 1)
ax12.set_title("Mh binning")
ax12.set(xlabel = "Mh (GeV)")
ax22.scatter(xbins, TMDxaffinityplus)
ax22.axhline(y=0, color="gray", lw = 1)
ax22.set_title("x binning")
ax22.set(xlabel = "x")
ax32.scatter(zbins, TMDzaffinityplus)
ax32.axhline(y=0, color="gray", lw = 1)
ax32.set_title("z_h binning")
ax32.set(xlabel = "z_h")
ax42.scatter(qTQbins, TMDqTQaffinityplus)
ax42.axhline(y=0, color="gray", lw = 1)
ax42.set_title("qTQ binning")
ax42.set(xlabel = "qTQ")
# ax52.scatter(Q2bins, TMDQ2affinityplus)
# ax52.axhline(y=0, color="gray", lw = 1)
# ax52.set_title("Q2 binning")
# ax52.set(xlabel = "Q2")
plot.savefig("Plots/May_31_plots/TF_Affinity.jpeg")