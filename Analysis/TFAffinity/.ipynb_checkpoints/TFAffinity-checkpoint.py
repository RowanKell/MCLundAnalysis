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

d_plus = RDataFrame("tree_MC", "../../OutputFiles/Slurm/March_2/run_1/file_*.root")

#Bins (each has 8 including 0)
Mhbins = np.linspace(0,1.3,8)
pTbins = np.linspace(0.1,0.8,8)
xbins = np.array([0,0.1,0.13,0.16,0.19,0.235,0.3,0.5])
zbins = np.array([0,0.35,0.43,0.49,0.55,0.62,0.7,0.83])
Q2bins = np.array([0,1.2,1.8,2.3,3.1,4.3,7,11.1])
qTbins = np.linspace(0,0.7,8)

varName = np.array(["x", "z", "Q2", "pT", "R0max", "R1max", "R2max"])

px = [0 for i in range(7)]
pz = [0 for i in range(7)]
pMh=[0 for i in range(7)]
ppT=[0 for i in range(7)]
pxQ2 = [0 for i in range(7)]
pq = [0 for i in range(7)]

pxcut = [0 for i in range(7)]
pzcut = [0 for i in range(7)]
pMhcut=[0 for i in range(7)]
ppTcut=[0 for i in range(7)]
pqcut=[0 for i in range(7)]

# start_time = time.time()
xformat = "x <= {} && x > {}"
zformat = "z <= {} && z > {}"
Mhformat = "Mh <= {} && Mh > {}"
pTformat = "pT <= {} && pT > {}"
R0format = "R0max <= {} && R0max > {}"
R1format = "R1max <= {} && R1max > {}"
R2format = "R2max <= {} && R2max > {}"
qdivformat = "q_TdivQ <= {} && q_TdivQ > {}"

#Piplus
#x bins
#i is the kinematic variable
#j is the bin num
d_plus_x = [0 for i in range(7)]
d_plus_Mh = [0 for i in range(7)]
d_plus_z =[0 for i in range(7)]
d_plus_z = [0 for i in range(7)]
d_plus_pT = [0 for i in range(7)]
d_plus_qT = [0 for i in range(7)]

pxcut = [0 for i in range(7)]
pzcut = [0 for i in range(7)]
pMhcut=[0 for i in range(7)]
ppTcut=[0 for i in range(7)]
pqcut=[0 for i in range(7)]
for i in range(7):
    d_plus_x[i] = d_plus.Filter(xformat.format(xbins[i + 1],xbins[i]))
    d_plus_Mh[i] = d_plus.Filter(Mhformat.format(Mhbins[i + 1],Mhbins[i]))
    d_plus_z[i] = d_plus.Filter(zformat.format(zbins[i + 1],zbins[i]))
    d_plus_pT[i] = d_plus.Filter(pTformat.format(pTbins[i + 1],pTbins[i]))
    d_plus_qT[i] = d_plus.Filter(qdivformat.format(qTbins[i + 1],qTbins[i]))
xarr =[0 for i in range(7)]
Mharr =[0 for i in range(7)]
zarr =[0 for i in range(7)]
pTarr =[0 for i in range(7)]
qTarr =[0 for i in range(7)]
for i in range(7):
    xarr[i] = (d_plus_x[i].AsNumpy(columns=["x", "z", "Q2", "pT"]))
    Mharr[i] = (d_plus_Mh[i].AsNumpy(columns=["x", "z", "Q2", "pT"]))
    zarr[i] = (d_plus_z[i].AsNumpy(columns=["x", "z", "Q2", "pT"]))
    pTarr[i] = (d_plus_pT[i].AsNumpy(columns=["x", "z", "Q2", "pT"]))
    qTarr[i] = (d_plus_qT[i].AsNumpy(columns=["x", "z", "Q2", "pT"]))
Mhbins = np.linspace(0.3,1.3,7)
xbins = np.array([0.1,0.13,0.16,0.19,0.235,0.3,0.5])
zbins = np.array([0.35,0.43,0.49,0.55,0.62,0.7,0.83])
def calculator(array, region, binType, binnedVariable = 0):
    R0max = 0.3
    R1max = 0.3
    R2max = 0.3
    if binType == "x":
        z = array[1]
        Q2 = array[2]
        pT = array[3]
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
TMDxaffinityplus = np.zeros(7)
TMDzaffinityplus = np.zeros(7)
TMDMhaffinityplus = np.zeros(7)
Currentxaffinityplus = np.zeros(7)
Currentzaffinityplus = np.zeros(7)
CurrentMhaffinityplus = np.zeros(7)

colxaffinityminus = np.zeros(7)
colzaffinityminus = np.zeros(7)
colMhaffinityminus = np.zeros(7)
TMDxaffinityminus = np.zeros(7)
TMDzaffinityminus = np.zeros(7)
TMDMhaffinityminus = np.zeros(7)
Currentxaffinityminus = np.zeros(7)
Currentzaffinityminus = np.zeros(7)
CurrentMhaffinityminus = np.zeros(7)

region = "collinear"
region2 = "tmd"
region3 = "current"
for i in range(7):
    colzaffinityplus[i] = calculator(zarr[i], region, "z", zbins[i])
    colxaffinityplus[i] = calculator(xarr[i], region, "x", xbins[i])
    colMhaffinityplus[i] = calculator(Mharr[i], region, "Mh")
    TMDzaffinityplus[i] = calculator(zarr[i], region2, "z", zbins[i])
    TMDxaffinityplus[i] = calculator(xarr[i], region2, "x", xbins[i])
    TMDMhaffinityplus[i] = calculator(Mharr[i], region2, "Mh")
    Currentzaffinityplus[i] = calculator(zarr[i], region3, "z", zbins[i])
    Currentxaffinityplus[i] = calculator(xarr[i], region3, "x", xbins[i])
    CurrentMhaffinityplus[i] = calculator(Mharr[i], region3, "Mh")
print("done")