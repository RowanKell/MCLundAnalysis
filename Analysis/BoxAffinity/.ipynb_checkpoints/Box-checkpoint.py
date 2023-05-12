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

d_plus = RDataFrame("tree_MC", "../../OutputFiles/Slurm/April_21/Run_1/file_*.root")

#Bins (each has 8 including 0)
Mhbins = np.linspace(0,1.3,8)
pTbins = np.linspace(0.1,0.8,8)
xbins = np.array([0,0.1,0.13,0.16,0.19,0.235,0.3,0.5])
zbins = np.array([0,0.35,0.43,0.49,0.55,0.62,0.7,0.83])
qTQbins = np.linspace(0,0.7,8)
Q2bins = np.array([0,1,1.4,2,2.8,4,5.6,7.9,11.1])

varName = np.array(["x", "z", "Q2", "pT", "R0max", "R1max", "R2max"])

px = [0 for i in range(7)]
pz = [0 for i in range(7)]
pMh=[0 for i in range(7)]
ppT=[0 for i in range(7)]
pQ2 = [0 for i in range(8)]
pqTQ = [0 for i in range(7)]

pxcut = [0 for i in range(7)]
pzcut = [0 for i in range(7)]
pMhcut=[0 for i in range(7)]
ppTcut=[0 for i in range(7)]
pQ2cut=[0 for i in range(8)]
pqTQcut=[0 for i in range(7)]

# start_time = time.time()
xformat = "x <= {} && x > {}"
zformat = "z <= {} && z > {}"
Mhformat = "Mh <= {} && Mh > {}"
pTformat = "pT <= {} && pT > {}"
R0format = "R0max <= {} && R0max > {}"
R1format = "R1max <= {} && R1max > {}"
R2format = "R2max <= {} && R2max > {}"
qdivformat = "q_TdivQ <= {} && q_TdivQ > {}"
Q2format = "Q2 <= {} && Q2 > {}"

#Piplus
#x bins
#i is the kinematic variable
#j is the bin num
for i in range(7):
    pxcut[i] = d_plus.Filter(xformat.format(xbins[i + 1],xbins[i])).Filter("R0 < 0.3").Filter("R2 < 0.3").Filter("R1 < 0.3").Count()
    px[i] = d_plus.Filter(xformat.format(xbins[i + 1],xbins[i])).Count()
    pMhcut[i] = d_plus.Filter(Mhformat.format(Mhbins[i + 1],Mhbins[i])).Filter("R0 < 0.3").Filter("R2 < 0.3").Filter("R1 < 0.3").Count()
    pMh[i] = d_plus.Filter(Mhformat.format(Mhbins[i + 1],Mhbins[i])).Count()
    pzcut[i] = d_plus.Filter(zformat.format(zbins[i + 1],zbins[i])).Filter("R0 < 0.3").Filter("R2 < 0.3").Filter("R1 < 0.3").Count() 
    pz[i] = d_plus.Filter(zformat.format(zbins[i + 1],zbins[i])).Count()
    ppTcut[i] = d_plus.Filter(pTformat.format(pTbins[i + 1],pTbins[i])).Filter("R0 < 0.3").Filter("R2 < 0.3").Filter("R1 < 0.3").Count()
    ppT[i] = d_plus.Filter(pTformat.format(pTbins[i + 1],pTbins[i])).Count()
    pqTQ[i] = d_plus.Filter(qdivformat.format(qTQbins[i + 1],qTQbins[i])).Count()
    pqTQcut[i] = d_plus.Filter(qdivformat.format(qTQbins[i + 1],qTQbins[i])).Filter("R0 < 0.3").Filter("R2 < 0.3").Filter("R1 < 0.3").Count()
for i in range(8):
    pQ2[i] = d_plus.Filter(Q2format.format(Q2bins[i + 1],Q2bins[i])).Count()
    pQ2cut[i] = d_plus.Filter(Q2format.format(Q2bins[i + 1],Q2bins[i])).Filter("R0 < 0.3").Filter("R2 < 0.3").Filter("R1 < 0.3").Count()
        
for i in range(7):
    pxcut[i] = pxcut[i].GetValue()
    px[i] = px[i].GetValue()
    pMhcut[i] = pMhcut[i].GetValue()
    pMh[i] = pMh[i].GetValue()
    pzcut[i] = pzcut[i].GetValue()
    pz[i] = pz[i].GetValue()
    ppTcut[i] = ppTcut[i].GetValue()
    ppT[i] = ppT[i].GetValue()
    pqTQcut[i] = pqTQcut[i].GetValue()
    pqTQ[i] = pqTQ[i].GetValue()
for i in range(8):
    pQ2cut[i] = pQ2cut[i].GetValue()
    pQ2[i] = pQ2[i].GetValue()
    
    
pxval = [0 for i in range(7)]
pzval = [0 for i in range(7)]
pMhval=[0 for i in range(7)]
ppTval=[0 for i in range(7)]
pqTQval = [0 for i in range(7)]
pQ2val = [0 for i in range(8)]
Mhbinsno0 = np.linspace(0.3,1.3,7)
pTbinsno0 = np.linspace(0.2,0.8,7)
xbinsno0 = np.array([0.1,0.13,0.16,0.19,0.235,0.3,0.5])
zbinsno0 = np.array([0.35,0.43,0.49,0.55,0.62,0.7,0.83])
qTQbinsno0 = np.linspace(0.1,0.7,7)
Q2binsno0 = np.array([1,1.4,2,2.8,4,5.6,7.9,11.1])


for i in range(7):
    if(px[i] == 0):
        pxval[i] = 0
    else: pxval[i] = pxcut[i] / px[i]
    if(pz[i] == 0):
        pzval[i] = 0
    else: pzval[i] = pzcut[i] / pz[i]
    if(ppT[i] == 0):
        ppTval[i] = 0
    else: ppTval[i] = ppTcut[i] / ppT[i]
    if(pMh[i] == 0):
        pMhval[i] = 0
    else: pMhval[i] = pMhcut[i] / pMh[i]
    if(pqTQ[i] == 0):
        pqTQval[i] = 0
    else: pqTQval[i] = pqTQcut[i] / pqTQ[i]
for i in range(8):
    if(pQ2[i] == 0):
        pQ2val[i] = 0
    else: pQ2val[i] = pQ2cut[i] / pQ2[i]
    
    
fig2, ((ax42, ax22, ax32), (ax52, ax62, ax72)) = plot.subplots(2, 3, figsize = (15, 10), dpi=60)
# fig2, ((ax42, ax22, ax32), (ax52, ax62, ax72), (ax82, ax92, ax02)) = plot.subplots(3, 3, figsize = (15, 12), dpi=60)
fig2.suptitle("BOX Pi+ Affinity in the TMD region: Rmax = 0.3")
# fig2.ylable("Affinity")
# ax12.set(ylabel = "Affinity")
# ax12.scatter(pTbinsno0, ppTval, c = 'r', marker = "+")
# # fig2.legend(shadow = "true", title = "Key", fontsize = 12)
# ax12.axhline(y=0, color="gray", lw = 1)
# ax12.set_title("pT binning")
# ax12.set(xlabel = "pT (GeV)")
ax22.scatter(xbinsno0, pxval, c = 'r', marker = '+')
ax22.axhline(y=0, color="gray", lw = 1)
ax22.set_title("x binning")
ax22.set(xlabel = "x")
ax32.scatter(zbinsno0, pzval, c = 'r', marker = '+')
ax32.axhline(y=0, color="gray", lw = 1)
ax32.set_title("z_h binning")
ax32.set(xlabel = "z_h")
ax42.set(ylabel = "Affinity")
ax42.scatter(Mhbinsno0, pMhval, c = 'r', marker = "+")
ax42.axhline(y=0, color="gray", lw = 1)
ax42.set_title("Mh binning")
ax42.set(xlabel = "Mh (GeV)")
ax52.scatter(qTQbinsno0, pqTQval, c = 'r', marker = "+")
ax52.axhline(y=0, color="gray", lw = 1)
ax52.set_title("qTdivQ binning")
ax52.set(xlabel = "qTdivQ")
ax62.scatter(Q2binsno0, pQ2val, c = 'r', marker = "+")
ax62.axhline(y=0, color="gray", lw = 1)
ax62.set_title("Q2 binning")
ax62.set(xlabel = "Q2")
# ax72.scatter(qTQbinsno0, pqTQ, c = 'r', marker = "+")
# ax72.axhline(y=0, color="gray", lw = 1)
# ax72.set_title("qTQ counts")
# ax72.set(xlabel = "qTQ")
# ax82.scatter(xbinsno0, px, c = 'r', marker = "+")
# ax82.axhline(y=0, color="gray", lw = 1)
# ax82.set_title("x counts")
# ax82.set(xlabel = "x")
# ax92.scatter(zbinsno0, pz, c = 'r', marker = "+")
# ax92.axhline(y=0, color="gray", lw = 1)
# ax92.set_title("z counts")
# ax92.set(xlabel = "z")
# ax02.scatter(Q2binsno0, pQ2, c = 'r', marker = "+")
# ax02.axhline(y=0, color="gray", lw = 1)
# ax02.set_title("Q2 counts")
# ax02.set(xlabel = "Q2")
plot.savefig("Plots/Box_April_28.jpeg")
