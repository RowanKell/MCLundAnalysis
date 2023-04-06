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
for i in range(7):
    pxcut[i] = d_plus.Filter(xformat.format(xbins[i + 1],xbins[i])).Filter("R0max < 0.3").Filter("R2max < 0.3").Filter("R1max < 0.3").Count()
    px[i] = d_plus.Filter(xformat.format(xbins[i + 1],xbins[i])).Count()
    pMhcut[i] = d_plus.Filter(Mhformat.format(Mhbins[i + 1],Mhbins[i])).Filter("R0max < 0.3").Filter("R2max < 0.3").Filter("R1max < 0.3").Count()
    pMh[i] = d_plus.Filter(Mhformat.format(Mhbins[i + 1],Mhbins[i])).Count()
    pzcut[i] = d_plus.Filter(zformat.format(zbins[i + 1],zbins[i])).Filter("R0max < 0.3").Filter("R2max < 0.3").Filter("R1max < 0.3").Count() 
    pz[i] = d_plus.Filter(zformat.format(zbins[i + 1],zbins[i])).Count()
    ppTcut[i] = d_plus.Filter(pTformat.format(pTbins[i + 1],pTbins[i])).Filter("R0max < 0.3").Filter("R2max < 0.3").Filter("R1max < 0.3").Count()
    ppT[i] = d_plus.Filter(pTformat.format(pTbins[i + 1],pTbins[i])).Count()
    pqcut[i] = d_plus.Filter(qdivformat.format(qTbins[i + 1],qTbins[i])).Filter("R0max < 0.3").Filter("R2max < 0.3").Filter("R1max < 0.3").Count()
    pq[i] = d_plus.Filter(qdivformat.format(qTbins[i + 1],qTbins[i])).Count()
        
        
for i in range(7):
    pxcut[i] = pxcut[i].GetValue()
    px[i] = px[i].GetValue()
    pMhcut[i] = pMhcut[i].GetValue()
    pMh[i] = pMh[i].GetValue()
    pzcut[i] = pzcut[i].GetValue()
    pz[i] = pz[i].GetValue()
    ppTcut[i] = ppTcut[i].GetValue()
    ppT[i] = ppT[i].GetValue()
    pqcut[i] = pqcut[i].GetValue()
    pq[i] = pq[i].GetValue()
    
    
pxval = [0 for i in range(7)]
pzval = [0 for i in range(7)]
pMhval=[0 for i in range(7)]
ppTval=[0 for i in range(7)]
pqval = [0 for i in range(7)]
Mhbinsno0 = np.linspace(0.3,1.3,7)
pTbinsno0 = np.linspace(0.2,0.8,7)
xbinsno0 = np.array([0.1,0.13,0.16,0.19,0.235,0.3,0.5])
zbinsno0 = np.array([0.35,0.43,0.49,0.55,0.62,0.7,0.83])
qTbinsno0 = np.linspace(0.1,0.7,7)


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
    if(pq[i] == 0):
        pqval[i] = 0
    else: pqval[i] = pqcut[i] / pq[i]
    
    
fig2, (ax42, ax22, ax32) = plot.subplots(1, 3, figsize = (15, 4), dpi=60)
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
# ax42.scatter(qTbinsno0, pqval, c = 'r', marker = "+")
# ax42.axhline(y=0, color="gray", lw = 1)
# ax42.set_title("qTdivQ binning")
# ax42.set(xlabel = "qTdivQ")
plot.savefig("Plots/Box3.jpeg")
