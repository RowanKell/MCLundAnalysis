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
import sys
logging.disable(logging.WARNING) 
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import argparse
parser = argparse.ArgumentParser(description ='Remove ratio cut')
parser.add_argument('--no_cut',type=str,default="",help='cut to remove e.g. \'R0\'')
args = parser.parse_args()
no_cut = args.no_cut

dir_prefix = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/"

driver = False

calculate_product = False # set true to calculate dihadron affinity as product of single pion affinities
single_pion = True #calculate the affinity of single pion
dihadron_double_cut = False #calculate affinity by cutting on both R1_p and R1_m
dihadron_single_cut = False

if((calculate_product and single_pion) or (calculate_product and dihadron_double_cut) or (dihadron_double_cut and single_pion)):
    print(f"two affinity calc flags are set. 0 or 1 of the following flags may be set: calculate_product;single_pion;dihadron_double_cut;")
    sys.exit(1)
if(single_pion):
    print("calculating single pion box affinity")
    no_cut_text = "single_pion_box_calc_test_500k_R2_adjust"
elif(calculate_product):
    print("calculating product box affinity")
    no_cut_text = "all_files_product"
elif(dihadron_double_cut):
    print("calculating dihadron double cut box affinity")
    no_cut_text = "all_files_dihadron_double_cut"
elif(dihadron_single_cut):
    print("calculating dihadron single cut box affinity")
    no_cut_text = "all_files_dihadron_single_cut"
no_cut_file_name = no_cut_text + ".svg"


# d = RDataFrame("tree_MC",dir_prefix + "OutputFiles/Slurm_Spring_24/April_13/Run_1_dihadron/file_*.root")
# d = RDataFrame("tree_MC",dir_prefix + "OutputFiles/Slurm_Spring_24/April_13/Run_1_single_pion/file_*.root")
# d = RDataFrame("tree_MC",dir_prefix + "/OutputFiles/Files_Spring_24/April_21/Run_1_single_pion/file_0_test.root")
# d = RDataFrame("tree_MC",dir_prefix + "OutputFiles/Files_Spring_24/April_21/Run_1_dihadron/file_0_test.root")
# d = RDataFrame("tree_MC",dir_prefix + "/OutputFiles/Files_Spring_24/May_7/Run_1_single_pion/file_0_test_fix.root")
# d = RDataFrame("tree_MC",dir_prefix + "/OutputFiles/Slurm_Spring_24/May_15/Run_1_single_pion/file_0.root")
d = RDataFrame("tree_MC",dir_prefix + "/OutputFiles/Files_Spring_24/May_23/for_box_500k.root")
# d = RDataFrame("tree_driver",dir_prefix + "../sidisregions/rowan_dev/root_files/driver_ratios_May23_500k.root")
print(f"column names: {d.GetColumnNames()}")
#Bins (each has 8 including 0)
Mhbins = np.linspace(0.25,1.6,8)
pTbins = np.linspace(0.1,0.8,8)
xbins = np.array([0,0.1,0.13,0.16,0.19,0.235,0.3,0.5])
zbins = np.array([0,0.35,0.43,0.49,0.55,0.62,0.7,0.83])
qTQbins = np.array([0,0.1,0.3,0.5,0.8,1.1,1.5,2,2.5,3])

Q2bins = np.array([0,1,1.4,2,2.8,4,5.6,7.9,11.1])

Mhbinsno0 = np.delete(Mhbins, 0)
pTbinsno0 = np.linspace(0.2,0.8,7)
xbinsno0 = np.array([0.1,0.13,0.16,0.19,0.235,0.3,0.5])
zbinsno0 = np.array([0.35,0.43,0.49,0.55,0.62,0.7,0.83])
qTQbinsno0 = np.array([0.1,0.3,0.5,0.8,1.1,1.5,2,2.5,3])
Q2binsno0 = np.array([1,1.4,2,2.8,4,5.6,7.9,11.1])

varName = np.array(["x", "z", "Q2", "pT", "R0", "R1", "R2"])

def calculate_Affinity(d, calculate_product, pi_num):
    print(f"IN CALCULATE_AFFINITY: calculate_product: {calculate_product} | pi_num: {pi_num}")

    px = [0 for i in range(7)]
    pz = [0 for i in range(7)]
    pMh=[0 for i in range(7)]
    ppT=[0 for i in range(7)]
    pQ2 = [0 for i in range(8)]
    pqTQ = [0 for i in range(9)]

    pxcut = [0 for i in range(7)]
    pzcut = [0 for i in range(7)]
    pMhcut=[0 for i in range(7)]
    ppTcut=[0 for i in range(7)]
    pQ2cut=[0 for i in range(8)]
    pqTQcut=[0 for i in range(9)]
    if(driver):#if using driver, branch names are diff
        xformat = "x.x_t <= {} && x.x_t > {}"
        zformat = "z.z_t <= {} && z.z_t > {}"
        Mhformat = "Mh.Mh_t <= {} && Mh.Mh_t > {}"
        pTformat = "pT.pT_t <= {} && pT.pT_t > {}"
        R0format = "R0.R0_t < 0.3"
        R2format = "R2.R2_t < 0.3"
        if(calculate_product):
            if(pi_num == 1):
                R1format = "R1_p < 0.3"
            else:
                R1format = "R1_m < 0.3"
        elif(single_pion):
            R1format = "R1.R1_t < 0.3" #if single pion, use R1 for first pion
        elif(dihadron_double_cut):
            R1format = "R1_p < 0.3 && R1_m < 0.3" #if single pion, use R1 for first pion
        else:
            R1format = "R1.R1_t < 0.3"
        print(f"Setting R1format to {R1format}")
        qdivformat = "qTQ_hadron.qTQ_hadron_t <= {} && qTQ_hadron.qTQ_hadron_t > {}"
        Q2format = "Q2.Q2_t <= {} && Q2.Q2_t > {}"
    else:
        xformat = "x <= {} && x > {}"
        zformat = "z <= {} && z > {}"
        Mhformat = "Mh <= {} && Mh > {}"
        pTformat = "pT <= {} && pT > {}"
        R0format = "R0 < 0.3"
        R2format = "R2_adjust < 0.3"
#         R2format = "R2.R2_t < 0.3"
        if(calculate_product):
            if(pi_num == 1):
                R1format = "R1_p < 0.3"
            else:
                R1format = "R1_m < 0.3"
        elif(single_pion):
            R1format = "R1 < 0.3" #if single pion, use R1 for first pion
        elif(dihadron_double_cut):
            R1format = "R1_p < 0.3 && R1_m < 0.3" #if single pion, use R1 for first pion
        else:
            R1format = "R1 < 0.3"
        print(f"Setting R1format to {R1format}")
        qdivformat = "qTQ_hadron <= {} && qTQ_hadron > {}"
        Q2format = "Q2 <= {} && Q2 > {}"
    
    for i in range(7):
        pxcut[i] = d.Filter(xformat.format(xbins[i + 1],xbins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count()
        px[i] = d.Filter(xformat.format(xbins[i + 1],xbins[i])).Count()
        if(not single_pion): #only do Mh binning for dihadrons
            pMhcut[i] = d.Filter(Mhformat.format(Mhbins[i + 1],Mhbins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count()
            pMh[i] = d.Filter(Mhformat.format(Mhbins[i + 1],Mhbins[i])).Count()
#             print(f"Calculating Mh affinity (should be dihadron or product then)")
        pzcut[i] = d.Filter(zformat.format(zbins[i + 1],zbins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count() 
        pz[i] = d.Filter(zformat.format(zbins[i + 1],zbins[i])).Count()
        ppTcut[i] = d.Filter(pTformat.format(pTbins[i + 1],pTbins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count()
        ppT[i] = d.Filter(pTformat.format(pTbins[i + 1],pTbins[i])).Count()
    for i in range(8):
        pQ2[i] = d.Filter(Q2format.format(Q2bins[i + 1],Q2bins[i])).Count()
        pQ2cut[i] = d.Filter(Q2format.format(Q2bins[i + 1],Q2bins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count()
    for i in range(9):
        pqTQ[i] = d.Filter(qdivformat.format(qTQbins[i + 1],qTQbins[i])).Count()
        pqTQcut[i] = d.Filter(qdivformat.format(qTQbins[i + 1],qTQbins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count()

    for i in range(7):
        pxcut[i] = pxcut[i].GetValue()
        px[i] = px[i].GetValue()
        if(not single_pion):
            pMhcut[i] = pMhcut[i].GetValue()
            pMh[i] = pMh[i].GetValue()
        pzcut[i] = pzcut[i].GetValue()
        pz[i] = pz[i].GetValue()
        ppTcut[i] = ppTcut[i].GetValue()
        ppT[i] = ppT[i].GetValue()
    for i in range(8):
        pQ2cut[i] = pQ2cut[i].GetValue()
        pQ2[i] = pQ2[i].GetValue()
    for i in range(9):
        pqTQcut[i] = pqTQcut[i].GetValue()
        pqTQ[i] = pqTQ[i].GetValue()


    pxval = [0 for i in range(7)]
    pzval = [0 for i in range(7)]
    if(not single_pion):
        pMhval=[0 for i in range(7)]
    ppTval=[0 for i in range(7)]
    pqTQval = [0 for i in range(9)]
    pQ2val = [0 for i in range(8)]


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
        if(not single_pion):
            if(pMh[i] == 0):
                pMhval[i] = 0
            else: pMhval[i] = pMhcut[i] / pMh[i]
    for i in range(8):
        if(pQ2[i] == 0):
            pQ2val[i] = 0
        else: pQ2val[i] = pQ2cut[i] / pQ2[i]
    for i in range(9):
        if(pqTQ[i] == 0):
            pqTQval[i] = 0
        else: pqTQval[i] = pqTQcut[i] / pqTQ[i]
    
    if(single_pion):
        Affinity_arrs = [pxval,pzval,pqTQval,pQ2val]
        print(f"Creating array without Mh affinity")
    else:
        Affinity_arrs = [pxval,pzval,pqTQval,pQ2val,pMhval]
    return Affinity_arrs
if(calculate_product):
    print(f"beginning product affinity calculation")
    plus_arrs = calculate_Affinity(d, calculate_product, 1)
    minus_arrs = calculate_Affinity(d, calculate_product, 2)
    prod_xval = [0 for i in range(7)]
    prod_zval = [0 for i in range(7)]
    prod_Mhval=[0 for i in range(7)]
    prod_qTQval = [0 for i in range(9)]
    prod_Q2val = [0 for i in range(8)]

    product_arrs = [prod_xval,prod_zval,prod_qTQval,prod_Q2val,prod_Mhval]
    bin_lengths = [7,7,9,8,7]

    #x, z, Mh, qTQ, Q2
    for variable in range(5):
        for bin_num in range(bin_lengths[variable]):
            product_arrs[variable][bin_num] = plus_arrs[variable][bin_num] * minus_arrs[variable][bin_num]
    print(f"Finished calculating affinity arr for product")
elif(single_pion):
    print(f"beginning single pion affinity calculation")
    product_arrs = calculate_Affinity(d, calculate_product, 0)
    print(f"Finished calculating affinity arr for single pion")
else:
    print(f"beginning dihadron affinity calculation")
    product_arrs = calculate_Affinity(d, calculate_product, 0)
    print(f"Finished calculating affinity arr for dihadron")
    
fig2, ((ax42, ax22), (ax32,ax52)) = plot.subplots(2, 2, figsize = (10, 10), dpi=30)
if(calculate_product):
    fig2.suptitle("Product Dihadron BOX Affinity in the TMD region")
elif(single_pion):
    fig2.suptitle("Single Pion BOX Affinity in the TMD region")
elif(dihadron_double_cut):
    fig2.suptitle("Dihadron doublecut BOX Affinity in the TMD region")
else:
    fig2.suptitle("Dihadron BOX Affinity in the TMD region")
ax22.scatter(xbinsno0, product_arrs[0], c = 'r', marker = '+')
ax22.axhline(y=0, color="gray", lw = 1)
ax22.set_title("x binning")
ax22.set(xlabel = "x")
ax32.scatter(zbinsno0, product_arrs[1], c = 'r', marker = '+')
ax32.axhline(y=0, color="gray", lw = 1)
ax32.set_title("z_h binning")
ax32.set(xlabel = "z_h")
ax42.set(ylabel = "Affinity")
if(not single_pion):
    ax42.scatter(Mhbinsno0, product_arrs[4], c = 'r', marker = "+")
    ax42.axhline(y=0, color="gray", lw = 1)
    ax42.set_title("Mh binning")
    ax42.set(xlabel = "Mh (GeV)")
ax52.scatter(qTQbinsno0, product_arrs[2], c = 'r', marker = "+")
ax52.axhline(y=0, color="gray", lw = 1)
ax52.set_title("qTdivQ binning")
ax52.set(xlabel = "qTdivQ")
plot.savefig(dir_prefix+"Analysis/BoxAffinity/Plots_S24/May_14/Box_affinity_" + no_cut_file_name)
