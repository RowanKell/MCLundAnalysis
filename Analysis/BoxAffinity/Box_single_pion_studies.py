#Copied from Box.py, but dihadron treatments are removed for readability and simplicity

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

from datetime import date

logging.disable(logging.WARNING) 
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import os

def check_and_create_directory(directory_path):
    # Check if the directory exists
    if not os.path.exists(directory_path):
        print(f"Directory {directory_path} does not exist. Creating it...")
        # Create the directory
        os.makedirs(directory_path)
        print(f"Directory {directory_path} created successfully.")
    else:
        print(f"Directory {directory_path} already exists.")




def CalculateBoxAffinity(inRootFilePath, useDriver, plotFileName,multipleFiles,plot_title = ""):
    today = date.today()
    date_dir = today.strftime("%b_%d/")
    dir_prefix = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/"
    check_and_create_directory(dir_prefix + "Analysis/BoxAffinity/Plots_S24/" +date_dir)

    print("calculating single pion box affinity")

    treeName = "tree_driver" if useDriver else "tree_driver"
#     d = RDataFrame(treeName,dir_prefix + inRootFilePath + "*.root") #from LundAnalysis_single_pion.C
    useDriver = True
    d = RDataFrame(treeName,dir_prefix + "Analysis/BoxAffinity/root_files/July_9_old_R2_driver_6_files_low.root") #from LundAnalysis_single_pion.C
    # d = RDataFrame("tree_driver",dir_prefix + "../sidisregions/rowan_dev/root_files/driver_ratios_May23_500k.root") #from driver.py
    pTbins = np.linspace(0.1,0.8,8)
    xbins = np.array([0,0.1,0.13,0.16,0.19,0.235,0.3,0.5])
    zbins = np.array([0,0.35,0.43,0.49,0.55,0.62,0.7,0.83])
    qTQbins = np.array([0,0.1,0.3,0.5,0.8,1.1,1.5,2,2.5,3])

    pTbinsno0 = np.linspace(0.2,0.8,7)
    xbinsno0 = np.array([0.1,0.13,0.16,0.19,0.235,0.3,0.5])
    zbinsno0 = np.array([0.35,0.43,0.49,0.55,0.62,0.7,0.83])
    qTQbinsno0 = np.array([0.1,0.3,0.5,0.8,1.1,1.5,2,2.5,3])

    varName = np.array(["x", "z", "Q2", "pT", "R0", "R1", "R2"])

    
    def calculate_Affinity(d):

        px = [0 for i in range(7)]
        pz = [0 for i in range(7)]
        ppT=[0 for i in range(7)]
        pqTQ = [0 for i in range(9)]

        pxcut = [0 for i in range(7)]
        pzcut = [0 for i in range(7)]
        ppTcut=[0 for i in range(7)]
        pqTQcut=[0 for i in range(9)]
        if(useDriver):#if using driver, branch names are diff
            xformat = "x.x_t <= {} && x.x_t > {}"
            zformat = "z.z_t <= {} && z.z_t > {}"
            pTformat = "pT.pT_t <= {} && pT.pT_t > {}"
            R0format = "R0.R0_t < 0.3"
            R2format = "R2.R2_t < 0.3"
            R1format = "R1.R1_t < 0.3"
            qdivformat = "qTQ_calc.qTQ_calc_t <= {} && qTQ_calc.qTQ_calc_t > {}"
        else:
            xformat = "x <= {} && x > {}"
            zformat = "z <= {} && z > {}"
            pTformat = "pT <= {} && pT > {}"
            R0format = "R0 < 0.3"
            R2format = "R2_adjust < 0.3"
#             R2format = "R2 < 0.3"
            R1format = "R1 < 0.3"
            qdivformat = "qTQ_calc <= {} && qTQ_calc > {}"

        for i in range(7):
            pxcut[i] = d.Filter(xformat.format(xbins[i + 1],xbins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count()
            px[i] = d.Filter(xformat.format(xbins[i + 1],xbins[i])).Count()
            pzcut[i] = d.Filter(zformat.format(zbins[i + 1],zbins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count() 
            pz[i] = d.Filter(zformat.format(zbins[i + 1],zbins[i])).Count()
            ppTcut[i] = d.Filter(pTformat.format(pTbins[i + 1],pTbins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count()
            ppT[i] = d.Filter(pTformat.format(pTbins[i + 1],pTbins[i])).Count()
        for i in range(9):
            pqTQ[i] = d.Filter(qdivformat.format(qTQbins[i + 1],qTQbins[i])).Count()
            pqTQcut[i] = d.Filter(qdivformat.format(qTQbins[i + 1],qTQbins[i])).Filter(R2format).Filter(R1format).Filter(R0format).Count()

        for i in range(7):
            pxcut[i] = pxcut[i].GetValue()
            px[i] = px[i].GetValue()
            pzcut[i] = pzcut[i].GetValue()
            pz[i] = pz[i].GetValue()
            ppTcut[i] = ppTcut[i].GetValue()
            ppT[i] = ppT[i].GetValue()
        for i in range(9):
            pqTQcut[i] = pqTQcut[i].GetValue()
            pqTQ[i] = pqTQ[i].GetValue()


        pxval = [0 for i in range(7)]
        pzval = [0 for i in range(7)]
        ppTval=[0 for i in range(7)]
        pqTQval = [0 for i in range(9)]


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
        for i in range(9):
            if(pqTQ[i] == 0):
                pqTQval[i] = 0
            else: pqTQval[i] = pqTQcut[i] / pqTQ[i]
        Affinity_arrs = [pxval,pzval,pqTQval]
        return Affinity_arrs
    print(f"beginning single pion affinity calculation")
    product_arrs = calculate_Affinity(d)
    print(f"Finished calculating affinity arr for single pion")
    if(useDriver):
        marker_color = 'b'
        marker_shape = 'o'
    else:
        marker_color = 'r'
        marker_shape = '+'
    fig2, ((ax42, ax22), (ax32,ax52)) = plot.subplots(2, 2, figsize = (10, 10), dpi=50)
    if(plot_title == ""):
        fig2.suptitle("Single Pion BOX Affinity in the TMD region")
    else:
        fig2.suptitle(plot_title)
    ax22.scatter(xbinsno0, product_arrs[0], c = marker_color, marker = marker_shape)
    ax22.axhline(y=0, color="gray", lw = 1)
    ax22.set_title("x binning")
    ax22.set(xlabel = "x")
    ax32.scatter(zbinsno0, product_arrs[1], c = marker_color, marker = marker_shape)
    ax32.axhline(y=0, color="gray", lw = 1)
    ax32.set_title("z_h binning")
    ax32.set(xlabel = "z_h")
    ax42.set(ylabel = "Affinity")
    ax52.scatter(qTQbinsno0, product_arrs[2], c = marker_color, marker = marker_shape)
    ax52.axhline(y=0, color="gray", lw = 1)
    ax52.set_title("qTdivQ binning")
    ax52.set(xlabel = "qTdivQ")
    plot.savefig(dir_prefix+"Analysis/BoxAffinity/Plots_S24/" + date_dir + plotFileName)

# Binned = False #set true to use pre-binned kinematcs; false to use tree_MC and bin after with Box.py
# xlsxFileName = "xlsx/May23_test_500k"
# inRootFileName = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/May_23/for_box_500k.root"
# outRootFileName = "/root_files/May23_test.root"

# plotFileName = "driver_test_500k.pdf"
# CalculateBoxAffinity("Analysis/BoxAffinity" + outRootFileName, True, plotFileName)