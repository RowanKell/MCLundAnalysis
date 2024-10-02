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
    
    treeName = "tree_driver" if useDriver else "tree_MC"
    # print(f"\n\npath: {inRootFilePath}\n\n")
    if os.path.isdir(inRootFilePath):
        d = RDataFrame(treeName,inRootFilePath + "*.root") #from LundAnalysis_single_pion.C
        
    elif os.path.isfile(inRootFilePath):
        d = RDataFrame(treeName,inRootFilePath) #from LundAnalysis_single_pion.C
    else:
        print("error: not a file or directory")
        exit()

    #closure test
    # up_file = uproot.open(inRootFilePath + ":" + treeName).arrays(library = "pandas")
    # R0 = np.array(up_file['z'])
    # R1 = np.array(up_file['x'])
    # R2 = np.array(up_file['qTQ_hadron'])
    # fig_test, ax_test = plot.subplots(1,3,figsize = (10,3))
    # ax_test[0].hist(R0,bins = 100);
    # ax_test[1].hist(R1,bins = 100);
    # ax_test[2].hist(R2,bins = 100);
    # fig_test.tight_layout()
    # if(useDriver):
    #     fig_test.savefig("Plots_S24/Sep_18/closure_driver_kinematics.pdf")
    # else:
    #     fig_test.savefig("Plots_S24/Sep_18/closure_MC_kinematics.pdf")

    
    pTbins = np.linspace(0.1,0.8,8)
    xbins = np.array([0,0.1,0.13,0.16,0.19,0.235,0.3,0.5])
    zbins = np.array([0,0.35,0.43,0.49,0.55,0.62,0.7,0.83])
    qTQbins = np.array([0,0.1,0.3,0.5,0.8,1.1,1.5,2,2.5,3])

    pTbinsno0 = np.linspace(0.2,0.8,7)
    xbinsno0 = np.array([0.1,0.13,0.16,0.19,0.235,0.3,0.5])
    zbinsno0 = np.array([0.35,0.43,0.49,0.55,0.62,0.7,0.83])
    qTQbinsno0 = np.array([0.1,0.3,0.5,0.8,1.1,1.5,2,2.5,3])

    varName = np.array(["x", "z", "Q2", "pT_BF", "R0", "R1", "R2"])

    
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
            pTformat = "pT_BF.pT_BF_t <= {} && pT_BF.pT_BF_t > {}"
            R0format = "R0.R0_t < 0.3"
            R2format = "R2.R2_t < 0.3"
            R1format = "R1.R1_t < 0.3"
            qdivformat = "qTQ_hadron.qTQ_hadron_t <= {} && qTQ_hadron.qTQ_hadron_t > {}"
        else:
            xformat = "x <= {} && x > {}"
            zformat = "z <= {} && z > {}"
            pTformat = "pT_BF <= {} && pT_BF > {}"
            R0format = "R0 < 0.3"
            # R2format = "R2_adjust < 0.3"
            R2format = "R2 < 0.3"
            R1format = "R1 < 0.3"
            qdivformat = "qTQ_HF <= {} && qTQ_HF > {}"

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
        marker_shape = 'D'
        marker_size = 40
    else:
        marker_color = 'r'
        marker_shape = '*'
        marker_size = 70
    fig2, axs2 = plot.subplots(1, 3, figsize = (10, 3.5))
#     if(plot_title == ""):
#         fig2.suptitle("Single Pion BOX Affinity in the TMD region")
#     else:
#         fig2.suptitle(plot_title)
    axs2[0].scatter(xbinsno0, product_arrs[0],marker_size, c = marker_color, marker = marker_shape)
    axs2[0].axhline(y=0, color="gray", lw = 1)
#     axs2[0].set_title("x binning")
    axs2[0].set_xlabel("$x_{Bj}$", fontsize=14)
    axs2[0].set_ylabel("Affinity", fontsize=14)
    axs2[1].scatter(zbinsno0, product_arrs[1],marker_size, c = marker_color, marker = marker_shape)
    axs2[1].axhline(y=0, color="gray", lw = 1)
#     axs2[1].set_title("z_h binning")
    axs2[1].set_xlabel("$z_h$", fontsize=14)
    axs2[2].scatter(qTQbinsno0, product_arrs[2],marker_size, c = marker_color, marker = marker_shape)
    axs2[2].axhline(y=0, color="gray", lw = 1)
#     axs2[2].set_title("qTdivQ binning")
    axs2[2].set_xlabel("$q_T/Q$", fontsize=14)
    fig2.tight_layout()
    plot.savefig(dir_prefix+"Analysis/BoxAffinity/Plots_S24/" + date_dir + plotFileName)

# Binned = False #set true to use pre-binned kinematcs; false to use tree_MC and bin after with Box.py
# xlsxFileName = "xlsx/May23_test_500k"
# inRootFileName = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/May_23/for_box_500k.root"
# outRootFileName = "/root_files/May23_test.root"

# plotFileName = "driver_test_500k.pdf"
# CalculateBoxAffinity("Analysis/BoxAffinity" + outRootFileName, True, plotFileName)