#imports
from ROOT import TFile,TTree,TCanvas,TH1F,gStyle,TLatex,gPad,TLegend,TLorentzVector,TH2F,TLine,TF1,TBox,RDataFrame,TPad,TF2
import ROOT
import numpy as np
import awkward as awk
import uproot
import pandas as pd
import matplotlib.pyplot as plot

#don't worry about this
logging.disable(logging.WARNING) 
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

#set the names of the files and trees
filepath = "PATH/TO/FILE"
treename = "NAME_OF_TREE"
inFile = ROOT.TFile.Open(filepath,"READ")

# TODO: Define the kinematic variables in functions like this
# (Replace param1, param2 with the necessary components)
def Q2(param1, param2):
    # TODO: put calculation here
    return Q2

#grab trees
tree = inFile.Get(treename)
#initialize array for kinematic
Q2array = []
#loop over all events in tree
for entry_num in range(0, tree.GetEntries()):
    #select the current entry
    tree.GetEntry(entry_num)
    #use the getattr() function to get the value of a branch at the current entry
    selectedHadPx = getattr(tree, "SelectedHadPx")
    #calculate kinematics
    Q2 = Q2(param1, param2)
    #add kinematics to the arrays
    Q2array.append(Q2)