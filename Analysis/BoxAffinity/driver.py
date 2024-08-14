from __future__ import division
import numpy as np
import ratlib as rat
import pandas as pd
import params as par
from tools import save, load, lprint, checkdir
import time
import ROOT as root
from array import array
import os
import uproot as up
from tqdm import tqdm

# xlsxFileName = "xlsx/May23_test"
# inRootFileName = "/w/hallb-scshelf2102/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/May_23/for_box_1k.root"
# outRootFileName = "root_files/May23_test.root"

def get_affinity(params,R0max,R1max,R2max,R3max,R4max,R5max,R1pmax,R1min,R2min,R3min,R1pmin):
    
    M         = params['M'] 
    M_h       = params['M_h']
    x         = params['x_bj']
    z         = params['z_h']
    Q         = params['Q']
    qT        = params['q_t']
    xi        = params['xi']
    zeta      = params['zeta']
    dkT       = params['delta_k_t']
    kit       = params['k_i_t']
    ki        = params['M_ki']
    kf        = params['M_kf']
    phi_i     = params['phi_i']
    phi       = params['phi']
    phi_ki       = params['phi_ki']

    R0 = np.maximum(
            np.abs(rat.get_R01( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki)),
            np.abs(rat.get_R02( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki)),
            np.abs(rat.get_R03( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki))
            )

    R1 = np.abs(rat.get_R1( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki))
#     R2 = np.ones(len(R1))*(qT**2) / (Q**2)
#     print(f"R1: {R1} | R2: {R2} | type(R1): {type(R1)} | type(R2): {type(R2)}")
    R2 = np.abs(rat.get_R2( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki))

    R3 = np.abs(rat.get_R3( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki))
   
    R4 = np.maximum(
            np.abs(rat.get_R41( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki)),
            np.abs(rat.get_R42( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki)),
            np.abs(rat.get_R43( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki))
            )   

    R1p = np.abs(rat.get_R1p( M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki))
    
    
    R5 = np.abs(np.log(R2))

    
    size=len(R1) # Number of iterations for one datapoint (kinematic bin)

    partonic_affinity = 0

    current_affinity = 0

    tmd_affinity = 0
    tmd_np_affinity = 0

    collinear_affinity = 0
    collinear_loworder_affinity = 0
    collinear_highorder_affinity = 0
    matching_affinity = 0
    target_affinity = 0
    soft_affinity = 0
    unclassified_affinity = 0 
    
   
    for i in range(size):
        # For partonic interpretation we have R0 << 1
        partonic = R0[i]<R0max
        if partonic : partonic_affinity += 1./size

        # For current region interpretation we have R1 << 1
        current = R1[i]<R1max
        if current : current_affinity += 1./size
            
        # for TMD affinity we have partonic, current region, and R2 << 1 NOTICE we do not use R1p as it does not appear to be reliable for current region
        tmd = R2[i]<R2max 
        tmdregion = partonic and current and tmd
        if tmdregion : tmd_affinity += 1./size
            
        # for non perturbative tmd R5 is small (resummation is not very important)    
        tmdnp = tmd and R5[i]<R5max 
        if partonic and current and tmdnp : tmd_np_affinity += 1./size
   
        #collinear QCD = partonic, current region, and R2 >> 1 and R4 << 1
        collinear = R2[i]>R2min and R4[i]<R4max
        collinearregion = partonic and current and collinear
        if collinearregion : collinear_affinity += 1./size

        #collinear QCD low order = collinear and R3 << 1 (2->2 process)
        collinearLowOrder = collinear and R3[i]<R3max 
        if partonic and current and collinearLowOrder : collinear_loworder_affinity += 1./size

        #collinear QCD high order = collinear and R3 >> 1 (2->3 process)
        collinearHighOrder = collinear and R3[i]>R3min 
        if partonic and current and collinearHighOrder : collinear_highorder_affinity += 1./size

        #for target region affinity we have partonic, not current region, and R1p << 1
        notcurrent = R1[i]>R1min
        target = R1p < R1pmax
        targetregion = partonic and notcurrent and target
        if targetregion : target_affinity += 1./size

        # The central (or soft) region is characterized by the production of  hadrons that are not the products of hard scattering but do not associate in any obvious way to a quark or target direction either.
        softregion = partonic and not current and not target and not collinear
        if softregion : soft_affinity += 1./size
                       
        # Matching is the region in between of TMD (R2<<1) and Collinear (R2>>1) regions
        matching = R2[i] > R2max and R2[i] < R2min
        if partonic and current and matching : matching_affinity += 1./size

        unclassifiedregion = not tmdregion and not collinearregion and not softregion and not targetregion
        if unclassifiedregion : unclassified_affinity += 1./size            
 
    # yi and yf depend on partonic kinematics    
    yi = rat.get_yi(M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki)
    yf = rat.get_yf(M,M_h,x,z,Q,qT,xi,zeta,dkT,kit,ki,kf,phi_i,phi,phi_ki)
            
    return partonic_affinity,current_affinity,tmd_affinity,tmd_np_affinity,collinear_affinity,collinear_loworder_affinity,collinear_highorder_affinity,matching_affinity,soft_affinity,target_affinity,unclassified_affinity,R0,R1,R2,R3,R4,R1p,yi,yf

def process_kinematics(fname,R0max,R1max,R2max,R3max,R4max,R5max,R1pmax,R1min,R2min,R3min,R1pmin,size,highAff,useMCNP):
    """
    Compute affinity accross kinematics of fname file 
    """
    print("We open file ",fname) # Write what file we open from directory /expdata
    kinematic_list_short = ["pT","Q2","x","z","R2_adjust"]
    np_list = ["M_ki","M_kf","delta_k_T","ki_T"]

    file_dir = fname
    num_files = len([name for name in os.listdir(file_dir) if not os.path.isdir(name)])
    file_names = [name for name in os.listdir(file_dir) if not os.path.isdir(name)]
    file_names_low = ["file_0.root","file_1.root","file_2.root","file_3.root","file_4.root","file_5.root","file_6.root"]
    if(highAff):
        tree_ext = ":tree_high"
    else:
        tree_ext = ":tree_low"
    for i in range(num_files):
        if(i == 0):
#             print(f"file_name: {file_names[i]}")
            uproot_df = up.open(file_dir + file_names[i] + tree_ext)
            tab = uproot_df.arrays(kinematic_list_short,library="pd")
#             np_MC = uproot_df.arrays(np_list,library="pd")
        else:
            if(not highAff):
                if(i > 6):
                    break
            tab = pd.concat([tab,up.open(file_dir + file_names_low[i] + tree_ext).arrays(kinematic_list_short,library="pd")],ignore_index = True)
#             np_MC = pd.concat([np_MC,up.open(file_dir + file_names[i] + tree_ext).arrays(np_list,library="pd")],ignore_index = True)
            print(f"file name: {file_names_low[i]}")
#     tab = pd.read_excel(fname)
#     tab=tab.to_dict(orient='list')
    if(highAff):
        np_MC = pd.read_excel("xlsx/np/August_13_high.xlsx")
    else:
        np_MC = pd.read_excel("xlsx/np/August_13_low.xlsx")
        print(f"len of np_MC: {len(np_MC)}")

#     npts=len(tab[list(tab.keys())[0]]) #npts is the number of datapoints in the xlsx file
    npts = min(len(np_MC),len(tab[list(tab.keys())[0]]))
    
    #Rowan edits
    print(f"npts: {npts}")
    ratios = np.zeros((3,npts,size))
    tab['partonicaff']=0.0
    tab['currentaff']=0.0


    tab['tmdaff']=0.0
    tab['tmdnpaff']=0.0

    tab['collinearaff']=0.0
    tab['collinearloworderaff']=0.0
    tab['collinearhighorderaff']=0.0
    tab['matchaff']=0.0
    tab['softaff']=0.0
    tab['targetaff']=0.0
    tab['unclassifiedaff']=0.0

    tab['R0'] = 0.0
    tab['R1'] = 0.0
    tab['R1p'] = 0.0
    tab['R2'] = 0.0
    tab['R3'] = 0.0
    tab['R4'] = 0.0
    tab['R5'] = 0.0
     
 
    # Let us modify the data file and add qT, xN, zN:   
    tab['qT'] = 0.0
    tab['xN'] = 0.0
    tab['zN'] = 0.0
    tab['yp'] = 0.0
    tab['yh'] = 0.0
    tab['yhtarget'] = 0.0
    tab['yi'] = 0.0
    tab['yf'] = 0.0
    
    #ROWAN EDIT
    tab['qTQ'] = 0.0




    print("starting first loop")
#     for i in tqdm(range(npts)):
#         tab.loc[i,'qT'] = tab['pT'][i]/tab['z'][i] # We will modify it later in this cell line 129
#         tab.loc[i,'xN'] = tab['x'][i] # We will modify it later
#         tab.loc[i,'zN'] = tab['z'][i] # We will modify it later
#         tab.loc[i,'yp'] = tab['x'][i] # We will modify it later
#         tab.loc[i,'yh'] = tab['x'][i] # We will modify it later
#         tab.loc[i,'yhtarget'] = tab['x'][i] # We will modify it later (this is for target region)
#         #Rowan EDIT - seems they are initializing with garbage
#         tab.loc[i,'qTQ'] = 0





    # Rapidities yi, yf depend on external partonic variables


    print("starting second loop")
    for i in tqdm(range(npts)):
        lprint('%d/%d'%(i,npts))        
        x   = tab['x'][i]
        z   = tab['z'][i]
        if 'pT' not in tab.keys(): 
            pT  = tab['pT2'][i]**0.5
        if 'pT' in tab.keys(): 
            pT  = tab['pT'][i]           
        Q2  = tab['Q2'][i]
        
        # Rowan Edit: only working with pi+/pi- and proton
#         tar = tab['target'][i]
#         had = tab['hadron'][i]
        tar = "pi+"
        had = "proton"

        params={}
        
        params['M'] = par.M
        params['M_h'] = par.Mpip

        params['x_bj'] = x
        params['z_h']  = z
        

        
        #Our functions take pT/z as an argument TODO check it
        params['q_t']  = pT/z    
            
        params['Q']    = Q2**0.5

        if tar.startswith('p'):
            params['M']    = par.M  
        elif tar.startswith('n'):
            params['M']    = par.MN
        elif tar.startswith('d'):
            params['M']    = par.MD
            
        if had.startswith('pi'):  
            params['M_h']=par.Mpi  
            if had.endswith("+"):
                params['M_h']=par.Mpip
            elif had.endswith("-"):
                params['M_h']=par.Mpim
            elif had.endswith("0"):
                params['M_h']=par.Mpi0   
        if had.startswith('k'):  
            params['M_h']=par.Mk    
            if had.endswith("+"):
                params['M_h']=par.Mkp
            elif had.endswith("-"):
                params['M_h']=par.Mkm
            elif had.endswith("0"):
                params['M_h']=par.Mk0 
                
        if had.startswith('h'):  
            params['M_h']=par.Mpi
            if had.endswith("+"):
                params['M_h']=par.Mpip
            elif had.endswith("-"):
                params['M_h']=par.Mpim
            elif had.endswith("0"):
                params['M_h']=par.Mpi0

        # Let us modify the data file and add qT:
#        if 'qT' not in tab.keys():            
        M         = params['M'] 
        M_h       = params['M_h']
        x         = params['x_bj']
        z         = params['z_h']
        Q         = params['Q']
        qT        = params['q_t']
        xi        = x # Not important for qT
        zeta      = z # Not important for qT
        dkT       = 0
        kit       = 0
        ki        = 0
        kf        = 0
        phi_i     = 0
        phi       = 0
        phi_ki       = 0
        #xN,zN, and qT do not depend on partonic variables, let us add them here
        tab.loc[i,'qT'] = np.sqrt( rat.get_qT2( M,M_h,x,z,Q,qT,xi,zeta,dkT,ki,kit,kf,phi_i,phi,phi_ki) ) 
        tab.loc[i,'xN'] = rat.get_xN( M,M_h,x,z,Q,qT,xi,zeta,dkT,ki,kit,kf,phi_i,phi,phi_ki)  
        tab.loc[i,'zN'] = rat.get_zN( M,M_h,x,z,Q,qT,xi,zeta,dkT,ki,kit,kf,phi_i,phi,phi_ki)  
        tab.loc[i,'yp'] = rat.get_yp( M,M_h,x,z,Q,qT,xi,zeta,dkT,ki,kit,kf,phi_i,phi,phi_ki)  
        tab.loc[i,'yh'] = rat.get_yh( M,M_h,x,z,Q,qT,xi,zeta,dkT,ki,kit,kf,phi_i,phi,phi_ki)  
        tab.loc[i,'yhtarget'] = rat.get_yh_target( M,M_h,x,z,Q,qT,xi,zeta,dkT,ki,kit,kf,phi_i,phi,phi_ki)  
        tab.loc[i,'qTQ'] = tab['qT'][i] / Q



        if(useMCNP):
            params['delta_k_t'] = np.array([np_MC['delta_k_T'][i]])
            params['k_i_t']     = np.array([np_MC['ki_T'][i]])
            params['M_ki']      = np.array([np_MC['M_ki'][i]])
            params['M_kf']      = np.array([np_MC['M_kf'][i]])
            params['zeta']      = np.array([np_MC['zeta'][i]])
            params['xi']      = np.array([np_MC['xi'][i]])
            params['phi_ki']      = np.array([np_MC['theta_ki'][i]])
            params['phi_i']      = np.array([np_MC['theta_H'][i]])
            params['phi']      = np.array([np_MC['theta_deltak'][i]])
#         print(f"type of np_MC: {type(np.array([np_MC['delta_k_T'].values[0]]))}")
#         print(f"value of np_MC: {np.array([np_MC['delta_k_T'].values[0]])}")
               
        
        #gen_np_values(params,x,z,0,size=size)  # test here
        gen_np_values(params,x,z,2,useMCNP,size=size) # Main level=2
        partonic_affinity,current_affinity,tmd_affinity,tmd_np_affinity,collinear_affinity,collinear_loworder_affinity,collinear_highorder_affinity,matching_affinity,soft_affinity,target_affinity,unclassified_affinity,R0,R1,R2,R3,R4,R1p,yi,yf = get_affinity(params,R0max,R1max,R2max,R3max,R4max,R5max,R1pmax,R1min,R2min,R3min,R1pmin)

        tab.loc[i,'yi'] = np.mean(yi)
        tab.loc[i,'yf'] = np.mean(yf)

        ratios[0,i] = R0
        ratios[1,i] = R1
        ratios[2,i] = R2
        
        if i == 0:
            R0MAX = max(R0)
            R0MIN = min(R0)
            R1MAX = max(R1)
            R1MIN = min(R1)
            R2MAX = max(R2)
            R2MIN = min(R2)
            R3MAX = max(R3)
            R3MIN = min(R3)
            R4MAX = max(R4)
            R4MIN = min(R4)
            R1PMAX = R1p
            R1PMIN = R1p
        else:
            if max(R0)> R0MAX: R0MAX = max(R0)
            if min(R0)< R0MIN: R0MIN = min(R0)
            if max(R1)> R1MAX: R1MAX = max(R1)
            if min(R1)< R1MIN: R1MIN = min(R1)
            if max(R2)> R2MAX: R2MAX = max(R2)
            if min(R2)< R2MIN: R2MIN = min(R2)
            if max(R3)> R3MAX: R3MAX = max(R3)
            if min(R3)< R3MIN: R3MIN = min(R3)
            if max(R4)> R4MAX: R4MAX = max(R4)
            if min(R4)< R4MIN: R4MIN = min(R4)
            if R1p> R1PMAX: R1PMAX = R1p
            if R1p< R1PMIN: R1PMIN = R1p

 # Let us also save values for Rs here, we will save mean values of Rs obtained in Monte Carlo:
        tab.loc[i,'R0'] = np.mean(R0)
        tab.loc[i,'R1'] = np.mean(R1)
        tab.loc[i,'R1p'] = np.mean(R1p)
        tab.loc[i,'R2'] = np.mean(R2)
        tab.loc[i,'R3'] = np.mean(R3)
        tab.loc[i,'R4'] = np.mean(R4)
        tab.loc[i,'R5'] = np.abs(np.log(np.mean(R4)))

        tab.loc[i,'partonicaff'] = partonic_affinity
        tab.loc[i,'currentaff'] = current_affinity     
        tab.loc[i,'tmdaff'] = tmd_affinity
        tab.loc[i,'tmdnpaff'] = tmd_np_affinity
        tab.loc[i,'collinearaff'] = collinear_affinity
        tab.loc[i,'collinearloworderaff'] = collinear_loworder_affinity
        tab.loc[i,'collinearhighorderaff'] = collinear_highorder_affinity
        tab.loc[i,'matchaff'] = matching_affinity
        tab.loc[i,'softaff'] = soft_affinity
        tab.loc[i,'targetaff'] = target_affinity
        tab.loc[i,'unclassifiedaff'] = unclassified_affinity
        
        #ROWAN EDIT
        #Add qT/Q
#       np_list = ["M_ki","M_kf","delta_k_T","ki_T"]
        for np_idx in range(len(np_list)):
            tab.loc[i,np_list[np_idx]] = np_MC[np_list[np_idx]][i]
        tab.loc[i,"zeta"] = params["zeta"]
        tab.loc[i,"xi"] = params["xi"]
        tab.loc[i,"phi_ki"] = params['phi_ki']
        tab.loc[i,"phi_i"] = params['phi_i'] 
        tab.loc[i,"phi"] = params['phi']   
      
    print("\n R0max = %s, R0min = %s"%(R0MAX, R0MIN))
    print("\n R1max = %s, R1min = %s"%(R1MAX, R1MIN))
    print("\n R2max = %s, R2min = %s"%(R2MAX, R2MIN))
    print("\n R3max = %s, R3min = %s"%(R3MAX, R3MIN))
    print("\n R4max = %s, R4min = %s"%(R4MAX, R4MIN))
    print("\n R1pmax = %s, R1pmin = %s"%(R1PMAX, R1PMIN))

    tab=pd.DataFrame(tab)
    return tab, ratios

def write_affinity_to_excel(tab,output):
    print("We write output into a file ",output) # Write what file we write output into directory /data
    tab.to_excel(output)      
def gen_np_values(params,x,z,level,useMCNP,size=100):
    #params['xi']   = np.random.uniform(x,x+0.1,size) # Should we generate partonic variable in a narrow interval?
    #params['zeta'] = np.random.uniform(z,z+0.1,size)

    #params['xi']   = np.random.uniform(x,1,size) # Should we generate partonic variable in a wide interval?
    #params['zeta'] = np.random.uniform(z,1,size)    
    
    delta = 0.1
    
    deltax = x+delta
    if deltax >= 1:
        deltax = 1
        
    deltaz = z+delta
    if deltaz >= 1:
        deltaz = 1

    if   level== 0: lower,upper=0.,0.225
    elif level== 1: lower,upper=0.0,0.25
    elif level== 2: lower,upper=0.0,par.M 
    #elif level== 2: lower,upper=0.25,0.75
    elif level== 3: lower,upper=0.01,0.05

#     print(f"type of random: {type(np.abs(np.random.normal((upper+lower)/2,(upper-lower)/2,size)))}")
#     print(f"value of random: {np.abs(np.random.normal((upper+lower)/2,(upper-lower)/2,size))}")
        
    if(not useMCNP):
        params['delta_k_t'] =  np.abs(np.random.normal((upper+lower)/2,(upper-lower)/2,size))
        params['k_i_t']     =  np.abs(np.random.normal((upper+lower)/2,(upper-lower)/2,size))
        params['M_ki']      =  np.abs(np.random.normal((upper+lower)/2,(upper-lower)/2,size)) 
        params['M_kf']      =  np.abs(np.random.normal((upper+lower)/2,(upper-lower)/2,size))
        params['xi']   = np.random.uniform(x,deltax,size) # Should we generate partonic variable in a wide interval?
        params['zeta'] = np.random.uniform(z,deltaz,size)
        lower_phi = 0
        upper_phi = 2*np.pi
        params['phi']       = np.random.uniform(lower_phi,upper_phi,size)
        params['phi_i']     = np.random.uniform(lower_phi,upper_phi,size)
        params['phi_ki']     = np.random.uniform(lower_phi,upper_phi,size)
    
def main00(rootFileName,highAff,useMCNP):

    fnames = {rootFileName}
    # The "box" size ~ qT/Q < 0.3 
    R0max = 0.3

    R1max = 0.3
    R1min = R1max*3
    
    R2max = 0.3
    R2min = R2max*3
    
    R3max = 0.3
    R3min = R3max*3
    
    R4max = 0.3
    R5max = 0.3
 
    R1pmax= 0.3
    R1pmin= R1pmax*3

    
    #--compute affinity
    size=1
    tabs = []
    for fname in fnames:
        tab,ratios=process_kinematics(fname,R0max,R1max,R2max,R3max,R4max,R5max,R1pmax,R1min,R2min,R3min,R1pmin,size,highAff,useMCNP)
#         output = fname + "_w_affinity.xlsx"
        #write_affinity_to_excel(tab,output)
        tabs.append(tab)
    return tabs, ratios

def maxmin(rootFilePath, ratios):
    file = root.TFile.Open(rootFilePath, "UPDATE")
    tree_max = root.TTree("tree_max","tree_max")

    R0_hist_max = array('d',[0])
    R1_hist_max = array('d',[0])
    R2_hist_max = array('d',[0])

    R0_hist_min = array('d',[0])
    R1_hist_min = array('d',[0])
    R2_hist_min = array('d',[0])

    tree_max.Branch('R0_hist_max', R0_hist_max,'R0_hist_max/D')
    tree_max.Branch('R1_hist_max', R1_hist_max,'R1_hist_max/D')
    tree_max.Branch('R2_hist_max', R2_hist_max,'R2_hist_max/D')

    tree_max.Branch('R0_hist_min', R0_hist_min,'R0_hist_min/D')
    tree_max.Branch('R1_hist_min', R1_hist_min,'R1_hist_min/D')
    tree_max.Branch('R2_hist_min', R2_hist_min,'R2_hist_min/D')

    R0_hist_max[0] = max(ratios[0])
    R1_hist_max[0] = max(ratios[1])
    R2_hist_max[0] = max(ratios[2])

    R0_hist_min[0] = min(ratios[0])
    R1_hist_min[0] = min(ratios[1])
    R2_hist_min[0] = min(ratios[2])
    tree_max.Fill()
    tree_max.Write()
