# JLab CheatSheet

#### Rowan Kelleher

##### May 31 2022

## Login

#### Credentials:

Username: rojokell

Password: D%fd2tt57&CmJ5%1

#### Login Steps:

1. Type "ssh rojokell@login.jlab.org" into ubuntu terminal
2. Enter password
3. Then type "ssh rojokell@ifarm"
4. Enter password
5. type "cd workdir"
6. You are now in your work directory

## Root

#### Running Root

1. Type "clas12root"

## Jupyter-Notebook

#### Opening a notebook

1. Go to "https://jupyterhub.jlab.org/hub/spawn" in a browser
2. Start a server with CLAS12
3. To enter your work directory, use a symbolic link
4. Example: %cd workdir

#### Using Bash commands

1. Use "!" in front of bash commands while in Jupyter
2. use "%cd" for changing directory or just "cd"

#### Convert Tree into Array

1. imports: 

```python
import awkward as ak
import numpy as np
import uproot
import pandas as pd
from ROOT import TFile, TTree, TCanvas, TH1F, gStyle, TLatex, gPad, TLegend, TLorentzVector, TH2F, TLine, TF1, TBox, RDataFrame, TPad, TF2
import ROOT
```

2. Open the file with uproot and assign a variable to it:

```python
<uproot object name> = uproot.open("path/to/.root/file")
#Example:
up4_file = uproot.open("CLAS12Analysis/tutorials/output.root")
```

3. Create object from tree:

```python
<tree object name> = <uproot object name>["<tree name>"]
#Example:
up4_reco = up4_file["tree_reco"]
```

4. You can show the contents of the tree using show()

```python
#Example:
up4_reco.show()
```

5. Create array object:

```python
<array name> = <tree object name>["<kinematic variable>"].array(library="np")
#Example:
xarray = up4_reco["x"].array(library="np")

```

6. Reference: https://awkward-array.org/how-to-convert-uproot.html

#### Convert Array to excel file

1. Set file path:

```python
#Example: creating .xlsx file in your current directory
filepath = 'test.xlsx'
#Example: creating .xlsx file in a directory
filepath = 'mydirectory/test.xlsx'
```

2. Create DataFrame variable:

```python
<dataframe variable> = pd.DaraFrame(<array>)
#Example:
<df> = pd.DataFrame(xarray)
```

3. Create excel file

```python
<dataframe variable>.to_excel(<filepath>, index=False)
```

#### Creating Histogram from Tree

1. Create RDataFrame variable from tree:

```python
<file object name> = RDataFrame("<tree name>" , "path/to/.root/file")
#Example:
file = RDataFrame("tree_MC","CLAS12Analysis/tutorials/output.root") 
```

2. Create histogram object:

```python
<histogram object> = <file object name>.Histo1D(("<histogram name>","<histogram title", <number of bins>, <x axis min>, <x axis max>),"<variable>")
#Example: create histogram object of W with 100 bins, 0<W<0.6
h = file.Histo1D(("h","W Histogram", 100 ,0,0.6),"W")
```

3. Create Canvas:

```python
<canvas object> = TCanvas("<canvas name>","<canvas title>",<size option 1>,<size option 2>)
#Example:
c = TCanvas("c","c",1100,800)
```

4. Draw the histogram and canvas

```
h.Draw("hist")
c.Draw()
```

## Ifarm stuff

#### Hipo directory

Hipo files for CLAS12 are stored at:

```python
/cache/clas12/rg-a/production/recon/fall2018/torus*/pass1/v1/dst/train/nSidis/
```

*both directories "torus+1" and "torus-1" within fall2018 are used

#### Output directory for processed root files

Root files should be stored in:

```
/work/clas12/users/rojokell/CLAS12Analysis/data/*
```

*the end directory/file depends on the directory we pull from (torus+1/-1 and fall2018 vs 2019)

## Macros

#### rowanTestpipluspiminus.C

1. Copied from "pipi0_process.C"
2. Changed "gmat" to "rojokell" in lines 5,6,14,16
3. Changed final state settings on line 43

```python
#Original:
settings.addFinalState(22,2,false); // 2 or more gammas
#New:
settings.addFinalState(-211,1,true); //Exactly 1 pi-
```

4. Replaced minimum energy of gammas with minimum momentum of pi- to match pi+

```python
#Original:
settings.addPIDforEmin(22,0.6);     // Gammas must have minimum energy of 0.6 GeV
#New:
settings.addPIDforPmin(-211,1.25);  // Pi- must have minimum momentum of 1.25 Gev
```

5. Added maximum chisquared for Pi- to match pi+

6. Changed post processing method from pipi0 to pi+pi-

7. ```
   #Original
     settings.setPostProcessMethod("pipluspi0");
   #New
     settings.setPostProcessMethod("pipluspiminus");
   ```

8. **Seems to be good from here, just need to finish adjusting PostProcess.C and PostProcess.h**



#### Postprocess.C

1. Added post process id for pi+pi-

2. ```c++
     case PROCESS_ID::pipluspiminus:
       pipluspiminus(_tree_postprocess);
       break;
   ```

3. Duplicated pip0 PostProcess function and replaced the name with "pipluspiminus"

4. Commented out "double Mgg;" as photons are not important here (ine 338)

5. Commented out "_tree_postprocess->Branch("Mdiphoton",&Mgg)" (line 349)

6. Changed TLorentzVectors from "pi" and "pi0" to "piplus" "piminus" (line 373-4)

7. Commented out gamma vectors (lines 375-6)

8. Changed "pi" to "piplus" and "pi0" to "pi0" (line 380-1, 383-4, 387-8)

9. Commented out "zpair" as it seems to be referring to a photon pair

10. Commented out "E1" and "E2" and "beta1" and "beta2" on lines 396-400

11. Removed Identification for loop for photons

12. Added Identification loop for pi- (lines 411-421)

13. Added z calculation for piplus and piminus (replacing pi) (lines 433-437)

    1. Kept same formula which matches affinity paper formulas

14. Added xF calculation replacing pi xF calc (447-451)

15. Replacing photon identification with cuts/appendings not associated with photons

16. commented out if statement on line 499 (and end bracket on 509) to check if this changes stuff

#### Troubleshooting Error

1. Error: postprocess method is unrecognized
2. Known: rowanTestpipluspiminus.C is not issue
3. Changed pipluspiminus to pipiminus in PostProcess.C (lines 336, 50,) and PostProcess.h (line 69)
   1. This avoids giving the same name (pipluspiminus) to different things

#### makeJobs.sh

1. Replaced "gmat" with "rojokell" on lines 3,8,9,11,38

2. Changed path to clas12 to fit my work directory
   1. deleted "packages" in line 38:
   
3. ```
   /w/hallb-scshelf2102/clas12/users/rojokell/packages/clas12root
   ```

4. Changed "rootname" from may24_ to june1_

5. Changed process directory to tutorials instead of macros/dihadron_process

   1. my test macro is in tutorials

6. Changed process codename to "rowanTestpipluspiminus.C"

   1. This is my test macro that may or may not be finding pi plus and minuses

#### Attempting to get new kinematics into trees

Kinematics wanted for histos:

1. x_F
   1. Feynman X
2. z_h
   1. Zpair
3. P_h
   1. transverse hadron momentum
4. M_x
   1. Missing Mass

PostProcess.C changes

1. Need to define the variables underneath the function beginning like on line 60
2. Need to set _tree_postprocess->Branch("<branch name>",&<variable name>)
3. Need to calculate each variable
   1. **x_F**
      1. $\frac{2P_h \cdot q}{|q|W}$
   2. **z_h**
      1. $\frac{P \cdot P_h}{P \cdot q}$
   3. P_h perp
      1. $P_h sin\theta_{\gamma h}$
      1. ^ unneeded, .Perp() function from TLorentzVector is much better
   4. M_x
      1. Calculated from THayward
4. Set zpiplus and zpiminus and P_piplus to 0 on line 400
5. Set the 4 variables to 0 on line 415
6. initialize theta_gh on line 399
7. Redefined xF and P_h from gregs suggestions
   1. x_F is now caclulated using
   2. 1. $\frac{2P_h \cdot q}{|q|W}$
      2. This requires boosting - see Timothy Hayward's code in https://github.com/tbhayward/SIDIS_dihadron/blob/master/Dihadrons.java
   3. P_h is now just
   4. $P_h = dihadron.Perp(q)$

8. June 8:
   1. x_F did not work when I did:
      1. $x_F = 2 * (dihadron.Vect()).Dot(q.Vect()) / (q.Mag() * W);$


#### NEED TO DO TODAY:

1. Try to fix phi_h and phi_R
2. Try to make 3x4 histogram pages
3. fill slideshow with hisograms
4. add cuts at the end
5. Learn about all cuts and kinematics
   1. W
      1. Invariant Mass
      2. Total mass in the COM rest frame
      3. In other frames, momentum is nonzero meaning that the total mass/relativistic mass is greater than the invariant mass, but invariant mass remains the same
      4. can be calculated in rest frame using
      5. $M^2 = E^2 - P^2$
   2. Q^2
      1. Measure of resolution power
   3. x_F
      1. https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.23.1415
         1. Feynman's paper where he introduced this variable
      2. Ratio of logitudanal momentum P_h to the total available W
   4. x
      1. Measure of momentum fraction of struck quark
      2. lorentz invariant
   5. y
      1. Measure of inelasticity
      2. fractional energy loss of incoming particle
   6. z
      1. fraction of virtual photon energy v carried by the hadron in the rest frame of proton beam
   7. M_h
      1. Relativistic mass of combination of both particles (pi+ and pi-)
      2. unlike invariant mass, relativistic mass has a contribution from the energy of the system depending on the reference frame
   8. M_x
      1. Missing mass: mass of particles not included in the final state particles
      2. In this case, any mass not found in the pi+pi- dihadron
      3. Example:
         1. such as an electron-and-box, if the electron bounces at high speed inside the box. It is only the lack of total momentum in the system (the system momenta sum to zero) which allows the kinetic energy of the electron to be "weighed". If the electron is *stopped* and weighed, or the scale were somehow sent after it, it would not be moving with respect to the scale, and again the relativistic and rest masses would be the same for the single electron
   9. P_ht
      1. Transverse momentum
   10. theta
       1. Polar angle between the hadron $P_1$ in the pair center of mass frame and the direction of the pair, $P_h$, in the photon-target rest frame (gN)
   11. phi_R
       1. $\phi_R = \mathbf{\frac{(q \times l) \cdot R_{\perp}}{|(q \times l) \cdot R_{\perp}|} }arccos\left (\mathbf{\frac{(q \times l) \cdot (q \times R_{\perp})}{|q \times l|\cdot| q \times R_{\perp}|}}  \right )$
          1. https://arxiv.org/pdf/1702.07317.pdf
       2. ![](C:\Users\rowan\Downloads\phi_Rphi_h.PNG)
       3. Phi R is the angle between the lepton scattering plane and the dihadron plane
          1. The dihadron plane contains both of the hadrons (pions) after collision
       4. 1. 
   12. Phi_h
       1. Phi h is the angle between the scattering plane and the q x P_h plane
          1. The q x P_h plane is a plane with the normal as the cross product of q and Ph, or the orthogonal vector to those two vectors
       2. $\phi_h = \mathbf{\frac{(q \times l) \cdot P_h}{|(q \times l) \cdot P_h|} }arccos\left (\mathbf{\frac{(q \times l) \cdot (q \times P_h)}{|q \times l|\cdot| q \times P_h|}}  \right )$


## Slurm

#### Commands

1. Check queue: 

   1. ```
      squeue -u rojokell
      ```

2. Kill open jobs

   ```
   scancel -u rojokell
   ```

3. To run a job

   ```
   sbatch <filename>.sh
   ```

   

4. Check status online: https://scicomp.jlab.org/scicomp/slurmJob/activeJob



## Affinity

#### Kinematics needed

1. $R_0$
   1. $max(|\frac{k_i^2}{Q^2},\frac{k_f^2}{Q^2},\frac{\delta k^2_T}{Q^2}|)$

1. $R_1$
   1. $\frac{P_h \cdot k_f}{P_h \cdot k_i}$
2. $R_2$
   1. $\frac{|k^2|}{Q^2}$
3. $z_h$
   1. $\frac{P \cdot P_h}{P \cdot q}$
4. $x$
   1. $\frac{Q^2}{2P \cdot q}$
   2. Bjorken scaling variable
5. $Q^2$
   1. Momentum transfer squared
   1. $2 * E * E' * (1.0 - cth)$
   1. E = initial lepton energy
   1. E' = final lepton energy
   5. cth -
      1. $\frac{Pz}{\sqrt{Pz * Pz + Pt * Pt}}$
      2. where $Pt = \sqrt{Px * Px + Py * Py}$
   6. https://www.ippp.dur.ac.uk/~krauss/Lectures/IntroToParticlePhysics/Lecture7.pdf
6. $p_T$
   1. Hadron momentum

#### Variables needed for calculation

1. P
   1. Momentum of initial Hadron
2. q
   1. Momentum Transfer of Incident Lepton
3. $P_h$
   1. Momentum of final hadron
4. $k$
   1. $k_f - q$
5. $k_f$
   1. Hadronizing parton(s) momentum
6. $k_i$
   1. Initial parton momentum
7. $\delta k_T^2$
   1. "Characterizes the size of the intrinsic transverse momentum of the parton"

## Cuts Made

#### PostProcess.C cuts

1. Vertex Position Cut
   1. |vz(electron) - vz(pion)| < 20 cm
      1. Hadron vertex has to be within 20cm of electron vertex

#### rowanTestpipluspiminus.C

1. Channel selection
   1. Q2 range: 1 < Q2 < 100
      1. Mass of virtual photon
   2. W range: 2 < W < 100
      1. Hadronic system mass
   3. Energy fraction y cut: 0 < y < 0.8
      1. Limit radiative effects
2. Minimum momentum
   1. Both pions: momentum > 1.25 GeV
      1. limit radiative effects and low momentum
3. Electron vertex
   1. -8 < vz < 3
   2. reject electrons scattered off of the target window
4. Chi2max
   1. pi+: abs(chi2) < 2.64
   2. pi-: abs(chi2) < 2.79
      1. These separate pions from kaons

#### Histogram Plotting

1. Mass cut
   1. Missing Mass > 1.5 GeV
   2. Timothy found that asymmetries were changing as a function of missing mass below 1.5 GeV, indicatng contributions from exclusive regions
2. x-feynman > 0 for each pion
   1. Reduce target-fragmentation region

#### Fiducial Cuts

1. Electron
   1. v > 9 (in code: v < 9 returns false)
   2. w > 9
      1. These are electron calorimeter coordinates
2. theta: $\frac{5\pi}{180}$ < theta < $\frac{35\pi}{180}$

## Work Notes

#### June 14th 2022

1. Looking at MC hipo file
   1. /cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo
   2. In position 2914 MC::Lund:
      1. pids:
         1. 11 - electron
         2. 2212 - proton
         3. 22 - photon
         4. 2 - up quark
         5. 92 ?
         6. 113 - rho 0
         7. 213 - rho +
         8. 211pi+
         9. 111 pi0
         10. 2103 - (ud)1
         11. 2101 - (ud)0
         12. 3122 - lambda (uds)
         13. 321 - K+
         14. 323 - K*(892)+ 
      2. Mapping 
         1. electron - index 1
            1. photon - index 3
            2. electron - index 4
         2. proton - index 2 - LT
            1. up quark - idx 6 -LT 1
               1. up quark - idx 8 - LT 1
                  1. 92 particle idx 10 - LT 0
                     1. rho0 - idx 11 - LT 0
                        1. pi+ - idx 14 - LT 1
                        2. pi- - idx 15 - LT -1
                     2. rho+ idx 12 - LT - 1
                        1. pi+ - idx 16 - LT 1
                        2. pi0 - idx 17 - LT 0
                           1. photon - idx 18 - LT 0
                           2. photon - idx 19 - LT 0
                     3. photon - idx 13 - LT 
            2. proton - idx 7 - LT 1
         3. anti-up quark - idx 5 - LT -1
            1. anti-up quark idx 9 - LT -1



#### June 15 2022

1. Hipo Bank descriptions
   1. https://clasweb.jlab.org/wiki/index.php/CLAS12_DSTs

#### June 16 2022

1. Trying to trace back whats going on in SIDISKinematicsReco.C
   1. SIDISKinematicsReco()
      1. Takes in outfile name
      2. Used in processing macro (such as rowanTestpipluspiminus.C at line 29)
         1. Takes output file name from either macro file, or from makeJobs.sh script which feeds root file names to macro
   2. Init()
      1. Makes TFile
      2. Makes Variable maps
         1. Event variable map
            1. holds nParticles, x, y, Q2, W, nu, helicity
            2. All variables that are specific to an entire event, not a particle
         2. Particle map
            1. holds pid, px, py, pz, pt, p, E, evtgen_E, theta, eta, phi, v, pindex, beta, chi2, parentID, parentPID
         3. These maps seem to hold what will become branches in the eventual TTree?
      3. Creates MC TTree with event and particle branches
         1. branches are looped through the particle variable map
            1. it = _map_event.begin(); it! _map_event.end()
               1. Uses "it->first" and "it->second" for branch names
      4. Creates reconstructed TTree as copy of monte carlo tree
         1. Keeps same event and particle branches
      5. Loading in HipoFiles
         1. Checks to make sure that there are hipo files fed to the process
            1. Uses _settings.hipoFileStrings().size() to check for files
         2. uses `_chain.Add(_settings.hipoFileStrings().at(idx).c)str())` to add the files with a for loop
      6. Sets beam energy
         1. Takes value from _settings.electronBeamEnergy() which is set in process macro
         2. saves value as _electron_beam_energy
      7. Creates/configs clas12reader
         1. `_config_c12=chain.GetC12Reader();`
      8. End of Init() stuff
   3. InitHipo()
      1. Here the function makes PID cuts for all events
         1. creates vector finalStatePIDs from `_settings.getFinalStatePIDS()`
         2. for loop loops over finalStatePIDs.size()
            1. In the process macro, _settings.addFinalState takes in final state PIDs and numbers
            2. In the Settings.C, addFinalState pushes the pids back to _fPID which is returned by finalStatePIDs if theres no repeats
         3. sets pid variable to finalStatePIDS.at(idx)
         4. does `_config_c12->addExactPid(pid,npid)`
            1. makes the PID cut which is just manually typed in analysis with extra bins macro
      2. Adds banks to _config_c12 object
         1. looks at _settings.getEventRecoMethod() for "useRecKinematicsBank"
            1. if this is the setting, then it sets _ix, _iQ2, _iy, _inu, _iW equal to the bank
               1. Ex: `_ix = _config_c12->getBankOrder(_idx_RECKin,"x");`
            2. Appears to do nothing if the setting is set to something else
         2. In my macro, the setting is 
            1. ` settings.setEventRecoMethod(Settings::eventRecoMethod::useLargestPinFD);` 
   4. process_events()
      1. establish clas12 event parser
         1. `auto &_c12= _chain.C12ref();`
      2. Create fiducial cuts object
         1. if the setting for fiducial cuts is true, then it makes a cut based on the _c12 object
      3. Creates map for reco particles
         1. recoparticleMap
      4. Creates map for true particles
         1. particleMap
      5. parsing through reco data
         1. looks for _settings.doReco() to be true (in processing macro)
         2. skips cut particles
            1. `if(CollectParticlesFromReco( _c12, recoparticleMap )!=0)`
      6. parsing through MC data
         1. looks for _settings.doMC() to be true
         2. adds particle info:
            1. `CollectParticlesFromTruth(_c12,particleMap);`
               1. **Try to figure out where CollectParticlesFromTruth comes from, maybe ask greg**
               1. Defined at end
         3. connects MC and reco particles
      7. Writing to TTrees
         1. checks for doReco setting
         2. resets the branch map ResetBranchMap
         3. Writes: `WriteParticlesToTree(recoparticleMap);`
         4. fill tree with _tree_Reco->Fill()
         5. Does same with MC
   5. CollectParticlesFromTruth()
      1. creates MC particle pointer
         1. `mcparticles=_c12->mcparts()`
      2. Loops over all particles to get particles in MC::Lund
         1. Makes sure all particles are final state (This would need to be changed/something would need to be added if I want to look for quarks that fragment into the dihadron)
      3. Makes a new SIDISParticle *sp
      4. sets all variables with mcparticles->getPid() for example
      5. uses sp pointer to set the properties to the SIDISParticle:: objects
      6. Adds SIDISParticle to map
   6. CollectParticlesFromReco()

#### June 17 2022

Trying to see the connections between SIDISKinematicsReco.C and AnalysisWithExtraBanks.C

1. Adding Hipo files to chain
   1. SKR uses a for loop to add files to the chain with _chain.Add()
   2. AWE just manually adds one hipo file with chain.Add()
2. Both config with config_c12=chain.GetC12Reader()
3. PID selection
   1. SKR takes in final state PIDs from settings::getFinalStatePIDs() ()

How to save vectors to TTree

1. https://root-forum.cern.ch/t/trees-with-vectors/8276/2

Vectors guide

1. https://www.codeguru.com/cplusplus/c-tutorial-a-beginners-guide-to-stdvector-part-1/

#### June 21 2022

Trying to get LundAnalysis.C to save kinematics from hipofile to tree

1. rn, tree is made in root file, but the kinematics are all off
2. Need to make sure that within the while loop we calculate 1 of each kinematic (electron pid, electron momentum, etc) or a particle that contains each value (maybe another branch or like something)

Update 4:53pm

1. The first loop, the event loop, in this macro runs through every event in the hipo file
   1. The second loop, the entry loop, runs through every entry (i think particle) in the event for the MC::Lund bank
2. This means that we want to find which particle each entry is, set particle vectors through if statements
3. Then, in the first loop we can then calculate event variables like Q2, dihadronPt, etc that we can then save to the ttree

Idea

1. Check each particle for pid
   1. if its an up (for example), add its info to a vector through push_back and make a vector for identification using index
   2. then after all particles are run through, search for particles with granddaughters that are pions, and then fill the TTree with this data

Need to do now:

1. Update TTree setting so that we have what we want
   1. current branches are not gonna work - we need P, Ph, Pt, Q2, etc
2. Add calculations (call kinematics) at the end of the first loop to calculate event variables
3. create if statements to figure out which quarks are the ones we care about
   1. use the parent and daughter info from MC92 to figure out what quarks are resulting in what hadrons
      1. add vquarkpid
   2. Use this info to assign parton momentum
      1. vquarkindex should hold the number in place each are so that we can know where the rest of the info is in each vector (place 1, 3, 9, etc)
   3. save this info to the tree
   4. remember to erase the vectors to avoid building billions of entries

For tomorrow - Need to figure out what quarks I want to use for parton momentum

1. it seems that often a quark will fragment into a hadron the decays into a pion, but not fragment directly into the pion - do we want these events? and if so, where do we get parton momentum? from quark that fragments into the non-pion hadron?
2. Finish coding on line 277 - need to write for loops to find quarks, as well as just write kinematic calculations to save to TTree.

#### June 22 2022

1. Need to finish for loop with parton if statements at line 360 (end of file)

#### June 25 2022

1. Add if statement concerning pid=92 where the parent is a quark
   1. Seems like parent of quark being 0 originally may not matter, but will keep for extra security
   2. need to connect the initial quark with pid=92 (they are separated by the second quark listed)
   3. We can probably look at the pid of the parent of pid=92 and check to see if it has parent id of same as quark whose parent is 0
   4. may look like:
      1. if(pid == 92){pid92part = true}
      2. if(pid92part == true)
      3. {
      4. if(MC92parent == quark0parent){fragmentindex = index}
      5. }

## Meeting Notes

#### June 10, 2022

1. Dihadron Cross-Section

## Next Steps

#### Getting a macro ready to process hipo files for pi+pi- endstates

##### Already done:

1. Attempted to make a usable macro by copying "pipi0_process.C" from Greg's macros
   1. See "rowanTestpipluspiminus.C" to see changes made

##### Need to do:

1. Finish altering rowanTestpipluspiminus.C to process pi+pi- endstates

   1. "pipi0_process.C" has a setting "settings.setPostProcessMethod("pipluspi0");"

   2. This seems to be a function that is defined in a different macro, probably in /src

   3. The method "ppluspi0" appears to be defined in /src/PostProcess.C on line 53

   4. I think I will need to define a new method in the function

      ```python
      int PostProcess::Process(TTree *_tree-postprocess)
      ```

   5. Since pipluspi0 is defined similarly, I anticipate needing to include:

      ```
      case PROCESS_ID::pipluspiminus:
      	pipluspiminus(_tree_postprocess);
      	break;
      ```

   6. Furthermore, I need to copy and then alter the PostProcess::pipi0 object to work for pi+pi-

      1. I'm not yet sure what needs to be changed as I do not know what a lot of the variables refer to, so I will probably need to consult Greg

#### Specific changes needed for PostProcess.C

1. Figure out what _piPID refers to on line 411
2. figure out what at(i) does
3. find how to get z of Pi+ and Pi- like on 421
4. Find how to get xF of Pi+ and Pi- like 426
5. Maybe need to initialize an object for mass to replace Mgg on line 338, but it seems like Mdihadron covers this
6. 

## Guesses about Greg's code

#### Kinematics.C

1. Defines a bunch of variables
   1. Q2, x, Px, Py, Pz, Pt, P, E, cth, y, nu, W, th, eta, phi, phi_h, phi_R, com_th
2. Calculates each of these variables by taking in data from hipo files



#### Constants.h

1. Defines a bunch of constants for use in calculations

#### PostProcess.C

1.  

#### makeJobs.sh

1. Takes hipo files from a given hipo directory
2. Takes process name (pi+pi- or pi+pi0)
3. Sends a job to compute cluster
   1. processes a bunch of hipo files according to the process you provide it
   2. Outputs the resulting root files into output directory

Objects within file:

1. workdir - your work directory
2. hipodir - directory containing hipo files you want to pull from
3. outputdir - directory you want to dump root files into
4. rootname - name of file that will be output like "{rootname}{i}.root"
5. processdir - directory containing process macros
6. processcodename - name of process you want to run (analysis type)

## Misc Notes

1. Fall2018 the electron beam was 10.6 GeV
2. Spring2019 electron beam was 10.2 GeV
3. Need to redo "make install" after changing things within CLAS12root analysis
4. Inbending is -
5. Outbending is +
6. June 6 3:50pm
   1. spring2019 had no failures
   2. fall 2018 may have had failures, maybe re-run
7. HEP conference room password
   1. 7 1 3
      1. 7 being 5 and 2 at the same time

8. JLab Collaboration Meeting - on s
   1. In november - make sure I'm registered at least a week before


## Bug fixes

1. Makefile error (error when make install)

   1. Follow gregs help:

   2. **it looks like your [Makefile.am](http://makefile.am/) changed**

      noinst_PROGRAMS = testexternalsBUILT_SOURCES = \
       testexternals.ctestexternals_LDADD = \
       [libclas12ana.la](http://libclas12ana.la/)testexternals.c:
      	echo "//*** this is a generated file. Do not commit, do not edit" > $@
      	echo "int main()" >> $@
      	echo "{" >> $@
      	echo " return 0;" >> $@
      	echo "}" >> $@

      

      **comment all of these**

2. 