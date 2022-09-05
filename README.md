# MCLundAnalysis
GitHub Repository for processing Monte-Carlo generated data from CLAS12, and analyzing that data through the calculation of affinity
## LundAnalysis.C
1. This file is the core of this program. It can read hipo files and output the kinematics needed for calculating affinity
1. How it works:
    1. The clas12reader is employed to parse through the MC::Lund bank of a hipo file, and reads through event by event, and then particle by particle within each event.
    1. For each particle, the pid, parentID and DaughterID are checked in order to identify what the particles are and where they came from.
    1. Cuts are made on the kinematics (1 < Q2 < 100) as well as on the type of event. You can alter the program to select the type of endstate you want, but currently the program is written for pi+pi- endstates.
    1. If there is not a suitable pi+pi- pair in the event, the event is skipped
    1. After passing cuts, certain kinematics are saved to vectors which are used to calculate means for each bin
    1. These means are saved to TTrees depending on which variable is being binned, and then the trees are written to a root file

## AffinityHistos.ipynb
1. This notebook allows for the conversion of TTrees into numpy arrays which can be fed into the affinity tensorflow model to calculate affinity values
