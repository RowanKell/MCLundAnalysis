#include <iostream>
#include <vector>
#include <cmath>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

/*NEED TO CALCULATE AND SAVE TO TREE:
1. P (Momentum of initial Hadron - Proton target)
2. q (momentum transfer of incident lepton - momentum of virtual photon)
3. Ph (momentum of final hadron - fragment hadron (dihadron))
4. k (kf - q)
5. kf (hadronizing parton(s) momentum) (really just Pquarkf)
6. ki (initial partons momentum) (really just Pquarki)
7. delta k^2
8. zh 
9. Bjorken x
10. Q^2
11. pT (transverse hadron momentum)
*/
double Pfunc(double Px, double Py, double Pz)
{
    return sqrt(Px*Px + Py*Py + Pz*Pz);
}
double Efunc(double M, double P)
{
    return sqrt(M * M + P * P);
}
double Ptfunc(double Px, double Py)
{
    return sqrt(Px*Px + Py*Py);
}

double cthfunc(double Px, double Py, double Pz)
{
    double Pt = sqrt(Px*Px + Py*Py);
    return Pz / sqrt(Pz * Pz + Pt * Pt);
}

double Q2func(double E1, double E2, double cth)
{
    return 2 * E1 * E2 * (1.0 - cth);
}

double yfunc(double E1,double E2)
{
    return (E1 - E2) / E1;
}
double sfunc(double M1,double M2,double E) //proton mass is M1, electron mass is M2, E is electron beam energy
{
    return M1 * M1 + M2 * M2 + 2 * M1 * E;
}
int LundAnalysis()
{
    
    gROOT->ProcessLine("#include <vector>");
    
    auto hipofile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo";
    auto rootfile = "testfile.root";
    
    TFile *f = TFile::Open(rootfile,"RECREATE");
    
    HipoChain chain;
    
    //Add file to HipoChain
    chain.Add(hipofile);
    auto config_c12 = chain.GetC12Reader();
    
    //Set PID cuts
    config_c12->addExactPid(11,1);    //exactly 1 electron
    config_c12->addExactPid(211,1);    //exactly 1 pi+
    config_c12->addExactPid(-211,1);    //exactly 1 pi-
    config_c12->addExactPid(2212,1);    //exactly 1 proton

    //Constants 
    double electron_beam_energy = 10.6; //(fall2018)
    double electronMass = 0.000511;
    double protonMass = 0.938272;
    double pipluspid = 211;
    double piminuspid =-211;
    double uppid = 2;
    double downpid =1;
    double antiuppid = -2;
    double antidownpid = -1;
    
    //Initialization
    int qcount;
    int MCindex;
    
    int pid;
    double id;
    double px;
    double py;
    double pz;
    double daughter;
    double parent;
    double mass;
    double P;
    double E;
    double cth;
    double Q2;
    double x;
    double pt;
    double z_h;
    double zpiplus;
    double zpiminus;
    double R0;
    double R1;
    double R2;
    double qparent;
    double diparent;
    double MC92parent;
    double s;
    double y;
    double init_qE;
    double final_qE;
    double Mdihadron;
    double MC92index;
    double protonE;
    
//    double energy;

    //Initializing particle vectors
    TLorentzVector q;
    TLorentzVector init_electron;
    TLorentzVector init_target;
    TLorentzVector electron;
    TLorentzVector piplus;
    TLorentzVector piminus;
    TLorentzVector dihadron;
    
    TLorentzVector k;
    TLorentzVector kf;
    TLorentzVector ki;
    TLorentzVector deltak;
    
    //Add MC::Lund bank for taking Lund data
    auto idx_MCLund= config_c12->addBank("MC::Lund");
    //Add a few items
    auto iPid=config_c12->getBankOrder(idx_MCLund,"pid");
    auto ipx=config_c12->getBankOrder(idx_MCLund,"px"); 
    auto ipy=config_c12->getBankOrder(idx_MCLund,"py");
    auto ipz=config_c12->getBankOrder(idx_MCLund,"pz");
    auto idaughter=config_c12->getBankOrder(idx_MCLund,"daughter");
    auto iparent=config_c12->getBankOrder(idx_MCLund,"parent");
    auto imass=config_c12->getBankOrder(idx_MCLund,"mass");
//    auto ivx=config_c12->getBankOrder(idx_MCLund,"vx");
//    auto ivy=config_c12->getBankOrder(idx_MCLund,"vy");
//    auto ivz=config_c12->getBankOrder(idx_MCLund,"vz");
    auto iE=config_c12->getBankOrder(idx_MCLund,"energy");
//    auto ipy=config_c12->getBankOrder(idx_MCLund,"type");
    
    //Creating vectors to fill using push_back in loop
    
    //hadron vectors
    std::vector<float> vpid;
    std::vector<float> vpx;
    std::vector<float> vpy;
    std::vector<float> vpz;
    std::vector<float> vdaughter;
    std::vector<float> vparent;
    std::vector<float> vmass;
    //quark vectors
    std::vector<float> vquarkindex;
    std::vector<float> vquarkPx;
    std::vector<float> vquarkPy;
    std::vector<float> vquarkPz;
    std::vector<float> vquarkparent;
    std::vector<float> vquarkdaughter;
    std::vector<float> vquarkmass;
    //MC92 particle (pid = 92) vectors
    std::vector<float> vMC92pid;
    std::vector<float> vMC92parent;
    std::vector<float> vMC92daughter;
//    std::vector<float> vMC92index; //keeps track of position in vector
    
//    std::vector<float> venergy;
    //Making new MC tree
    TTree *t = new TTree("tree_MC","Tree with MC data");
    
    //Assigning variables to branches
/*    t->Branch("pid",&pid);
    t->Branch("px",&px);
    t->Branch("py",&py);
    t->Branch("pz",&pz);
    t->Branch("daughter",&daughter);
    t->Branch("parent",&parent);
    t->Branch("mass",&mass);
*/
    t->Branch("z_h",&z_h);
//    t->Branch("R0",&R0);
//    t->Branch("R1",&R1);
//    t->Branch("R2",&R2);
    t->Branch("x",&x);
    t->Branch("pt",&pt);
    t->Branch("Q2",&Q2);
    t->Branch("Mdihadron",&Mdihadron); //dihadron mass
    t->Branch("qindex",&qindex);
    t->Branch("MC92index",&MC92index);
    
//    t->Branch("qparent",&qparent);
//    t->Branch("diparent",&diparent);
//    t->Branch("MC92parent",&MC92parent);
//    t->Branch("",&);
//    t->Branch("",&);
//    t->Branch("energy",&energy);
    
    //Tell the user that the loop is starting
    cout << "Start Event Loop" << endl;
        
    //now get reference to (unique)ptr for accessing data in loop
    //this will point to the correct place when file changes
    //
    //This line comes from AnalysisWithExtraBanks.C
    auto& c12=chain.C12ref();
    
    //Loop over all events in the file
    while(chain.Next()==true){
        if(c12->getDetParticles().empty())
            continue;
        
        qcount = 0;
        MC92index = 0;
        
        vpid.clear();
        vpx.clear();
        vpy.clear();
        vpz.clear();
        vdaughter.clear();
        vparent.clear();
        vmass.clear();
        
        vquarkPx.clear();
        vquarkPy.clear();
        vquarkPz.clear();
        vquarkparent.clear();
        vquarkdaughter.clear();
        vquarkmass.clear();
	vquarkindex.clear();
        
        vMC92pid.clear();
        vMC92parent.clear();
        vMC92daughter.clear();
//        vMC92index.clear();
        
//      venergy.clear();
        //Loop over MC::Lund entries in this event using its ID = idx_MCLund
        //Get PID from its id = iPid
        for(auto imc=0;imc<c12->getBank(idx_MCLund)->getRows();imc++){
/*            pid = c12->getBank(idx_MCLund)->getInt(iPid,imc);
            px = c12->getBank(idx_MCLund)->getInt(iPid,imc);
            py = c12->getBank(idx_MCLund)->getInt(iPid,imc);
            pz = c12->getBank(idx_MCLund)->getInt(iPid,imc);
            daughter = c12->getBank(idx_MCLund)->getInt(iPid,imc);
            parent = c12->getBank(idx_MCLund)->getInt(iPid,imc);
            mass = c12->getBank(idx_MCLund)->getInt(iPid,imc);
            energy = c12->getBank(idx_MCLund)->getInt(iPid,imc);
*/
            auto mcparticles = c12->mcparts();
            
            id = mcparticles->getIndex(imc);
            pid = mcparticles->getPid(imc);
            px = mcparticles->getPx(imc);
            py = mcparticles->getPy(imc);
            pz = mcparticles->getPz(imc);
            daughter = mcparticles->getDaughter(imc);
            parent = mcparticles->getParent(imc);
            mass = mcparticles->getMass(imc);
            P = Pfunc(px,py,pz);
            E = Efunc(mass,P);
//            energy = mcparticles->getEnergy(imc);

            //Kinematics
            // 


            //Setting scattered electron
            if(pid==11 && parent==1){
                electron.SetPxPyPzE(px,py,pz,E);
            }
            //pi+
            else if(pid==pipluspid){
                piplus.SetPxPyPzE(px,py,pz,E);
            }
            //pi-
            else if(pid==piminuspid){
                piminus.SetPxPyPzE(px,py,pz,E);
            }
            //inital up
            else if(pid==uppid){
	        qcount += 1;
                vquarkindex.push_back(id);
                vquarkPx.push_back(px);
                vquarkPy.push_back(py);
                vquarkPz.push_back(pz);
                vquarkparent.push_back(parent);
                vquarkdaughter.push_back(daughter);
                vquarkmass.push_back(mass);
            }
            //down
            else if(pid==downpid){
	        qcount += 1;
                vquarkindex.push_back(id);
                vquarkPx.push_back(px);
                vquarkPy.push_back(py);
                vquarkPz.push_back(pz);
                vquarkparent.push_back(parent);
                vquarkdaughter.push_back(daughter);
                vquarkmass.push_back(mass);
            }
            //anti-up
            else if(pid==antiuppid){
	        qcount += 1;
                vquarkindex.push_back(id);
                vquarkPx.push_back(px);
                vquarkPy.push_back(py);
                vquarkPz.push_back(pz);
                vquarkparent.push_back(parent);
                vquarkdaughter.push_back(daughter);
                vquarkmass.push_back(mass);
            }    
            //anti-down
            else if(pid==antidownpid){
	        qcount += 1;
                vquarkindex.push_back(id);
                vquarkPx.push_back(px);
                vquarkPy.push_back(py);
                vquarkPz.push_back(pz);
                vquarkparent.push_back(parent);
                vquarkdaughter.push_back(daughter);
                vquarkmass.push_back(mass);
            }
            //MCParticle
            else if(pid==92){
                MC92index += 1;
                vMC92parent.push_back(parent);
                vMC92daughter.push_back(daughter);
//                vMC92index.push_back(index);
            }
        }
        //Kinematics

        //Setting inital beam and target particles
        init_electron.SetPxPyPzE(0, 0, sqrt(electron_beam_energy * electron_beam_energy - electronMass * electronMass), electron_beam_energy);
        protonE = Efunc(0,protonMass);
        init_target.SetPxPyPzE(0, 0, 0, protonE);
        
        dihadron = piplus + piminus;
        Mdihadron = dihadron.M();
        q = init_electron - electron; //virtual photon
        cth = cthfunc(electron.Px(),electron.Py(),electron.Pz());
        Q2 = Q2func(electron_beam_energy,electron.E(),cth); //Momentum transfer
        zpiplus = (init_target * piplus) / (init_target * q);
        zpiminus = (init_target * piminus) / (init_target * q);
        z_h = zpiplus + zpiminus;
        s = sfunc(protonMass, electronMass, electron_beam_energy);
        y = yfunc(electron_beam_energy,E);
        x = Q2/s/y; // Bjorken x
        pt = Ptfunc(dihadron.Px(), dihadron.Py()); //hadron transverse momentum
//        init_qE = E();
//        final_qE = E();
//        ki.SetPxPyPzE(); //parton initial
//        kf.SetPxPyPzE(); //parton final
//        deltak = ;
//        k = kf - q; //parton

        //For loop for finding quarks that fragment from proton and into hadron
        for(int i = 0; i<qindex; i++)
        {
            if(vquarkparent[i] == 0)
	      {
		initparent[i] = true;
              } else
	      {
		initparent[i] = false;
	      }
	    if(
        }
        
        t->Fill();
    }
    f->Write();
    delete f;
    
    return 0;
}
