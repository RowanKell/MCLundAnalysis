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

//
//     Functions for calculating kinematics and ratios
//

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

double PtVectfunc(TLorentzVector lv)
{
    double Px = lv.Px();
    double Py = lv.Py();
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

double R0func(TLorentzVector ki, TLorentzVector kf, TLorentzVector deltak, double Q2)
{
    double init;
    double final;
    double delta;
    
    double mag;
    
    init = abs((ki * ki) / (Q2));
    final = abs((kf * kf) / (Q2));
    delta = abs((deltak * deltak) / (Q2));
    
    //Selecting the max of the ratios
    if(init > final && init > delta){mag = init;}
    else if(final > init && final > delta){mag = final;}
    else if(delta > init && delta > final){mag = delta;}
    return mag;
}

double R1func(TLorentzVector Ph, TLorentzVector ki, TLorentzVector kf)
{
    return (Ph * ki) / (Ph * kf);
}

double R2func(TLorentzVector k, double Q2)
{
    return abs(k * k) / Q2;
}

double Mxfunc(TLorentzVector q, TLorentzVector target, TLorentzVector hadron1, TLorentzVector hadron2)
{
    TLorentzVector lv_Mx;
    lv_Mx = q + target - hadron1 - hadron2;
    return lv_Mx.M();
}
double xFfunc(TLorentzVector p, TLorentzVector q, double W)
{
    return 2 * (p.Vect().Dot(q.Vect())) / (q.Vect().Mag() * W);
}

double nufunc(double E1, double E2)
{
  return E1-E2;
}

double Wfunc(double Q2, double mT, double nu)
{
    return sqrt(-Q2 + pow(mT,2) + 2 * mT * nu);
}

double thetafunc(double pt, double pz)
{
    return abs(atan(pt/pz));
}

double LightConeMinus(TLorentzVector lv)
{
    return (lv.E() - lv.Pz()) / (sqrt(2));
}

double LightConePlus(TLorentzVector lv)
{
    return (lv.E() + lv.Pz()) / (sqrt(2));
}
// 
//    Main body of analysis function
//

int PartonicMomentumTest()
{
    
    gROOT->ProcessLine("#include <vector>");
    
    auto hipofile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo";
    auto rootfile = "../OutputFiles/Lund_8_3/Aug_3_file3.root";
    
    TFile *f = TFile::Open(rootfile,"RECREATE");
    
    HipoChain chain;
    
    //Add file to HipoChain
    chain.Add(hipofile);
    auto config_c12 = chain.GetC12Reader();
    
    //Set PID cuts
    config_c12->addExactPid(11,1);    //exactly 1 electron
//    config_c12->addExactPid(211,1);    //exactly 1 pi+
//    config_c12->addExactPid(-211,1);    //exactly 1 pi-
    config_c12->addExactPid(2212,1);    //exactly 1 proton

    //Constants 
    double electron_beam_energy = 10.6; //(fall2018)
    double electronMass = 0.000511;
    double protonMass = 0.938272;
    double pipluspid = 211;
    double piminuspid =- 211;
    double uppid = 2;
    double downpid = 1;
    double antiuppid = -2;
    double antidownpid = -1;
    
    //Initialization
    int qcount;
    int MCindex;
    int MC92count;
    int qindex;
    int vdiquarksize;
    int pioncount;
    int vhadroncount;
    
    int pid;
    double id;
    double px;
    double py;
    double pz;
    double dihadronpt;
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
    double s;
    double y;
    double init_qE;
    double final_qE;
    double Mdihadron;
    double MC92index;
    double protonE;
    double MC92px;
    double MC92py;
    double MC92pz;
    double MC92mass;
    double MC92parent;
    double MC92daughter;
    double Ptarget;
    double Pdihadron;
    double MC92P;
    double MC92energy;
    double quarkinitP;
    double quarkinitE;  
    double energy;
    
    double kim;
    double kie;
    double kip;
    double kfm;
    double kfe;
    double kfp;
    double qdaughtercount;
    
    double qdiff;
    
    double M_x; //Missing Mass
    double xFpiplus;
    double xFpiminus;
    double nu;
    double W;
    double theta;
    
    bool Mxcut;
    bool xFcut;

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
    
    TLorentzVector lv_p1_gN;
    TLorentzVector lv_p2_gN;
    TLorentzVector lv_q_gN;
    TLorentzVector gN;
    TVector3 gNBoost;
    TVector3 gNBoostNeg;
    //For checking momentum conservation / if the lundstring contains momentum for all hadrons
    TLorentzVector diquark;
    TLorentzVector lundstring;
    TLorentzVector photon;
    TLorentzVector proton;
    
    //Breit frame variables
    TLorentzVector Breit;
    TLorentzVector Breit_target;
    TVector3 BreitBoost;
    TLorentzVector kfBreit;
    double kfBreitTran;
    TLorentzVector dihadronBreit;
    double dihadronBreitTran;
    
    //Photon Frame variables
    TLorentzVector PFFrame;
    TVector3 PFBoost;
    TLorentzVector qPF;
    TVector3 qPFVect;
    TVector3 qPFVectUnit;
    TVector3 zAxis(0,0,1);
    double PFAngle;
    TVector3 PFAxis;
    
    double qPFMinus;
    TLorentzVector dihadronPF;
    double dihadronPFMinus;
    
    double z_N;
    TLorentzVector q_T;
    
    
    float pxdiff;
    float pydiff;
    float pzdiff;
    float ediff;
    
    int piplusparent;
    int piminusparent;
    int pifirstparent;
    int pisecondparent;
    
    bool piparent;
    int event_count;
    
    // Checking ki for gauss fit
    TLorentzVector PFki;
    double PFkix;
    double PFkiy;
    double PFkiz;
    double PFkit;
    TLorentzVector PFkf;
    double PFkfx;
    double PFkfy;
    double PFkfz;
    double PFkft;
    
    double qPFx;
    double qPFy;
    double qPFz;
    double qPFpt;
    
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
//    auto iE=config_c12->getBankOrder(idx_MCLund,"energy");
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
    std::vector<float> vquarkenergy;
    //MC92 particle (pid = 92) vectors
//    std::vector<float> vMC92pid;
//    std::vector<float> vMC92parent;
//    std::vector<float> vMC92daughter;
    std::vector<bool> initparent;
    std::vector<float> vinitquarkindex;
//    std::vector<float> vMC92index; //keeps track of position in vector
    std::vector<int> vdiquarklist;
    std::vector<int> vhadronlist;
    
    //Vectors for comparing momentum
    std::vector<float> vpxcompare;
    std::vector<float> vpycompare;
    std::vector<float> vpzcompare;
    std::vector<float> vecompare;
    std::vector<float> vpidcompare;
    std::vector<float> vdaughtercompare;
    std::vector<float> vparentcompare;
    std::vector<int> vhadronparent;
    
    int qparent0;
    
//    std::vector<float> venergy;
    //Making new MC tree
    TTree *t = new TTree("tree_MC","Tree with MC data");


    t->Branch("PFkix",&PFkix); //Photon frame partonic momentum for checking ki values
    t->Branch("PFkiy",&PFkiy);
    t->Branch("PFkiz",&PFkiz);
    t->Branch("PFkit",&PFkit);
    t->Branch("PFkfx",&PFkfx);
    t->Branch("PFkfy",&PFkfy);
    t->Branch("PFkfz",&PFkfz);
    t->Branch("PFkft",&PFkft);
    t->Branch("qparent0",&qparent0);
    t->Branch("qPFx",&qPFx);
    t->Branch("qPFy",&qPFy);
    t->Branch("qPFz",&qPFz);
    t->Branch("qPFpt",&qPFpt);
    
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
        qdaughtercount = 0;
        MC92count = 0;
        pioncount = 0;
        vhadroncount = 0;
        event_count += 1;
        qparent0 = 0;
        
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
        vquarkenergy.clear();
        
        vhadronparent.clear();
        
        vdiquarklist = {1103, 2101, 2103, 2203, 3101, 3103, 3201, 3203, 3303, 4101, 4103, 4201, 4203, 4301, 4303, 4403, 5101, 5103, 5201, 5203, 5301, 5303, 5401, 5403, 5503};
        vdiquarksize = vdiquarklist.size();
        
        vhadronlist = {-3122, -211, 111, 211, 1114, 2114, 2212, 2214, 2224, 3112, 3114, 3122, 3214, 3222, 3224, 3312, 3324, -323, -313, -213, 113, 213, 221, 223, 310, 313, 323, 331, 333};
        
        //Loop over MC::Lund entries in this event using its ID = idx_MCLund
        //Get PID from its id = iPid
        for(auto imc=0;imc<c12->getBank(idx_MCLund)->getRows();imc++){
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

            //Kinematics
            // 

            if(pioncount == 0){
                //pi+
                if(pid==pipluspid){
                    piplus.SetPxPyPzE(px,py,pz,E);
                    pifirstparent = parent;
                    pioncount += 1;
                }
                //pi-
                else if(pid==piminuspid){
                    piminus.SetPxPyPzE(px,py,pz,E);
                    pifirstparent = parent;
                    pioncount += 1;
                }
            }
            if(pioncount == 1){
                //pi+
                if(pid==pipluspid && parent == pifirstparent){
                    piplus.SetPxPyPzE(px,py,pz,E);
                    pisecondparent = parent;
                    pioncount += 1;
                }
                //pi-
                else if(pid==piminuspid && parent == pifirstparent){
                    piminus.SetPxPyPzE(px,py,pz,E);
                    pisecondparent = parent;
                    pioncount += 1;
                }
            }
            //Setting scattered electron
            if(pid==11 && parent==1){
                electron.SetPxPyPzE(px,py,pz,E);
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
                vquarkenergy.push_back(E);
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
                vquarkenergy.push_back(E);
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
                vquarkenergy.push_back(E);
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
                vquarkenergy.push_back(E);
            }
            //MCParticle
            else if(pid==92 or pid == 91){
                MC92count += 1;
                MC92parent = parent;
                MC92daughter = daughter;
                MC92px = px;
                MC92py = py;
                MC92pz = pz;
                MC92mass = mass;
                MC92index = id;
                MC92energy = E;
                lundstring.SetPxPyPzE(px,py,pz,E);
//                vMC92index.push_back(index);
            }
            else if(std::count(vhadronlist.begin(), vhadronlist.end(), pid) && parent == 2) {
                vhadronparent.push_back(parent);
                vhadroncount += 1;
            }
            else if(std::count(vdiquarklist.begin(), vdiquarklist.end(), pid) && parent == 2){
                diquark.SetPxPyPzE(px,py,pz,E);
            }
            else if(pid == 22){
                photon.SetPxPyPzE(px,py,pz,E);
            }
            else if(id == 2){
                proton.SetPxPyPzE(px,py,pz,E);
            }
        }
        
        
        //Setting inital beam and target particles
        init_electron.SetPxPyPzE(0, 0, sqrt(electron_beam_energy * electron_beam_energy - electronMass * electronMass), electron_beam_energy);
        protonE = Efunc(0,protonMass);
        init_target.SetPxPyPzE(0, 0, 0, protonE);
        Ptarget = init_target.P();
        
        dihadron = piplus + piminus;
        Mdihadron = dihadron.M();
        Pdihadron = dihadron.P();
        q = init_electron - electron; //virtual photon
        qdiff = (q - photon).P();
        cth = cthfunc(electron.Px(),electron.Py(),electron.Pz());
        Q2 = Q2func(electron_beam_energy,electron.E(),cth); //Momentum transfer
        zpiplus = (init_target * piplus) / (init_target * q);
        zpiminus = (init_target * piminus) / (init_target * q);
        z_h = zpiplus + zpiminus;
        s = sfunc(protonMass, electronMass, electron_beam_energy);
        y = yfunc(electron_beam_energy,E);
        x = Q2/s/y; // Bjorken x
        pt = PtVectfunc(dihadron); //hadron transverse momentum

        //For loop for finding quarks that fragment from proton and into hadron
        for(int i = 0; i<qcount; i++) //quark is from proton target 
        {
            if(vquarkparent[i] == 0)
	      {
                kf.SetPxPyPzE(vquarkPx[i],vquarkPy[i],vquarkPz[i],vquarkenergy[i]);
                qparent0 += 1;
          }
        }
        
        

        
        // Breit Frame Kinematics for delta k
        Breit = q;
        Breit_target.SetPxPyPzE(0,0,0,2 * x * protonMass); // E^2 = M^2 + P^2 --> P = 0 so E = M = 2 * x * protonmass
        Breit += Breit_target;
        BreitBoost = Breit.BoostVector();
        BreitBoost = -1 * BreitBoost;
        kfBreit = kf;
        kfBreit.Boost(BreitBoost);
//        kfBreitTran = Ptfunc(kfBreit.Px(),kfBreit.Py()); //kfbT in delta k calculation
        
        dihadronBreit = dihadron;
        dihadronBreit.Boost(BreitBoost);
        dihadronBreitTran = PtVectfunc(dihadronBreit); //PBbT in qT part of delta k calculation
        
        PFFrame = q + init_target;
        PFBoost = PFFrame.BoostVector();
        PFBoost = -1 * PFBoost;
        qPF = q;
        qPF.Boost(PFBoost);
        qPFVect = qPF.Vect();
        qPFVectUnit = qPFVect.Unit();
        PFAngle = qPFVectUnit.Angle(zAxis);
        PFAxis = qPFVectUnit.Cross(zAxis);
        //To rotate -> vector.Rotate(PFAngle,PFAxis);
        
        //
        //Photon Frame
        //
        
        //Dihadron
        dihadronPF = dihadron;
        dihadronPF.Boost(PFBoost);
        dihadronPF.Rotate(PFAngle,PFAxis);
//        dihadronPFMinus = LightConeMinus(dihadronPF);
        //Virtual Photon
        qPF.Rotate(PFAngle,PFAxis);
        qPFpt = PtVectfunc(qPF);
        qPFx = qPF.Px();
        qPFy = qPF.Py();
        qPFz = qPF.Pz();
//        qPFMinus = LightConeMinus(qPFMinus);
        //z_N and q_T
        z_N = dihadronPFMinus / qPFMinus;
//        q_T = -1 * dihadronBreitTran / z_N;
        //ki, k, and delta k
//        deltak = kfBreitTran - (-1 * z_N * q_T); 
        ki = kf - q;
        k = kf - q;
        
        R0 = R0func(ki, kf, deltak, Q2);
        R1 = R1func(dihadron, ki, kf);
        R2 = R2func(k, Q2);
            
        kim = ki.M();
        kie = ki.E();
        kip = ki.P();
        
        kfm = kf.M();
        kfe = kf.E();
        kfp = kf.P();
        
        PFki = ki;
        PFki.Rotate(PFAngle,PFAxis);
        PFkf = kf;
        PFkf.Rotate(PFAngle,PFAxis);
        PFkix = PFki.Px();
        PFkiy = PFki.Py();
        PFkiz = PFki.Pz();
        PFkit = PtVectfunc(PFki);
        
        PFkfx = PFkf.Px();
        PFkfy = PFkf.Py();
        PFkfz = PFkf.Pz();
        PFkft = PtVectfunc(PFkf);
        if(pifirstparent != pisecondparent || pioncount != 2){
            continue;
        }
        t->Fill();
    }
    f->Write();
    delete f;
    
    return 0;
}
