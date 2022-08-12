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

class MCParticle
{
    public:
    
    //Lund bank variables
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
    
    E = Efunc(mass, P);
    //TLorentzVector
    TLorentzVector lv;
    lv.SetPxPyPzE(px,py,pz,E)
    
    //Calculated kinematics
    double Pt;
    TVector2 PtVect;
    
    //Calculations
    P = Pfunc(px, py, pz);
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
}

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

TVector2 PtVectfunc(TLorentzVector lv)
{
    TVector2 Pt;
    double Px = lv.Px();
    double Py = lv.Py();
    Pt.SetX(Px);
    Pt.SetY(Py);
    return Pt;
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

double R0func(TLorentzVector ki, TLorentzVector kf, TVector2 deltak, double Q2)
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

int LundAnalysis()
{
    
    gROOT->ProcessLine("#include <vector>");
    
    auto hipofile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo";
    auto rootfile = "OutputFiles/Lund_8_6/file3.root";
    
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
    //
    //Initialization
    //

    //Counters
    int qcount;
    int MCindex;
    int MC92count;
    int qindex;
    int vdiquarksize;
    int pioncount;
    int vhadroncount;

    //MC::Lund bank entries
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

    //Calculated SIDIS kinematics
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

    //Dihadron kinematics
    double Mdihadron;
    double dihadronpt;
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
    
    //Cut Kinematics
    double M_x; //Missing Mass
    double xFpiplus;
    double xFpiminus;
    double nu;
    double W;
    double theta;
    
    bool Mx_cut;
    bool xF_cut;
    bool momentum_cut;
    bool W_cut;
    bool y_cut;
    bool vz_cut;
    bool Q2_cut;

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
    TVector2 deltak; // Transverse light cone vector - (V_x,V_y)
    double deltakx;
    double deltaky;
    double deltak2;
    
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
    TVector2 kfBreitTran;
    TLorentzVector dihadronBreit;
    TVector2 dihadronBreitTran;
    
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
    TVector2 q_T;
    
    
    float pxdiff;
    float pydiff;
    float pzdiff;
    float ediff;
    
    int piplusparent;
    int piminusparent;
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
    
    double ki2;
    double kf2;
    
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
    
//    std::vector<float> venergy;
    //Making new MC tree
    TTree *t = new TTree("tree_MC","Tree with MC data");

    t->Branch("z_h",&z_h);
    t->Branch("x",&x);
    t->Branch("pt",&pt);
    t->Branch("Q2",&Q2);
    t->Branch("Ph",&Pdihadron);
    t->Branch("P",&Ptarget); //target momentum
    t->Branch("Mdihadron",&Mdihadron); //dihadron mass
    t->Branch("MC92index",&MC92index); //index for Lund string
    t->Branch("R0",&R0); //initial parton momentum
    t->Branch("R1",&R1); //final parton momentum
    t->Branch("R2",&R2);
    t->Branch("pxdiff",&pxdiff);
    t->Branch("pydiff",&pydiff);
    t->Branch("pzdiff",&pzdiff);
    t->Branch("ediff",&ediff);
    t->Branch("PFkix",&PFkix); //Photon frame partonic momentum for checking ki values
    t->Branch("PFkiy",&PFkiy);
    t->Branch("PFkiz",&PFkiz);
    t->Branch("PFkit",&PFkit);
    t->Branch("PFkfx",&PFkfx);
    t->Branch("PFkfy",&PFkfy);
    t->Branch("PFkfz",&PFkfz);
    t->Branch("PFkft",&PFkft);
    t->Branch("deltakx",&deltakx);
    t->Branch("deltaky",&deltaky);
    t->Branch("deltak2",&deltak2);
    t->Branch("ki2",&ki2);
    t->Branch("kf2",&kf2);
    
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
        
        W_cut = false;
        xF_cut = false;
        Mx_cut = false;
        momentum_cut = false;
        vz_cut = false;
        
        qcount = 0;
        qdaughtercount = 0;
        MC92count = 0;
        pioncount = 0;
        vhadroncount = 0;
        event_count += 1;
        
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
            vz = mcparticles->getVz(imc);

            //Kinematics
            // 


            //Setting scattered electron
            if(pid==11 && parent==1){
                electron.SetPxPyPzE(px,py,pz,E);
                vz_electron = vz;
            }
            //pi+
            else if(pid==pipluspid){
                piplus.SetPxPyPzE(px,py,pz,E);
                piplusparent = parent;
                pioncount += 1;
                vz_piplus = vz;
            }
            //pi-
            else if(pid==piminuspid){
                piminus.SetPxPyPzE(px,py,pz,E);
                piminusparent = parent;
                pioncount += 1;
                vz_piminus = vz
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
            else if(pid==92){
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
        
        //Skipping events with multiple quarks as I can't extract momentum from these events yet
        if((piplusparent != MC92index) || (piminusparent != MC92index)) {
            piparent = false;
        }
        else {piparent = true;}
        for(int i = 0; i < vhadroncount; i++) {
            if(vhadronparent[i] == 2) {
                continue;
            }
        }
        if(
            (qcount != 2) || 
            (pioncount > 2) || 
            (piparent == false) || 
            (MC92count != 1)
           ){
            continue;
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
        pt = Ptfunc(dihadron.Px(), dihadron.Py()); //hadron transverse momentum

        //For loop for finding quarks that fragment from proton and into hadron
        for(int i = 0; i<qcount; i++) //quark is from proton target 
        {
            if(vquarkparent[i] == 0)
	      {
                kf.SetPxPyPzE(vquarkPx[i],vquarkPy[i],vquarkPz[i],vquarkenergy[i]);
          }
        }
        
        
        //Cut Kinematics
        M_x = Mxfunc(q,init_target,piplus,piminus);
        
        nu = nufunc(electron_beam_energy,electron.E());
        W = Wfunc(Q2,protonMass,nu);
        
        if((W > 2) && (W < 100)) {
            W_cut = true;
        }
        else {W_cut = false;}
        
        if((Q2 > 1) && (Q2 < 100)) {
            Q2_cut = true;
        }
        else {Q2_cut = false;}
        
        if(abs(vz_electron-vz_piplus)<20 && abs(vz_electron-vz_piminus)<20) {
            vz_cut = true;
        }
        else {vz_cut = false}
        
        if((y > 0) && (y < 0.8)) {
            y_cut = true;
        }
        else {y_cut = false;}
        dihadronpt = Ptfunc(dihadron.Px(),dihadron.Py());
        theta = thetafunc(dihadronpt,dihadron.Pz());
        
        if((piplus.P() > 1.25) && (piplus.P() > 1.25)) {
            momentum_cut = true;
        }
        else {momentum_cut = false;}
        
        gN = q;
        gN += init_target;
        gNBoost = gN.BoostVector();
        gNBoostNeg = -gNBoost;
        
        lv_p1_gN = piplus;
        lv_p2_gN = piminus;
        lv_p1_gN.Boost(gNBoostNeg);
        lv_p2_gN.Boost(gNBoostNeg);
        
        lv_q_gN = q;
        lv_q_gN.Boost(gNBoostNeg);
        
        xFpiplus = xFfunc(lv_p1_gN,lv_q_gN,W);
        xFpiminus = xFfunc(lv_p2_gN,lv_q_gN,W);
        if((xFpiplus > 0) && (xFpiminus > 0)) {
            xF_cut = true;
        }
        
        pxdiff = proton.Px() + photon.Px() - diquark.Px() - kf.Px();
        pydiff = proton.Py() + photon.Py() - diquark.Py() - kf.Py();
        pzdiff = proton.Pz() + photon.Pz() - diquark.Pz() - kf.Pz();
        ediff = proton.E() + photon.E() - diquark.E() - kf.E();
        
        // Breit Frame Kinematics for delta k
        Breit = q;
        Breit_target.SetPxPyPzE(0,0,0,2 * x * protonMass); // E^2 = M^2 + P^2 --> P = 0 so E = M = 2 * x * protonmass
        Breit += Breit_target;
        BreitBoost = Breit.BoostVector();
        BreitBoost = -1 * BreitBoost;
        kfBreit = kf;
        kfBreit.Boost(BreitBoost);
        kfBreitTran = PtVectfunc(kfBreit); //kfbT in delta k calculation - needs to be a transverse light cone vector of form (V_x, V_y)
        
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
        dihadronPFMinus = LightConeMinus(dihadronPF);
        //Virtual Photon
        qPF.Rotate(PFAngle,PFAxis);
        qPFMinus = LightConeMinus(qPF);
        //z_N and q_T
        z_N = dihadronPFMinus / qPFMinus;
        q_T = -1 * dihadronBreitTran / z_N;
        //ki, k, and delta k
        deltak = kfBreitTran - (-1 * z_N * q_T); 
        deltakx = deltak.Px();
        deltaky = deltak.Py();
        deltak2 = deltak * deltak;
        ki = kf - q;
        ki2 = abs(ki * ki);
        kf2 = abs(kf * kf);
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
        PFkit = Ptfunc(PFkix,PFkiy);
        
        PFkfx = PFkf.Px();
        PFkfy = PFkf.Py();
        PFkfz = PFkf.Pz();
        PFkft = Ptfunc(PFkfx,PFkfy);
        if((Mx_cut == true) && (xF_cut == true) && (momentum_cut == true) && (W_cut == true) && (y_cut == true) && (vz_cut == true) && (Q2_cut == true) {
            t->Fill();
        }
        else{continue;}
    }
    f->Write();
    delete f;
    
    return 0;
}
