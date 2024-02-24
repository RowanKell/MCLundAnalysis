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
//    Class Stuff
//

class MCParticle
{
    public:
    
    //Lund bank variables
    int pid = 0;
    int id = 0;
    double px = 0;
    double py = 0;
    double pz = 0;
    int daughter = 0;
    int parent = 0;
    double mass = 0;
    double P = 0;
    double E = 0;
    double vz = 0;
    //TLorentzVector
    TLorentzVector lv;
    
    //Calculations
    double Pt = 0;
    TVector2 PtVect;
    
    void inputPxPyPzM(double _px, double _py, double _pz, double _m);
    
    void SetParentDaughter(double _parent,double _daughter);
    
    void fillParticle(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz);
    
    void Calculate();
};

void MCParticle::inputPxPyPzM(double _px, double _py, double _pz, double _m)
{
    px = _px;
    py = _py;
    pz = _pz;
    mass = _m;
    
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
    lv.SetPxPyPzE(px,py,pz,E);
    
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
    
}
void MCParticle::SetParentDaughter(double _parent,double _daughter)
{
    parent = _parent;
    daughter = _daughter;
}
void MCParticle::Calculate()
{
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
    lv.SetPxPyPzE(px,py,pz,E);
    
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
}
void MCParticle::fillParticle(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz)
{
    id = _id;
    pid = _pid;
    px = _px;
    py = _py;
    pz = _pz;
    daughter = _daughter;
    parent = _parent;
    mass = _mass;
    vz = _vz;
    
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
    lv.SetPxPyPzE(px,py,pz,E);
    
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
}

class MultiParticle : public MCParticle
{
    public:
    
    vector<int> v_id;
    vector<int> v_pid;
    vector<double> v_px;
    vector<double> v_py;
    vector<double> v_pz;
    vector<int> v_daughter;
    vector<int> v_parent;
    vector<double> v_mass;
    vector<double> v_P;
    vector<double> v_E;
    vector<double> v_vz;
    vector<TLorentzVector> v_lv;
    vector<double> v_Pt;
    vector<TVector2> v_PtVect;
    
    void update(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz, double _momentum, double _energy, TLorentzVector _lv, double _pt, TVector2 _ptvect)
    {
        v_id.push_back(_id);
        v_pid.push_back(_pid);
        v_px.push_back(_px);
        v_py.push_back(_py);
        v_pz.push_back(_pz);
        v_daughter.push_back(_daughter);
        v_parent.push_back(_parent);
        v_mass.push_back(_mass);
        v_vz.push_back(_vz);
    
        v_P.push_back(_momentum);
        v_E.push_back(_energy);
    
        v_lv.push_back(_lv);
    
        v_Pt.push_back(_pt);
        v_PtVect.push_back(_ptvect);
    }
    
};

class Pion : public MultiParticle
{
    public:
    
    int select_id = -999;
    
};
class Quark : public MultiParticle
{
    public:
    
    int initial_id = -999;
    int final_id = -999;
};

// 
//    Main body of analysis function
//

int ClassTest3()
{
    cout << "Entering LundAnalysis()";
    gROOT->ProcessLine("#include <vector>");
    
    auto hipofile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo";
    auto rootfile = "/../OutputFiles/Lund_8_15/file1.root";
    
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

    //MC::Lund bank entries
    int pid;
    int id;
    double px;
    double py;
    double pz;
    int daughter;
    int parent;
    double mass;
    double vz;

    //Calculated SIDIS kinematics
    double cth;
    double Q2;
    double x;
    double pt;
    double z_h = -999;
    double zpiplus = -999;
    double zpiminus = -999;
    
    //Affinity ratios
    double R0;
    double R1;
    double R2;
    int qparent;
    int diparent;
    double s;
    double y;

    //Dihadron kinematics
    double Mdihadron;
    double dihadronpt;
    double Pdihadron;

    //Intializing MCParticles
    MCParticle electron;
    MCParticle proton;
    MCParticle photon;
    MCParticle Lund;
    
    Pion piplus;
    Pion piminus;
    
    Quark quark;
    
    MultiParticle diquark;
    MultiParticle Hadron;
    
    //MCParticle variables
    double P;
    double E;
    
    //Initializing particle vectors
    TLorentzVector q;
    TLorentzVector init_electron;
    TLorentzVector init_target;
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
//    TLorentzVector diquark;
    TLorentzVector lundstring;
//    TLorentzVector photon;
//    TLorentzVector proton;
    
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
    
    // Checking ki for gauss fit
    TLorentzVector PFki;
    double PFkix = -999;
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
    
    //Checking for low momentum quarks
    double kffrac;
    
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
    auto ivz = config_c12->getBankOrder(idx_MCLund,"vz");
    
    //Creating vectors to fill using push_back in loop
    
    vector<int> vdiquarklist;
    vector<int> vhadronlist;
    vector<int> vquarklist;
    vquarklist = {-8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8};
    vdiquarklist = {1103, 2101, 2103, 2203, 3101, 3103, 3201, 3203, 3303, 4101, 4103, 4201, 4203, 4301, 4303, 4403, 5101, 5103, 5201, 5203, 5301, 5303, 5401, 5403, 5503};
    vhadronlist = {-3122, -211, 111, 211, 1114, 2114, 2212, 2214, 2224, 3112, 3114, 3122, 3214, 3222, 3224, 3312, 3324, -323, -313, -213, 113, 213, 221, 223, 310, 313, 323, 331, 333};
    
    vector<int> extra_pid;
    
//    std::vector<float> venergy;
    //Making new MC tree
    TTree *t = new TTree("tree_MC","Tree with MC data");

    t->Branch("z_h",&z_h);
    t->Branch("x",&x);
    t->Branch("pt",&pt);
    t->Branch("Q2",&Q2);
    t->Branch("Ph",&Pdihadron);
    t->Branch("Mdihadron",&Mdihadron); //dihadron mass
    t->Branch("R0",&R0); //initial parton momentum
    t->Branch("R1",&R1); //final parton momentum
    t->Branch("R2",&R2);
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
    t->Branch("extra",&extra_pid);
    t->Branch("kffrac", &kffrac);
    
    //Tell the user that the loop is starting
    cout << "Start Event Loop" << endl;
    
    int event_count = 0;
    int pi_event_count = 0;
        
    //now get reference to (unique)ptr for accessing data in loop
    //this will point to the correct place when file changes
    //
    //This line comes from AnalysisWithExtraBanks.C
    auto& c12=chain.C12ref();
    
    //Loop over all events in the file
    while(chain.Next()==true){
        event_count += 1;
        if(c12->getDetParticles().empty()) {
            cout << "continuing past event #" << event_count << '\n';
            continue;
    }
        if(pi_event_count > 10) {
            cout << "Breaking at event #" << event_count << '\n';
            break;
        }

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

            //
            //Kinematics
            // 


            //Setting scattered electron
            if(pid==11 && parent==1){
                electron.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
            }
            //pi+
            else if(pid==pipluspid){
                piplus.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                piplus.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, piplus.P, piplus.E, piplus.lv, piplus.Pt, piplus.PtVect);
            }
            //pi-
            else if(pid==piminuspid){
                piminus.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                piminus.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, piplus.P, piplus.E, piplus.lv, piplus.Pt, piplus.PtVect);
            }
            //all quarks
            else if(std::count(vquarklist.begin(), vquarklist.end(), pid)){
                quark.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                quark.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, quark.P, piplus.E, piplus.lv, piplus.Pt, piplus.PtVect);
            }
            //MCParticle
            else if(pid==92 || pid == 91){
                Lund.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
            }
            else if(std::count(vhadronlist.begin(), vhadronlist.end(), pid)) {
                Hadron.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                Hadron.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, piplus.P, piplus.E, piplus.lv, piplus.Pt, piplus.PtVect);
            }
            else if(std::count(vdiquarklist.begin(), vdiquarklist.end(), pid)){
                diquark.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                diquark.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz, piplus.P, piplus.E, piplus.lv, piplus.Pt, piplus.PtVect);
            }
            else if(pid == 22){
                photon.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
            }
            else if(id == 2){
                proton.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
            }
            else{extra_pid.push_back(pid);}
        }
        
        if(quark.v_id.size() > 2) {
            continue;
        }
        
        //Skipping events with multiple quarks as I can't extract momentum from these events yet
        for(int i = 0; i < piplus.v_id.size(); i++) {
            if(piplus.v_parent[i] == Lund.id) {
                piplus.select_id = i;
            }
        }
        for(int i = 0; i < piminus.v_id.size(); i++) {
            if(piminus.v_parent[i] == Lund.id) {
                piminus.select_id = i;
            }
        }
        //Skip non-dipion events
        if((piplus.select_id || piminus.select_id) == -999) {
            continue;
        }
        pi_event_count += 1;
        //Setting inital beam and target particles
        init_electron.SetPxPyPzE(0, 0, sqrt(electron_beam_energy * electron_beam_energy - electronMass * electronMass), electron_beam_energy);
        init_target.SetPxPyPzE(0, 0, 0, proton.E);
        
        dihadron = piplus.lv + piminus.lv;
        Mdihadron = dihadron.M();
        Pdihadron = dihadron.P();
        q = init_electron - electron.lv; //virtual photon
        cth = cthfunc(electron.px,electron.py,electron.pz);
        Q2 = Q2func(electron_beam_energy,electron.E,cth); //Momentum transfer
        zpiplus = (init_target * piplus.lv) / (init_target * q);
        zpiminus = (init_target * piminus.lv) / (init_target * q);
        z_h = zpiplus + zpiminus;
        s = sfunc(protonMass, electronMass, electron_beam_energy);
        y = yfunc(electron_beam_energy,electron.E);
        x = Q2/s/y; // Bjorken x
        pt = Ptfunc(dihadron.Px(), dihadron.Py()); //hadron transverse momentum

        //For loop for finding quarks that fragment from proton and into hadron
        for(int i = 0; i<quark.v_pid.size(); i++) //quark is from proton target 
        {
            if(quark.v_parent[i] == 0)
	      {
                quark.final_id = i;
          }
        }
        quark.fillParticle(quark.v_id[quark.final_id], quark.v_pid[quark.final_id], quark.v_px[quark.final_id], quark.v_py[quark.final_id], 
                           quark.v_pz[quark.final_id], quark.v_daughter[quark.final_id], quark.v_parent[quark.final_id], quark.v_mass[quark.final_id], quark.v_vz[quark.final_id]);
        kf = quark.lv;
        ki = kf - q;
        //Cut Kinematics
        
        dihadronpt = Ptfunc(dihadron.Px(),dihadron.Py());
        
        gN = q;
        gN += init_target;
        gNBoost = gN.BoostVector();
        gNBoostNeg = -gNBoost;
        
        lv_p1_gN = piplus.lv;
        lv_p2_gN = piminus.lv;
        lv_p1_gN.Boost(gNBoostNeg);
        lv_p2_gN.Boost(gNBoostNeg);
        
        lv_q_gN = q;
        lv_q_gN.Boost(gNBoostNeg);
        
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
        
        kffrac = kf.Pz() / q.Pz();
        
        ki2 = abs(ki * ki);
        kf2 = abs(kf * kf);
        k = kf - q;
        
        R0 = R0func(ki, kf, deltak, Q2);
        R1 = R1func(dihadron, ki, kf);
        R2 = R2func(k, Q2);
        
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
        
        cout << "z_h: " << z_h << '\n';
        cout << "PFkix: " << PFkix << '\n';
        cout << "deltakx: " << deltakx << '\n';
        cout << "Pdihadron: " << Pdihadron << '\n';
    }
    f->Write();
    delete f;
    
    return 0;
}
