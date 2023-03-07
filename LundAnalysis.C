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
    double M2 = 0;
    if(M < 0) {M2 = M * (-M);}
    else {M2 = M * M;}
    return sqrt(M2 + P * P);
}

double EVirtualfun(double M, double P)
{
    return sqrt(M * (-M) + P * P);
}
double Ptfunc(double Px, double Py)
{
    return sqrt(Px*Px + Py*Py);
}
double Ptfunc(TLorentzVector lv) {
    double x = lv.Px();
    double y = lv.Py();
    return sqrt(x * x + y * y);
}
double Ptfunc(TVector2 v) {
    double x = v.X();
    double y = v.Y();
    return sqrt(x * x + y * y);
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
    
    double mag = -999;
    
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

double meanfunc(vector<double> v)
{
    double sum = 0;
    for(int i = 0; i < v.size(); i++) {
        sum += v[i];
    }
    return sum / v.size();
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
    
    void setVectors();
};

void MCParticle::inputPxPyPzM(double _px, double _py, double _pz, double _m)
{
    px = _px;
    py = _py;
    pz = _pz;
    mass = _m;
    
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
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
void MCParticle::setVectors()
{
    lv.SetPxPyPzE(px,py,pz,E);
    
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
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
    vector<double> v_vz;
    
    
    void update(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz)
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
    }
    
};

class Pidi : public MultiParticle
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
class Diquark : public MultiParticle
{
    public:
    
    int select_id  = -999;
    
    void diquarkReset()
    {
        v_id.clear();
        v_pid.clear();
        v_px.clear();
        v_py.clear();
        v_pz.clear();
        v_daughter.clear();
        v_parent.clear();
        v_mass.clear();
        v_vz.clear();
    }
};
class BinVariable
{
    
    public:
    //Need: x, z_h, Q2, pT, R0, R1, R2
    vector<double> v_x;
    vector<double> v_z_h;
    vector<double> v_Q2;
    vector<double> v_pT;
    vector<double> v_R0;
    vector<double> v_R1;
    vector<double> v_R2;
    
    double xmean;
    double z_hmean;
    double Q2mean;
    double pTmean;
    double R0mean;
    double R1mean;
    double R2mean;
    //zFillVectors(z_h, Q2, pT, R0, R1, R2)
    void zFillVectors(double x, double Q2, double pT, double R0, double R1, double R2) {
        v_x.push_back(x);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //xFillVectors(z_h, Q2, pT, R0, R1, R2);
    void xFillVectors(double z_h, double Q2, double pT, double R0, double R1, double R2) {
        v_z_h.push_back(z_h);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //mhFillVectors(x, z_h, Q2, pT, R0, R1, R2);
    void mhFillVectors(double x, double z_h, double Q2, double pT, double R0, double R1, double R2) {
        v_x.push_back(x);
        v_z_h.push_back(z_h);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    
    //Methods for calculating mean
    void meanZ_h() {
        xmean = meanfunc(v_x);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanx() {
        z_hmean = meanfunc(v_z_h);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanmh() {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
};
// 
//    Main body of analysis function
//

int LundAnalysis(
                 const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3051_0.hipo",
//                 const char * rootfile = "OutputFiles/AffinityFiles/Files_10_17/noRcuts4.root"
//                    const char * rootfile = "OutputFiles/Separate_Test_10_20/file2.root"
                   const char * rootfile = "OutputFiles/March_4/file1.root"
//                 const char * rootfile = "OutputFiles/AffinityFiles/Files_9_16/TMD1.root"
//                 const char * rootfile = "OutputFiles/AffinityFiles/Files_9_12/collinear1.root"
)
{
    gROOT->ProcessLine("#include <vector>");
    //Below file is now disappeared...
//    auto hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo";
// Current files: defined in main function though
//    auto hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3051_0.hipo";
//    auto rootFile = "OutputFiles/AffinityFiles/Files_9_5/Exactfile2.root";
    
    TFile *f = TFile::Open(rootfile,"RECREATE");
    
    HipoChain chain;
    
    //Add file to HipoChain
    chain.Add(hipoFile);
    auto config_c12 = chain.GetC12Reader();
    
    //Set PID cuts
    config_c12->addExactPid(11,1);    //exactly 1 electron
//    config_c12->addAtLeastPid(211,1);    //exactly 1 pi+
//    config_c12->addAtLeastPid(-211,1);    //exactly 1 pi-
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

    // Debugging counters
    int zbin0count = 0;
    int zbin1count = 0;
    int zbin2count = 0;
    int zbin3count = 0;
    int zbin4count = 0;
    int zbin5count = 0;
    int zbin6count = 0;
    
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
    double Q2_calc;
    double x;
    double pt_lab;
    double z_h;
    double z_h_plus;
    double z_h_minus;
    double xF;
    
    double m_plus;
    double m_minus;

    double q_TdivQ;
    double q_TdivQplus;
    double q_TdivQminus;
    
    //Cut Kinematics
    double W;
    double xFpiplus;
    double xFpiminus;
    double nu;
    double Mx;
    
    //Affinity ratios
    double R0;
    double R1;
    double R1_plus;
    double R1_minus;
    double R2;
    int qparent;
    int diparent;
    double s;
    double y;

    //Dihadron kinematics
    double Mdihadron;
    double dihadronpt;
    double Pdihadron;

    
    //MCParticle variables
    double P;
    double E;
    
    //Initializing particle vectors
    TLorentzVector q;
    TLorentzVector q_calc;
    double qx;
    double qy;
    double qz;
    double qE;
    
    double q_gNx;
    double q_gNy;
    double q_gNz;
    double q_gNE;
    
    double photonx;
    double photony;
    double photonz;
    double photonE;
    
    TLorentzVector init_electron;
    TLorentzVector init_target;
    TLorentzVector dihadron;
    
    TLorentzVector k;
    TLorentzVector kf;
    TLorentzVector ki;
    TVector2 deltak; // Transverse light cone vector - (V_x,V_y)
    TVector2 deltak_gN;
    
    double ki2;
    double kf2;
    double deltak2;
    int R0check;
    
    double kTx;
    double kTy;
    double kTz;
    double kTE;
    
    double kix;
    double kiy;
    double kiz;
    double kiE;

    double kfx;
    double kfy;
    double kfz;
    double kfE;
    
    //gNframe
    TLorentzVector lv_p1_gN;
    TLorentzVector lv_p2_gN;
    TLorentzVector dihadron_gN;
    
    TLorentzVector lv_q_gN;
    
    TLorentzVector target_gN;
    
    TLorentzVector ki_gN;
    TLorentzVector kf_gN;
    TLorentzVector k_gN;
    
    double pt_gN;
    
    double pt_gN_plus;
    double pt_gN_minus;
    
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
    TLorentzVector plusBreit;
    TLorentzVector minusBreit;
    TVector2 dihadronBreitTran;
    TVector2 plusBreitTran;
    TVector2 minusBreitTran;
    
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
    TLorentzVector plusPF;
    TLorentzVector minusPF;
    double dihadronPFMinus;
    double plusPFMinus;
    double minusPFMinus;
    
    double z_N;
    double z_Nplus;
    double z_Nminus;
    TVector2 q_T;
    TVector2 q_Tplus;
    TVector2 q_Tminus;
    
    double Qdiff;
        // Bin objects for collecting kinematic variables
    
    BinVariable zbin0;
    BinVariable zbin1;
    BinVariable zbin2;
    BinVariable zbin3;
    BinVariable zbin4;
    BinVariable zbin5;
    BinVariable zbin6;
    
    BinVariable xbin0;
    BinVariable xbin1;
    BinVariable xbin2;
    BinVariable xbin3;
    BinVariable xbin4;
    BinVariable xbin5;
    BinVariable xbin6;
    
    BinVariable Mhbin0;
    BinVariable Mhbin1;
    BinVariable Mhbin2;
    BinVariable Mhbin3;
    BinVariable Mhbin4;
    BinVariable Mhbin5;
    BinVariable Mhbin6;
    
    vector<double> xbins{0.1,0.13,0.16,0.19,0.235,0.3,0.5};
    vector<double> zbins{0.35,0.43,0.49,0.55,0.62,0.7,0.83};
    vector<double> Mhbins;
    for(int i = 0;i < 7; i++) {
        Mhbins.push_back(0.3 + i / 6.);
    }
    //Vectors for calculating means
    vector<BinVariable> zbinv = {zbin0, zbin1, zbin2, zbin3, zbin4, zbin5, zbin6};
    vector<BinVariable> xbinv = {xbin0, xbin1, xbin2, xbin3, xbin4, xbin5, xbin6};       
    vector<BinVariable> Mhbinv = {Mhbin0, Mhbin1, Mhbin2, Mhbin3, Mhbin4, Mhbin5, Mhbin6};
    vector<string> vinfoString = {"0th bin", "1st bin", "2nd bin", "3rd bin", "4th bin", "5th bin", "6th bin"};
    
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
    
    int hash_count = 0;

//     //Making new MC tree for piplus
//     TTree *t_plus = new TTree("tree_MC_plus","Tree with MC data from pi+ hadron");

//     t_plus->Branch("z",&z_h_plus);
//     t_plus->Branch("x",&x);
//     t_plus->Branch("pT",&pt_gN_plus);
//     t_plus->Branch("Q2",&Q2);
//     t_plus->Branch("R0max",&R0); //initial parton momentum
//     t_plus->Branch("R1max",&R1_plus); //final parton momentum
//     t_plus->Branch("R2max",&R2);
//     t_plus->Branch("Mh",&Mdihadron);
//     t_plus->Branch("q_TdivQ",&q_TdivQ);
    
    
//     //Making new MC tree for piminus
//     TTree *t_minus = new TTree("tree_MC_minus","Tree with MC data from pi- hadron");
//     t_minus->Branch("z",&z_h_minus);
//     t_minus->Branch("x",&x);
//     t_minus->Branch("pT",&pt_gN_minus);
//     t_minus->Branch("Q2",&Q2);
//     t_minus->Branch("R0max",&R0); //initial parton momentum
//     t_minus->Branch("R1max",&R1_minus); //final parton momentum
//     t_minus->Branch("R2max",&R2);
//     t_minus->Branch("Mh",&Mdihadron);
//     t_minus->Branch("q_TdivQ",&q_TdivQ);
    
        
    //Making new MC tree for dihadron
    TTree *tree_MC = new TTree("tree_MC","Tree with MC data from dihadron");
    tree_MC->Branch("z",&z_h);
    tree_MC->Branch("x",&x);
    tree_MC->Branch("pT",&pt_gN);
    tree_MC->Branch("Q2",&Q2);
    tree_MC->Branch("Q2calc",&Q2_calc);
    tree_MC->Branch("R0max",&R0); //initial parton momentum
    tree_MC->Branch("R1max",&R1); //final parton momentum
    tree_MC->Branch("R1maxplus",&R1_plus);
    tree_MC->Branch("R2max",&R2);
    tree_MC->Branch("Mh",&Mdihadron);
    tree_MC->Branch("q_TdivQ",&q_TdivQ);
    tree_MC->Branch("R0check", &R0check);
    
//     tree_MC->Branch("ki2",&ki2);
//     tree_MC->Branch("kf2",&kf2);
//     tree_MC->Branch("deltak2",&deltak2);
    
//     tree_MC->Branch("qx",&qx);
//     tree_MC->Branch("qy",&qy);
//     tree_MC->Branch("qz",&qz);
//     tree_MC->Branch("qE",&qE);
    
//     tree_MC->Branch("q_gNx",&q_gNx);
//     tree_MC->Branch("q_gNy",&q_gNy);
//     tree_MC->Branch("q_gNz",&q_gNz);
//     tree_MC->Branch("q_gNE",&q_gNE);
    
//     tree_MC->Branch("photonx",&photonx);
//     tree_MC->Branch("photony",&photony);
//     tree_MC->Branch("photonz",&photonz);
//     tree_MC->Branch("photonE",&photonE);
    
//     tree_MC->Branch("kix",&kix);
//     tree_MC->Branch("kiy",&kiy);
//     tree_MC->Branch("kiz",&kiz);
//     tree_MC->Branch("kiE",&kiE);
    
//     tree_MC->Branch("kfx",&kfx);
//     tree_MC->Branch("kfy",&kfy);
//     tree_MC->Branch("kfz",&kfz);
//     tree_MC->Branch("kfE",&kfE);
    
//     tree_MC->Branch("kTx",&kTx);
//     tree_MC->Branch("kTy",&kTy);
    
    
    
    //Tell the user that the loop is starting
    cout << "Start Event Loop" << endl;
    int tree_count = 0;
    //now get reference to (unique)ptr for accessing data in loop
    //this will point to the correct place when file changes
    //
    //This line comes from AnalysisWithExtraBanks.C
    auto& c12=chain.C12ref();
    int event_count = 0;
    
    //Loop over all events in the file
    while(chain.Next()==true){
        
        event_count += 1;
        if(event_count > 100000) {
            break;
        }
        if(event_count == 1) {
            cout << '\n';
            cout << "\033[96m";
            cout << "\t\t\t\t" << " ~~~~~~~~~~~~" << '\n';
            cout << "\t\t\t\t" << "|Progress Bar|" << '\n';
            cout << "\t\t\t\t" << " ~~~~~~~~~~~~" << '\n';
            cout << "\t\t[";
        }
        if(event_count % 16388 == 0) {
            
            hash_count += 1;
//            cout << "\033[A" << "\033[A";
            cout << '\r' << "\t[";
            for (int i = 1; i < hash_count + 1;i++) {
                cout << '#';
            }
            for (int i = 0; i < 100 - hash_count; i++) {
                cout << ' ';
            }
            cout << ']' << ' ';
            cout << event_count / 16388. << '%';
            if(event_count / 16388. == 100) {
                cout << endl;
            }
            cout << flush;
        }
        //Break at event 100 for testing with shorter run time
/*        if(event_count >= 1000) {
            cout << "Breaking at event: " << event_count << '\n';
            break;
        }*/
        if(c12->getDetParticles().empty())
            continue;
        
        //Intializing MCParticles
        MCParticle electron;
        MCParticle proton;
        MCParticle photon;
        MCParticle Lund;

        Pidi piplus;
        Pidi piminus;

        Quark quark;

        Pidi diquark;
        
        MultiParticle Hadron;
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
                electron.setVectors();
            }
            //pi+
            else if(pid==pipluspid){
                piplus.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                piplus.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
            //pi-
            else if(pid==piminuspid){
                piminus.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                piminus.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
            //all quarks
            else if(std::count(vquarklist.begin(), vquarklist.end(), pid)){
                quark.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                quark.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
            //MCParticle
            else if(pid==92 || pid == 91){
                Lund.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                Lund.setVectors();
            }
            
            else if(std::count(vdiquarklist.begin(), vdiquarklist.end(), pid)){
                diquark.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                diquark.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
            else if(pid == 22 && parent == 1){
                photon.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                photon.setVectors();
            }
            else if(id == 2){
                proton.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                proton.setVectors();
            }
            else if(std::count(vhadronlist.begin(), vhadronlist.end(), pid)) {
                Hadron.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                Hadron.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
        }
        
        //Selecting pions that come from Lund particle
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
        
        //Selecting initial quark
        for(int i = 0; i < quark.v_id.size(); i++) {
            if(quark.v_parent[i] == 0) {
                quark.final_id = i;
            }
        }
        
        //Selecting diquark
        for(int i = 0; i < diquark.v_id.size(); i++) {
            if(diquark.v_parent[i] == 2) {
                diquark.select_id = i;
            }
        }
        
        //Skip non-dipion events
        if(piplus.select_id == -999) {
            continue;
        }
        else if(piminus.select_id == -999) {
            continue;
        }
        if(diquark.select_id != -999) {
            diquark.fillParticle(diquark.v_id[diquark.select_id], diquark.v_pid[diquark.select_id], diquark.v_px[diquark.select_id], diquark.v_py[diquark.select_id], 
                           diquark.v_pz[diquark.select_id], diquark.v_daughter[diquark.select_id], diquark.v_parent[diquark.select_id], diquark.v_mass[diquark.select_id], diquark.v_vz[diquark.select_id]);
            diquark.setVectors();
        }
        piplus.fillParticle(piplus.v_id[piplus.select_id], piplus.v_pid[piplus.select_id], piplus.v_px[piplus.select_id], piplus.v_py[piplus.select_id], 
                           piplus.v_pz[piplus.select_id], piplus.v_daughter[piplus.select_id], piplus.v_parent[piplus.select_id], piplus.v_mass[piplus.select_id], piplus.v_vz[piplus.select_id]);
        piplus.setVectors();
        piminus.fillParticle(piminus.v_id[piminus.select_id], piminus.v_pid[piminus.select_id], piminus.v_px[piminus.select_id], piminus.v_py[piminus.select_id], 
                           piminus.v_pz[piminus.select_id], piminus.v_daughter[piminus.select_id], piminus.v_parent[piminus.select_id], piminus.v_mass[piminus.select_id], piminus.v_vz[piminus.select_id]);
        piminus.setVectors();
        
        quark.fillParticle(quark.v_id[quark.final_id], quark.v_pid[quark.final_id], quark.v_px[quark.final_id], quark.v_py[quark.final_id], 
                           quark.v_pz[quark.final_id], quark.v_daughter[quark.final_id], quark.v_parent[quark.final_id], quark.v_mass[quark.final_id], quark.v_vz[quark.final_id]);
        quark.setVectors();
       
        //Setting inital beam and target particles
        init_electron.SetPxPyPzE(0, 0, sqrt(electron_beam_energy * electron_beam_energy - electronMass * electronMass), electron_beam_energy);
        init_target.SetPxPyPzE(0, 0, 0, proton.E);
        
        
        dihadron = piplus.lv + piminus.lv;
        m_plus = piplus.mass;
        m_minus = piminus.mass;
        Mdihadron = dihadron.M();
        Pdihadron = dihadron.P();
//         q = init_electron - electron.lv; //virtual photon
        q = photon.lv;
        q_calc = init_electron - electron.lv;
        qx = q_calc.Px();
        qy = q_calc.Py();
        qz = q_calc.Pz();
        qE = q_calc.E();
        
        photonx = photon.lv.Px();
        photony = photon.lv.Py();
        photonz = photon.lv.Pz();
        photonE = photon.lv.E();
        
        
        //Missing mass
        Mx = Mxfunc(q, init_target, piplus.lv, piminus.lv);

        cth = cthfunc(electron.px,electron.py,electron.pz);
        Q2 = -(q * q);
        Q2_calc = Q2func(electron_beam_energy,electron.E,cth); //Momentum transfer

        z_h_plus = (init_target * piplus.lv) / (init_target * q);
        z_h_minus = (init_target * piminus.lv) / (init_target * q);
        z_h = z_h_plus + z_h_minus;
        s = sfunc(protonMass, electronMass, electron_beam_energy);
        y = yfunc(electron_beam_energy,electron.E);
        x = Q2/s/y; // Bjorken x
        pt_lab = Ptfunc(dihadron.Px(), dihadron.Py()); //hadron transverse momentum
        
        kf = quark.lv;
        ki = kf - q;
        //Cut Kinematics
        
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
        
        q_gNx = lv_q_gN.Px();
        q_gNy = lv_q_gN.Py();
        q_gNz = lv_q_gN.Pz();
        q_gNE = lv_q_gN.E();
        
        //Need dihadron in gN frame for pT
        dihadron_gN = dihadron;
        dihadron_gN.Boost(gNBoostNeg);
        pt_gN = Ptfunc(dihadron_gN);
        pt_gN_plus = Ptfunc(lv_p1_gN);
        pt_gN_minus = Ptfunc(lv_p2_gN);
        
        //Need target in gN
        target_gN = init_target;
        target_gN.Boost(gNBoostNeg);
        
        //Need partonic in gN
        ki_gN = ki;
        ki_gN.Boost(gNBoostNeg);
        
        kf_gN = kf;
        kf_gN.Boost(gNBoostNeg);
        
        //Feynman x
        xFpiplus = xFfunc(lv_p1_gN,lv_q_gN,W);
        xFpiminus = xFfunc(lv_p2_gN,lv_q_gN,W);
        
        //nu and W
        nu = nufunc(electron_beam_energy,electron.E);
        W = Wfunc(Q2,protonMass,nu);
        
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
        
        plusBreit = piplus.lv;
        plusBreit.Boost(BreitBoost);
        plusBreitTran = PtVectfunc(plusBreit);
        minusBreit = piminus.lv;
        minusBreit.Boost(BreitBoost);
        minusBreitTran = PtVectfunc(minusBreit);
        
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
        
        plusPF = piplus.lv;
        plusPF.Boost(PFBoost);
        plusPF.Rotate(PFAngle,PFAxis);
        plusPFMinus = LightConeMinus(plusPF);
        minusPF = piminus.lv;
        minusPF.Boost(PFBoost);
        minusPF.Rotate(PFAngle,PFAxis);
        minusPFMinus = LightConeMinus(minusPF);
        //Virtual Photon
        qPF.Rotate(PFAngle,PFAxis);
        qPFMinus = LightConeMinus(qPF);
        //z_N and q_T
        z_N = dihadronPFMinus / qPFMinus;
        z_Nplus = plusPFMinus / qPFMinus;
        z_Nminus = minusPFMinus / qPFMinus;
        q_T = -1 * dihadronBreitTran / z_N;
        q_Tplus = -1 * plusBreitTran / z_Nplus;
        q_Tminus = -1 * minusBreitTran / z_Nminus;
        

        //q_T / Q for plotting
        q_TdivQ = Ptfunc(q_T) / sqrt(Q2);
        q_TdivQplus = Ptfunc(q_Tplus) / sqrt(Q2);
        q_TdivQminus = Ptfunc(q_Tminus) / sqrt(Q2);
        

        //ki, k, and delta k
        deltak = kfBreitTran - (-1 * z_N * q_T); 
        
        k = kf - q;
        k_gN = k;
        k_gN.Boost(gNBoostNeg);
//         These ratios are calculated in lab frame
//        R0 = R0func(ki, kf, deltak, Q2);
//        R1 = R1func(dihadron, ki, kf);
//        R2 = R2func(k, Q2);
        
        //Ratios in gN frame
        R0 = R0func(ki_gN, kf_gN, deltak, Q2);
        kix = ki_gN.Px();
        kiy = ki_gN.Py();
        kiz = ki_gN.Pz();
        kiE = ki_gN.E();
        
        kfx = kf_gN.Px();
        kfy = kf_gN.Py();
        kfz = kf_gN.Pz();
        kfE = kf_gN.E();
        
        kTx = deltak_gN.Px();
        kTy = deltak_gN.Py();
        double ki2 = abs(ki_gN * ki_gN);
        double kf2 = abs(kf_gN * kf_gN);
        double deltak2 = abs(deltak * deltak);
        if(deltak2 > ki2 && deltak2 > kf2) {
            R0check = 0;//DeltaK is biggest
        }
        else if(ki2 > kf2) {
            R0check = 1;//ki is biggest
        }
        else {
            R0check = 2;//kf is biggest
        }
        
        R1 = R1func(dihadron_gN, ki_gN, kf_gN);
        R1_plus = R1func(lv_p1_gN,ki_gN,kf_gN);
        R1_minus = R1func(lv_p2_gN,ki_gN,kf_gN);
        
        R2 = R2func(k_gN, Q2);
        xF = xFpiplus + xFpiminus;
        
        //CUTS:
        //Region cuts:
        
/*
        //TMD and Collinear region selection:
        if(R0 > 0.3 || R1 > 0.3) {
            continue;
        }
        //TMD:
        if(R2 > 0.3) {
            continue;
        }
        //Collinear
//        if(R2 < 0.9) {
//            continue;
//        }
         
*/
        //Missing mass
        if(Mx <= 1.5) {
            continue;
        }
        
        //Feynman x
        if(xFpiplus <= 0 || xFpiminus <= 0) {
            continue;
        }
        
        //Vertex Position
        if(abs(electron.vz - piplus.vz) >= 20) {
            continue;
        }
        if(abs(electron.vz - piminus.vz) >= 20) {
            continue;
        }
        if(electron.vz <= -8 || electron.vz >= 3) {
            continue;
        }
        
        //Channel Selection
        
        //Virtual photon mass / momentum transfer
        if(Q2 <= 1 || Q2 >= 100) {
            continue;
        }
        //Hadronic system mass
        if(W <= 2 || W >= 100) {
            continue;
        }
        //Energy fraction
        if(y <= 0 || y >= 0.8) {
            continue;
        }
        if(piplus.P <= 1.25 || piminus.P <= 1.25) {
            continue;
        }

        tree_count += 1;
        tree_MC->Fill();
//         t_plus->Fill();
//         t_minus->Fill();
        //Need: x, z, Q2, pT, R0, R1, R2
        //zbins:
	
        for(int i = 0; i < zbins.size(); i++) {
            if(z_h <= zbins[i]) {
                zbinv[i].zFillVectors(x, Q2, pt_gN, R0, R1, R2);
                break;
            }
        }
        //Mh bins
        for(int i = 0; i < Mhbins.size(); i++) {
            if(Mdihadron <= Mhbins[i]) {
                Mhbinv[i].mhFillVectors(x, z_h, Q2, pt_gN, R0, R1, R2);
                break;
            }
        }
        //x bins
        for(int i = 0; i < xbins.size(); i++) {
            if(x <= xbins[i]) {
                xbinv[i].xFillVectors(z_h, Q2, pt_gN, R0, R1, R2);
                break;
            }
        }
        //print out tree count every 100 to give update to user
        if(tree_count % 100 == 0) {
	//            cout << "Tree_count: " << tree_count << '\n';
        }
	
    }
    cout << "\033[0m" << "\033[49m";
    cout << "Final tree_count: " << tree_count << '\n';
    
    //Making new Affinity trees
    TTree *t_z_h = new TTree("tree_z_h_bins","Tree with mean values binned by z_h affinity calculations");
    TTree *t_x = new TTree("tree_x_bins","Tree with mean values binned by x affinity calculations");
    TTree *t_Mh = new TTree("tree_Mh_bins","Tree with mean values binned by Mh affinity calculations");
    
    string infoString;
    Double_t z_h_t;
    Double_t x_t;
    Double_t Q2_t;
    Double_t pT_t;
    Double_t R0_t;
    Double_t R1_t;
    Double_t R2_t;
    
    t_z_h->Branch("Name",&infoString);
    t_z_h->Branch("x", &x_t);
    t_z_h->Branch("Q2", &Q2_t);
    t_z_h->Branch("pT", &pT_t);
    t_z_h->Branch("R0", &R0_t);
    t_z_h->Branch("R1", &R1_t);
    t_z_h->Branch("R2", &R2_t);
    
    t_x->Branch("Name",&infoString);
    t_x->Branch("z_h", &z_h_t);
    t_x->Branch("Q2", &Q2_t);
    t_x->Branch("pT", &pT_t);
    t_x->Branch("R0", &R0_t);
    t_x->Branch("R1", &R1_t);
    t_x->Branch("R2", &R2_t);
    
    t_Mh->Branch("Name",&infoString);
    t_Mh->Branch("x", &x_t);
    t_Mh->Branch("z_h", &z_h_t);
    t_Mh->Branch("Q2", &Q2_t);
    t_Mh->Branch("pT", &pT_t);
    t_Mh->Branch("R0", &R0_t);
    t_Mh->Branch("R1", &R1_t);
    t_Mh->Branch("R2", &R2_t);
    
    //Calculating means
    //Setting zbin means
    for(int i = 0; i < vinfoString.size(); i++) {
        zbinv[i].meanZ_h();
        infoString = vinfoString[i];
        x_t = zbinv[i].xmean;
        Q2_t = zbinv[i].Q2mean;
        pT_t = zbinv[i].pTmean;
        R0_t = zbinv[i].R0mean;
        R1_t = zbinv[i].R1mean;
        R2_t = zbinv[i].R2mean;
        t_z_h->Fill();
        }
    for(int i = 0; i < vinfoString.size(); i++) {
        xbinv[i].meanx();
        infoString = vinfoString[i];
        z_h_t = xbinv[i].z_hmean;
        Q2_t = xbinv[i].Q2mean;
        pT_t = xbinv[i].pTmean;
        R0_t = xbinv[i].R0mean;
        R1_t = xbinv[i].R1mean;
        R2_t = xbinv[i].R2mean;
        t_x->Fill();
        }
    for(int i = 0; i < vinfoString.size(); i++) {
        Mhbinv[i].meanmh();
        infoString = vinfoString[i];
        x_t = Mhbinv[i].xmean;
        z_h_t = Mhbinv[i].z_hmean;
        Q2_t = Mhbinv[i].Q2mean;
        pT_t = Mhbinv[i].pTmean;
        R0_t = Mhbinv[i].R0mean;
        R1_t = Mhbinv[i].R1mean;
        R2_t = Mhbinv[i].R2mean;
        t_Mh->Fill();
        }
    
    f->Write();
    delete f;
    
    return 0;
}
