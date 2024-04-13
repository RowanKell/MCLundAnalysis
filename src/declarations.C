#include "BinVariable.C"
#pragma once
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
double z_h_1;
double z_h_2;
double xF;

double m_1;
double m_2;

double qTQ;
double qTQ_lab;
double qTQ_hadron;
double q_TdivQ1;
double q_TdivQ2;

//Cut Kinematics
double W;
double xFpi1;
double xFpi2;
double nu;
double Mx;

//Affinity ratios
double R0;
double R1;
double R1_p;
double R1_m;
double R1lab;
double R1breit;
double R2;

//Numerator and Denominator for R1
double R1num;
double R1denom;
double R1num_lab;
double R1denom_lab;
double R1num_Breit;
double R1denom_Breit;

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

TLorentzVector ki_Breit;
TLorentzVector kf_Breit;
TLorentzVector k_Breit;
TLorentzVector q_Breit;
TLorentzVector proton_Breit;

double pt_gN;

double pt_gN_1;
double pt_gN_2;

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
TLorentzVector Breit1;
TLorentzVector Breit2;
TVector2 dihadronBreitTran;
TVector2 BreitTran1;
TVector2 BreitTran2;

//Hadron Frame Variables
TLorentzVector q_hadron;
TLorentzVector Hadron_frame;
TVector3 HadronBoost;

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
TLorentzVector PF1;
TLorentzVector PF2;
double dihadronPFMinus;
double PF1Minus;
double PF2Minus;

double z_N;
double z_N1;
double z_N2;
TVector2 q_T;
TVector2 q_T1;
TVector2 q_T2;
TVector2 q_T_lab;
TVector2 q_T_hadron;

TVector2 q_T_frac;
double qTQfrac;

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

BinVariable Q2bin0;
BinVariable Q2bin1;
BinVariable Q2bin2;
BinVariable Q2bin3;
BinVariable Q2bin4;
BinVariable Q2bin5;
BinVariable Q2bin6;
BinVariable Q2bin7;

BinVariable qTQbin0;
BinVariable qTQbin1;
BinVariable qTQbin2;
BinVariable qTQbin3;
BinVariable qTQbin4;
BinVariable qTQbin5;
BinVariable qTQbin6;
BinVariable qTQbin7;
BinVariable qTQbin8;

vector<double> xbins{0.1,0.13,0.16,0.19,0.235,0.3,0.5};
vector<double> zbins{0.35,0.43,0.49,0.55,0.62,0.7,0.83};
double z_h_cut_val;
int j_start; //Used to either run over duplicate pion pairs or to not run over duplicates
vector<double> Q2bins{1,1.4,2,2.8,4,5.6,7.9,11.1};
vector<double> qTQbins{0.1,0.3,0.5,0.8,1.5,2,2.5,3,4};
vector<double> Mhbins{0.4429,0.6357,0.8286,1.0214,1.2143,1.4071,1.6};


//Vectors for calculating means
vector<BinVariable> zbinv = {zbin0, zbin1, zbin2, zbin3, zbin4, zbin5, zbin6};
vector<BinVariable> xbinv = {xbin0, xbin1, xbin2, xbin3, xbin4, xbin5, xbin6};       
vector<BinVariable> Mhbinv = {Mhbin0, Mhbin1, Mhbin2, Mhbin3, Mhbin4, Mhbin5, Mhbin6};
vector<BinVariable> Q2binv = {Q2bin0, Q2bin1, Q2bin2, Q2bin3, Q2bin4, Q2bin5, Q2bin6, Q2bin7};
vector<BinVariable> qTQbinv = {qTQbin0, qTQbin1, qTQbin2, qTQbin3, qTQbin4, qTQbin5, qTQbin6,qTQbin7,qTQbin8};
vector<string> vinfoString = {"0th bin", "1st bin", "2nd bin", "3rd bin", "4th bin", "5th bin", "6th bin", "7th bin", "8th bin"};

//Creating vectors to fill using push_back in loop

vector<int> vquarklist = {-8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8};
vector<int> vdiquarklist = {1103, 2101, 2103, 2203, 3101, 3103, 3201, 3203, 3303, 4101, 4103, 4201, 4203, 4301, 4303, 4403, 5101, 5103, 5201, 5203, 5301, 5303, 5401, 5403, 5503};
vector<int> vhadronlist = {-3122, -211, 111, 211, 1114, 2114, 2212, 2214, 2224, 3112, 3114, 3122, 3214, 3222, 3224, 3312, 3324, -323, -313, -213, 113, 213, 221, 223, 310, 313, 323, 331, 333};

int hash_count = 0;
int declarations() {
//     cout << "qTQbins:\n";
//     for(int i = 0;i < 9; i++) {
//         cout << i << "th bin: " <<  qTQbins[i] << "\n";
//     }
//     for(int i = 0;i < 7; i++) {
//         Mhbins.push_back(0.25 + (i+1) / 6.);
//     }
    return 0;
}