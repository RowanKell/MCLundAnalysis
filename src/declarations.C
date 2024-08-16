#include "BinVariable.C"
#include <random>
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
double R2_adjust;


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

double M_ki;
double M_kf;
double delta_k_T;
double ki_T;

double random_mean = 0.93891897 / 2;
double random_std = 0.93891897 / 2;

double r_M_ki;
double r_M_kf;
double r_delta_k_T;

std::default_random_engine generator_M_ki;
std::normal_distribution<double> distribution_M_ki(random_mean,random_std);

std::default_random_engine generator_M_kf;
std::normal_distribution<double> distribution_M_kf(random_mean,random_std);


std::default_random_engine generator_delta_k_T;
std::normal_distribution<double> distribution_delta_k_T(random_mean,random_std);

int R0check;


//


TLorentzVector lundstring;


//
//Breit frame variables
//

//SIDIS
TLorentzVector q_BF;
TLorentzVector proton_BF;
TLorentzVector pi1_BF;

//pT stuff
TVector2 pT_BF_vect;
double pT_BF;

//partons
TLorentzVector ki_BF;
TLorentzVector kf_BF;
TLorentzVector k_BF;
TVector2 kfT_BF_vect;
double kfT_BF;

//Boosting
TLorentzVector BF;
TLorentzVector BF_target;
TVector3 BFBoost;

//rotation
TLorentzVector q_BF_no_rot;
TVector3 q_BF_no_rot_VectUnit;
double BFAngle;
TVector3 BFAxis;

//misc
TLorentzVector proton_BF_no_rot;


//
//Hadron Frame Variables
//

//SIDIS
TLorentzVector pi1_HF;
TLorentzVector q_HF;

double pT_HF;

//qT stuff
double qT_HF;
double qTQ_HF;

//partons
TLorentzVector k_HF;
TLorentzVector ki_HF;
TLorentzVector kf_HF;

//Boosting
TLorentzVector Hadron_frame;
TVector3 HadronBoost;

//Rotation
TLorentzVector pi1_HF_no_rot;
TVector3 pi1_HF_no_rot_VectUnit;
double HFAngle;
TVector3 HFAxis;

//
//Photon Frame variables
//

//SIDIS Stuff
TLorentzVector P_PF;
TLorentzVector q_PF;
TLorentzVector pi1_PF;
double pT_PF;

//Boost
TLorentzVector PFFrame;
TVector3 PFBoost;

//Rotation
TVector3 qPFVectUnit;
TVector3 zAxis(0,0,1);
double PFAngle;
TVector3 PFAxis;
TLorentzVector q_PF_no_rot;
TVector3 q_PF_no_rot_VectUnit;

//Partons
TLorentzVector ki_PF;
TLorentzVector kf_PF;
TLorentzVector k_PF;

//
//other qTs
//

TVector2 qT_from_zN_vect;
double qT_from_zN;

TVector2 qT_from_pT_BF_vect;
double qT_from_pT_BF;

TVector2 qT_from_pT_PF_z_vect;
double qT_from_pT_PF_z;

TVector2 qT_from_pT_HF_z_vect;
double qT_from_pT_HF_z;

double qTQ_from_pT_BF;




double qPFMinus;
TLorentzVector dihadronPF;
TLorentzVector PF1;
TLorentzVector PF2;
double dihadronPFMinus;
double PF1Minus;
double PF2Minus;

double z_N;
double x_N;


//qT variables for single pion affinity qTQ binning
TLorentzVector pion_frame;
TVector3 pionBoost;
TLorentzVector q_pion;


TVector2 q_T_frac;

double Qdiff;

//Variables for individual qTQ bins
double qTQ_low_aff;
double qTQ_high_aff;

double qTQ_low_tol;
double qTQ_high_tol;

//qT from zN
double q_T_zN_val;

//sidis regions stuff
double xi;
double zeta;

double theta_ki;
double theta_H;
double theta_deltak;

double M;

double x_N_kin;
TLorentzVector PPF;
double P_B_b_minus;
double q_b_minus;

double qPF_T;
double PPF_T;

//One binnings
int low_bool;
int high_bool;

double Q2_low_aff;
double Q2_low_tol;
double pT_low_aff;
double pT_low_tol;
double x_low_aff;
double x_low_tol;
double z_low_aff;
double z_low_tol;

double Q2_high_aff;
double Q2_high_tol;
double pT_high_aff;
double pT_high_tol;
double x_high_aff;
double x_high_tol;
double z_high_aff;
double z_high_tol;

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
vector<double> qTQbins{0.1,0.3,0.5,0.8,1.1,1.5,2,2.5,3};
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

double R0_min = 9999;
double R1_min = 9999;
double R2_min = 9999;

double R0_max = 0;
double R1_max = 0;
double R2_max = 0;
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