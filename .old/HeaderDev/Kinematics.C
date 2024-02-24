#include "Kinematics.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
using namespace std;

double Kinematics::Pfunc(double Px, double Py, double Pz)
{
    return sqrt(Px*Px + Py*Py + Pz*Pz);
}

double Kinematics::Efunc(double M, double P)
{
    return sqrt(M * M + P * P);
}

double Kinematics::Ptfunc(double Px, double Py)
{
    return sqrt(Px*Px + Py*Py);
}
double Kinematics::Ptfunc(TLorentzVector lv) {
    double x = lv.Px();
    double y = lv.Py();
    return sqrt(x * x + y * y);
}

TVector2 Kinematics::PtVectfunc(TLorentzVector lv)
{
    TVector2 Pt;
    double Px = lv.Px();
    double Py = lv.Py();
    Pt.SetX(Px);
    Pt.SetY(Py);
    return Pt;
}

double Kinematics::cthfunc(double Px, double Py, double Pz)
{
    double Pt = sqrt(Px*Px + Py*Py);
    return Pz / sqrt(Pz * Pz + Pt * Pt);
}

double Kinematics::Q2func(double E1, double E2, double cth)
{
    return 2 * E1 * E2 * (1.0 - cth);
}

double Kinematics::yfunc(double E1,double E2)
{
    return (E1 - E2) / E1;
}
double Kinematics::sfunc(double M1,double M2,double E) //proton mass is M1, electron mass is M2, E is electron beam energy
{
    return M1 * M1 + M2 * M2 + 2 * M1 * E;
}

double Kinematics::R0func(TLorentzVector ki, TLorentzVector kf, TVector2 deltak, double Q2)
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

double Kinematics::R1func(TLorentzVector Ph, TLorentzVector ki, TLorentzVector kf)
{
    return (Ph * ki) / (Ph * kf);
}

double Kinematics::R2func(TLorentzVector k, double Q2)
{
    return abs(k * k) / Q2;
}

double Kinematics::Mxfunc(TLorentzVector q, TLorentzVector target, TLorentzVector hadron1, TLorentzVector hadron2)
{
    TLorentzVector lv_Mx;
    lv_Mx = q + target - hadron1 - hadron2;
    return lv_Mx.M();
}
double Kinematics::xFfunc(TLorentzVector p, TLorentzVector q, double W)
{
    return 2 * (p.Vect().Dot(q.Vect())) / (q.Vect().Mag() * W);
}

double Kinematics::nufunc(double E1, double E2)
{
  return E1-E2;
}

double Kinematics::Wfunc(double Q2, double mT, double nu)
{
    return sqrt(-Q2 + pow(mT,2) + 2 * mT * nu);
}

double Kinematics::thetafunc(double pt, double pz)
{
    return abs(atan(pt/pz));
}

double Kinematics::LightConeMinus(TLorentzVector lv)
{
    return (lv.E() - lv.Pz()) / (sqrt(2));
}

double Kinematics::LightConePlus(TLorentzVector lv)
{
    return (lv.E() + lv.Pz()) / (sqrt(2));
}

double Kinematics::meanfunc(vector<double> v)
{
    double sum = 0;
    for(int i = 0; i < v.size(); i++) {
        sum += v[i];
    }
    return sum / v.size();
}