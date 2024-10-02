#include <vector>
#include <cmath>
#include "TROOT.h"
#pragma once

//3-Momentum magnitude
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

//Transverse momentum magnitude
double Ptfunc(double Px, double Py)
{
    return sqrt(Px*Px + Py*Py);
}

//""
double Ptfunc(TLorentzVector lv) {
    double x = lv.Px();
    double y = lv.Py();
    return sqrt(x * x + y * y);
}
//"" (these three just overload the same func)
double Ptfunc(TVector2 v) {
    double x = v.X();
    double y = v.Y();
    return sqrt(x * x + y * y);
}

//Transverse momentum 2-vector
TVector2 PtVectfunc(TLorentzVector lv)
{
    TVector2 Pt;
    double Px = lv.Px();
    double Py = lv.Py();
    Pt.SetX(Px);
    Pt.SetY(Py);
    return Pt;
}

//Cosine of theta (polar angle)
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
double xfunc(double Q2, double s, double y){
  return Q2/(s-0.938272*0.938272)/y;
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
    return (Ph * kf) / (Ph * ki);
}

double R2func(TLorentzVector k, double Q2)
{
    return abs(k * k) / Q2;
}

double R2func_adjust(double qT, double Q2)
{
    return (qT * qT) / Q2;
}

double Mxfunc(TLorentzVector q, TLorentzVector target, TLorentzVector hadron1, TLorentzVector hadron2)
{
    TLorentzVector lv_Mx;
    lv_Mx = q + target - hadron1 - hadron2;
    return lv_Mx.M();
}

double Mxfunc(TLorentzVector q, TLorentzVector target, TLorentzVector hadron1)
{
    TLorentzVector lv_Mx;
    lv_Mx = q + target - hadron1;
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

double xNfunc(double x_Bj, double M, double Q2) {
    double num = 2 * x_Bj;
    double denom = 1 + pow((1 + (4 * pow(x_Bj,2)*pow(M,2)) / Q2),0.5);
    return num / denom;
}