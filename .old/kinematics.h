#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <iostream>
#include <vector>
#include <cmath>
#include "TROOT.h"

double Pfunc(double Px, double Py, double Pz);

double Efunc(double M, double P);

double EVirtualfun(double M, double P);
    
double Ptfunc(double Px, double Py);
    
double Ptfunc(TLorentzVector lv);
double Ptfunc(TVector2 v);

TVector2 PtVectfunc(TLorentzVector lv);

double cthfunc(double Px, double Py, double Pz);

double Q2func(double E1, double E2, double cth);

double yfunc(double E1,double E2);

double sfunc(double M1,double M2,double E); //proton mass is M1, electron mass is M2, E is electron beam energy

double R0func(TLorentzVector ki, TLorentzVector kf, TVector2 deltak, double Q2);

double R1func(TLorentzVector Ph, TLorentzVector ki, TLorentzVector kf);

double R2func(TLorentzVector k, double Q2);

double Mxfunc(TLorentzVector q, TLorentzVector target, TLorentzVector hadron1, TLorentzVector hadron2);

double xFfunc(TLorentzVector p, TLorentzVector q, double W);

double nufunc(double E1, double E2);

double Wfunc(double Q2, double mT, double nu);

double thetafunc(double pt, double pz);

double LightConeMinus(TLorentzVector lv);

double LightConePlus(TLorentzVector lv);

double meanfunc(vector<double> v);

#endif