#include <vector>
#include <cmath>
#include "TROOT.h"
#pragma once
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
        //qTQFillVectors(x, z_h, Q2, pT, R0, R1, R2);
    void qTQFillVectors(double x, double z_h, double Q2, double pT, double R0, double R1, double R2) {
        v_x.push_back(x);
        v_z_h.push_back(z_h);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
        //Q2FillVectors(z_h, Q2, pT, R0, R1, R2);
    void Q2FillVectors(double x, double z_h, double pT, double R0, double R1, double R2) {
        v_z_h.push_back(z_h);
        v_x.push_back(x);
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
    void meanqTQ() {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanQ2() {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
};