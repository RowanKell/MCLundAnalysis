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
    
    vector<double> v_x_1;
    vector<double> v_z_h_1;
    vector<double> v_Q2_1;
    vector<double> v_pT_1;
    
    vector<double> v_x_2;
    vector<double> v_z_h_2;
    vector<double> v_Q2_2;
    vector<double> v_pT_2;
    vector<double> v_R0;
    vector<double> v_R1;
    vector<double> v_R1_p;
    vector<double> v_R1_m;
    vector<double> v_R2;
    
    double xmean;
    double z_hmean;
    double Q2mean;
    double pTmean;
    
    double z_hmean_1;
    double pTmean_1;
    
    double z_hmean_2;
    double pTmean_2;
    double R0mean;
    double R1mean;
    double R1_p_mean;
    double R1_m_mean;
    double R2mean;
    
    /*
    First set of functions is for dihadron direct affinity box calculations
    
    Second set of functions overload the first and are for storing the R1 for each hadron as well
    */
    //zFillVectors(z_h, Q2, pT, R0, R1, R2)
    void zFillVectors(double x, double Q2, double pT, double z_h_1, double pT_1, double z_h_2, double pT_2, double R0, double R1, double R2) {
        v_x.push_back(x);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);

        v_z_h_1.push_back(z_h_1);
        v_pT_1.push_back(pT_1);

        v_z_h_2.push_back(z_h_2);
        v_pT_2.push_back(pT_2);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //Overload: z (use this for single_pion where only 1 z, pT, etc)
    void zFillVectors(double x, double Q2, double pT, double R0, double R1, double R2) {
        v_x.push_back(x);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //xFillVectors(z_h, Q2, pT, R0, R1, R2);
    void xFillVectors(double z_h, double Q2, double pT,double z_h_1, double pT_1,double z_h_2, double pT_2, double R0, double R1, double R2) {
        v_z_h.push_back(z_h);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);
        
        v_z_h_1.push_back(z_h_1);
        v_pT_1.push_back(pT_1);
        
        v_z_h_2.push_back(z_h_2);
        v_pT_2.push_back(pT_2);

        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //Overload: x
    void xFillVectors(double z_h, double Q2, double pT, double R0, double R1, double R2) {
        v_z_h.push_back(z_h);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //mhFillVectors(x, z_h, Q2, pT, R0, R1, R2);
    void mhFillVectors(double x, double z_h, double Q2, double pT, double z_h_1, double pT_1, double z_h_2, double pT_2, double R0, double R1, double R2) {
        v_x.push_back(x);
        v_z_h.push_back(z_h);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);
        

        v_z_h_1.push_back(z_h_1);
        v_pT_1.push_back(pT_1);
        

        v_z_h_2.push_back(z_h_2);
        v_pT_2.push_back(pT_2);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //overload: Mh
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
    void qTQFillVectors(double x, double z_h, double Q2, double pT, double z_h_1, double pT_1, double z_h_2, double pT_2, double R0, double R1, double R2) {
        v_x.push_back(x);
        v_z_h.push_back(z_h);
        v_Q2.push_back(Q2);
        v_pT.push_back(pT);
        

        v_z_h_1.push_back(z_h_1);
        v_pT_1.push_back(pT_1);
        

        v_z_h_2.push_back(z_h_2);
        v_pT_2.push_back(pT_2);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //Overload: qTQ
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
    void Q2FillVectors(double x, double z_h, double pT, double z_h_1, double pT_1, double z_h_2, double pT_2, double R0, double R1, double R2) {
        v_x.push_back(x);
        v_z_h.push_back(z_h);
        v_pT.push_back(pT);
        
        v_z_h_1.push_back(z_h_1);
        v_pT_1.push_back(pT_1);
        
        v_z_h_2.push_back(z_h_2);
        v_pT_2.push_back(pT_2);
        
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //Overload: Q2
    void Q2FillVectors(double x, double z_h, double pT, double R0, double R1, double R2) {
        v_x.push_back(x);
        v_z_h.push_back(z_h);
        v_pT.push_back(pT);
        v_R0.push_back(R0);
        v_R1.push_back(R1);
        v_R2.push_back(R2);
    }
    //Methods for calculating mean
    void meanZ_h() {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        
        z_hmean_1 = meanfunc(v_z_h_1);
        pTmean_1 = meanfunc(v_pT_1);
        
        z_hmean_2 = meanfunc(v_z_h_2);
        pTmean_2 = meanfunc(v_pT_2);
        
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanx() {
        z_hmean = meanfunc(v_z_h);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        
        z_hmean_1 = meanfunc(v_z_h_1);
        pTmean_1 = meanfunc(v_pT_1);
        
        z_hmean_2 = meanfunc(v_z_h_2);
        pTmean_2 = meanfunc(v_pT_2);
        
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanmh() {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        
        z_hmean_1 = meanfunc(v_z_h_1);
        pTmean_1 = meanfunc(v_pT_1);
        
        z_hmean_2 = meanfunc(v_z_h_2);
        pTmean_2 = meanfunc(v_pT_2);
        
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanqTQ() {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        
        z_hmean_1 = meanfunc(v_z_h_1);
        pTmean_1 = meanfunc(v_pT_1);
        
        z_hmean_2 = meanfunc(v_z_h_2);
        pTmean_2 = meanfunc(v_pT_2);
        
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanQ2() {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        pTmean = meanfunc(v_pT);
        
        z_hmean_1 = meanfunc(v_z_h_1);
        pTmean_1 = meanfunc(v_pT_1);
        
        z_hmean_2 = meanfunc(v_z_h_2);
        pTmean_2 = meanfunc(v_pT_2);
        
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    
// OVERLOADED mean FUNCTIONS FOR SINGLE_PION: quick fix - use a single argument in function to use this definition which doesn't involve z_h_1 etc
    void meanZ_h(int single_pion) {
        xmean = meanfunc(v_x);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanx(int single_pion) {
        z_hmean = meanfunc(v_z_h);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanmh(int single_pion) {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanqTQ(int single_pion) {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        Q2mean = meanfunc(v_Q2);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
    void meanQ2(int single_pion) {
        z_hmean = meanfunc(v_z_h);
        xmean = meanfunc(v_x);
        pTmean = meanfunc(v_pT);
        R0mean = meanfunc(v_R0);
        R1mean = meanfunc(v_R1);
        R2mean = meanfunc(v_R2);
    }
};