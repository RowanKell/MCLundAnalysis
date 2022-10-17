#include "BinVariable.h"
//zFillVectors(z_h, Q2, pT, R0, R1, R2)
BinVariable::BinVariable() {
    xmean = 0;
    z_hmean = 0;
    Q2mean = 0;
    pTmean = 0;
    R0mean = 0;
    R1mean = 0;
    R2mean = 0;
}
BinVariable::zFillVectors(double x, double Q2, double pT, double R0, double R1, double R2) {
    v_x.push_back(x);
    v_Q2.push_back(Q2);
    v_pT.push_back(pT);
    v_R0.push_back(R0);
    v_R1.push_back(R1);
    v_R2.push_back(R2);
}
//xFillVectors(z_h, Q2, pT, R0, R1, R2);
BinVariable::xFillVectors(double z_h, double Q2, double pT, double R0, double R1, double R2) {
    v_z_h.push_back(z_h);
    v_Q2.push_back(Q2);
    v_pT.push_back(pT);
    v_R0.push_back(R0);
    v_R1.push_back(R1);
    v_R2.push_back(R2);
}
//mhFillVectors(x, z_h, Q2, pT, R0, R1, R2);
BinVariable::mhFillVectors(double x, double z_h, double Q2, double pT, double R0, double R1, double R2) {
    v_x.push_back(x);
    v_z_h.push_back(z_h);
    v_Q2.push_back(Q2);
    v_pT.push_back(pT);
    v_R0.push_back(R0);
    v_R1.push_back(R1);
    v_R2.push_back(R2);
}

//Methods for calculating mean
BinVariable::meanZ_h() {
    xmean = meanfunc(v_x);
    Q2mean = meanfunc(v_Q2);
    pTmean = meanfunc(v_pT);
    R0mean = meanfunc(v_R0);
    R1mean = meanfunc(v_R1);
    R2mean = meanfunc(v_R2);
}
BinVariable::meanx() {
    z_hmean = meanfunc(v_z_h);
    Q2mean = meanfunc(v_Q2);
    pTmean = meanfunc(v_pT);
    R0mean = meanfunc(v_R0);
    R1mean = meanfunc(v_R1);
    R2mean = meanfunc(v_R2);
}
BinVariable::meanmh() {
    z_hmean = meanfunc(v_z_h);
    xmean = meanfunc(v_x);
    Q2mean = meanfunc(v_Q2);
    pTmean = meanfunc(v_pT);
    R0mean = meanfunc(v_R0);
    R1mean = meanfunc(v_R1);
    R2mean = meanfunc(v_R2);
}