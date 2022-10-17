#include "MCClasses.h"

MCParticle::MCParticle()
{
    //Lund bank variables
    pid = 0;
    id = 0;
    px = 0;
    py = 0;
    pz = 0;
    daughter = 0;
    parent = 0;
    mass = 0;
    P = 0;
    E = 0;
    vz = 0;
    
    //Calculations
    Pt = 0;
}

MCParticle::inputPxPyPzM(double _px, double _py, double _pz, double _m)
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
MCParticle::SetParentDaughter(double _parent,double _daughter)
{
    parent = _parent;
    daughter = _daughter;
}
MCParticle::Calculate()
{
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
    lv.SetPxPyPzE(px,py,pz,E);
    
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
}
MCParticle::setVectors()
{
    lv.SetPxPyPzE(px,py,pz,E);
    
    P = Pfunc(px, py, pz);
    E = Efunc(mass, P);
    
    Pt = Ptfunc(px, py);
    PtVect = PtVectfunc(lv);
}
MCParticle::fillParticle(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz)
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

MultiParticle::MultiParticle()
{
    //Lund bank variables
    pid = 0;
    id = 0;
    px = 0;
    py = 0;
    pz = 0;
    daughter = 0;
    parent = 0;
    mass = 0;
    P = 0;
    E = 0;
    vz = 0;
    
    //Calculations
    Pt = 0;
}
MultiParticle::update(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz)
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

Pidi::Pidi()
{
    select_id = -999; 
}
Quark::Quark()
{
    initial_id = -999;
    final_id = -999;
}
Diquark::Diquark()
{
    select_id  = -999;
}
Diquark::diquarkReset()
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

//zFillVectors(z_h, Q2, pT, R0, R1, R2)
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