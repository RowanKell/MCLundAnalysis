#include <vector>
#include <cmath>
#include "TROOT.h"
#pragma once
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

