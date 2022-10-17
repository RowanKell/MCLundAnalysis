#include "MCParticle.h"
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