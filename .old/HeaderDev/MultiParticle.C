#include "MultiParticle.h"
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