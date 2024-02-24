#ifndef MULTIPARTICLE_H
#define MULTIPARTICLE_H
#include <iostream>
#include <vector>
#include <cmath>
#include "TROOT.h"
#include "MCParticle.h"
class MultiParticle : public MCParticle
{
    public:
    
    vector<int> v_id;
    vector<int> v_pid;
    vector<double> v_px;
    vector<double> v_py;
    vector<double> v_pz;
    vector<int> v_daughter;
    vector<int> v_parent;
    vector<double> v_mass;
    vector<double> v_vz;
    
    
    void update(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz);
    
};
#endif