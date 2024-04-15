#include <vector>
#include <cmath>
#include "TROOT.h"
#include "MultiParticle.C"
#pragma once
class Diquark : public MultiParticle
{
    public:
    
    int select_id  = -999;
    
    void diquarkReset()
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
};