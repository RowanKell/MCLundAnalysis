#include <vector>
#include <cmath>
#include "TROOT.h"
#include "MultiParticle.C"
#pragma once
class Quark : public MultiParticle
{
    public:
    
    int initial_id = -999;
    int final_id = -999;
};