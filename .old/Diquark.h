#ifndef Diquark_h
#define Diquark_h
#include <iostream>
#include <vector>
#include <cmath>
#include "TROOT.h"
#include "MultiParticle.h"
class Diquark : public MultiParticle
{
    public:
    Diquark();
    int select_id;
    
    void diquarkReset();
};
#endif