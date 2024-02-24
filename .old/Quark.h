#ifndef Quark_h
#define Quark_h
#include "MultiParticle.h"
class Quark : public MultiParticle
{
    public:
    
    Quark();
    int initial_id;
    int final_id;
};
#endif