#include "Diquark.h"
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