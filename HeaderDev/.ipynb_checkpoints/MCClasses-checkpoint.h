#ifndef MCClasses_h
#define MCClasses_h


class MultiParticle : public MCParticle
{
    public:
    
    MultiParticle();
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

class Pidi : public MultiParticle
{
    public:
    Pidi();
    int select_id;
    
};
class Quark : public MultiParticle
{
    public:
    
    Quark();
    int initial_id;
    int final_id;
};
class Diquark : public MultiParticle
{
    public:
    Diquark();
    int select_id;
    
    void diquarkReset();
};
class BinVariable
{
    
    public:
    //Need: x, z_h, Q2, pT, R0, R1, R2
    vector<double> v_x;
    vector<double> v_z_h;
    vector<double> v_Q2;
    vector<double> v_pT;
    vector<double> v_R0;
    vector<double> v_R1;
    vector<double> v_R2;
    
    double xmean;
    double z_hmean;
    double Q2mean;
    double pTmean;
    double R0mean;
    double R1mean;
    double R2mean;
    //zFillVectors(z_h, Q2, pT, R0, R1, R2)
    void zFillVectors(double x, double Q2, double pT, double R0, double R1, double R2);
    //xFillVectors(z_h, Q2, pT, R0, R1, R2);
    void xFillVectors(double z_h, double Q2, double pT, double R0, double R1, double R2);
    //mhFillVectors(x, z_h, Q2, pT, R0, R1, R2);
    void mhFillVectors(double x, double z_h, double Q2, double pT, double R0, double R1, double R2);
    
    //Methods for calculating mean
    void meanZ_h();
    void meanx();
    void meanmh();
};
#endif