#ifndef MCParticle_h
#define MCParticle_h

class MCParticle
{
    public:
    
    //Lund bank variables
    int pid;
    int id;
    double px;
    double py;
    double pz;
    int daughter;
    int parent;
    double mass;
    double P;
    double E;
    double vz;
    //TLorentzVector
    TLorentzVector lv;
    
    //Calculations
    double Pt;
    TVector2 PtVect;
    
    MCParticle();
    
    void inputPxPyPzM(double _px, double _py, double _pz, double _m);
    
    void SetParentDaughter(double _parent,double _daughter);
    
    void fillParticle(int _id, int _pid, double _px, double _py, double _pz, int _daughter, int _parent, double _mass, double _vz);
    
    void Calculate();
    
    void setVectors();
};

#endif