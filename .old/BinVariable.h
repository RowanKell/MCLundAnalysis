#ifndef BINVARIABLE_H
#define BINVARIABLE_H

class BinVariable {
    
    public:
    
    BinVariable();
    void zFillVectors(double x, double Q2, double pT, double R0, double R1, double R2);
    void xFillVectors(double z_h, double Q2, double pT, double R0, double R1, double R2);
    void mhFillVectors(double x, double z_h, double Q2, double pT, double R0, double R1, double R2);
    void meanZ_h();
    void meanx();
    void meanmh();
};
#endif