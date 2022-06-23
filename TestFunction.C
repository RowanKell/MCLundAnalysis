#include <iostream>
#include <cmath>
using namespace std;
double P(double Px, double Py, double Pz)
{
    return sqrt(Px * Px + Py * Py + Pz * Pz);
}

double Q2(double E1, double E2, double cth)
{
    return 2 * E1 * E2 * (1.0 - cth);
}

int TestFunction()
{
    double myPx = 20;
    double myPy = 15;
    double myPz = 10;
    double myP = P(myPx, myPy, myPz);
    double myE1 = 2;
    double myE2 = 3;
    double mycth = 0.5;
    double myQ2 = Q2(myE1, myE2, mycth);
    std::cout << "This is myP: " << myP << endl;
    std::cout << "This is myQ2: " << myQ2 << endl;
    return 0;
}