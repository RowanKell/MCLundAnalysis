#include <iostream>
#include <vector>
#include "TROOT.h"
#include "TTree.h"
using namespace std;
int LundAnalysis{
    
    gROOT->ProcessLine("#include <vector>");
    
    hipofile = "/path/to/hipo/";
    
    HipoChain chain;
    
    //Add file to HipoChain
    chain.Add(hipofile);
    auto config_c12 = chain.GetC12Reader();
    
    //Set PID cuts
    config_c12->addExactPid(11,1);    //exactly 1 electron
    config_c12->addExactPid(211,1);    //exactly 1 pi+
    config_c12->addExactPid(-211,1);    //exactly 1 pi-
    config_c12->addExactPid(2212,1);    //exactly 1 proton
    
    //Add MC::Lund bank for taking Lund data
    auto idx_MCLund= config_c12->addBank("MC::Lund");
    //Add a few items
    auto ipx=config_c12->getBankOrder(idx_MCLund,"px"); 
    auto ipy=config_c12->getBankOrder(idx_MCLund,"py");
    auto ipz=config_c12->getBankOrder(idx_MCLund,"pz");
    auto idaughter=config_c12->getBankOrder(idx_MCLund,"daughter");
    auto iparent=config_c12->getBankOrder(idx_MCLund,"parent");
    auto imass=config_c12->getBankOrder(idx_MCLund,"mass");
//    auto ivx=config_c12->getBankOrder(idx_MCLund,"vx");
//    auto ivy=config_c12->getBankOrder(idx_MCLund,"vy");
//    auto ivz=config_c12->getBankOrder(idx_MCLund,"vz");
    auto iE=config_c12->getBankOrder(idx_MCLund,"energy");
//    auto ipy=config_c12->getBankOrder(idx_MCLund,"type");
    
    //Creating vectors to fill using push_back in loop
    std::vector<float> px;
    std::vector<float> py;
    std::vector<float> pz;
    std::vector<float> daughter;
    std::vector<float> parent;
    std::vector<float> mass;
    std::vector<float> energy;
    //Making new MC tree
    TTree *t = new TTree("tree_MC","Tree with MC data");
    
    //Assigning variables to branches
    t->Branch("px",&px);
    t->Branch("py",&py);
    t->Branch("pz",&pz);
    t->Branch("daughter",&daughter);
    t->Branch("parent",&parent);
    t->Branch("mass",&mass);
    t->Branch("energy",&energy);
    
}