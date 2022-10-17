#include "MCClasses.h"
#include "Kinematics.h"
#include "BinVariable.h"
#include "MCParticle.h"
#include "MultiParticle.h"
#include "Pidi.h"
#include "Quark.h"
#include "Diquark.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;
int LundAnalysisH(
                 const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3051_0.hipo",
//                 const char * rootfile = "../OutputFiles/AffinityFiles/Files_9_16/file1.root"
                 const char * rootfile = "../OutputFiles/Test10_13/file1.root"
//                 const char * rootfile = "OutputFiles/AffinityFiles/Files_9_16/TMD1.root"
//                 const char * rootfile = "OutputFiles/AffinityFiles/Files_9_12/collinear1.root"
)
{
    gROOT->ProcessLine("#include <vector>");
    //Below file is now disappeared...
//    auto hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3301_3.hipo";
    
// Current files: defined in main function though
//    auto hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3051_0.hipo";
//    auto rootFile = "OutputFiles/AffinityFiles/Files_9_5/Exactfile2.root";
    
    TFile *f = TFile::Open(rootfile,"RECREATE");
    
    HipoChain chain;
    
    //Add file to HipoChain
    chain.Add(hipoFile);
    auto config_c12 = chain.GetC12Reader();
    
    //Set PID cuts
    config_c12->addExactPid(11,1);    //exactly 1 electron
//    config_c12->addAtLeastPid(211,1);    //exactly 1 pi+
//    config_c12->addAtLeastPid(-211,1);    //exactly 1 pi-
    config_c12->addExactPid(2212,1);    //exactly 1 proton

    //Constants 
    double electron_beam_energy = 10.6; //(fall2018)
    double electronMass = 0.000511;
    double protonMass = 0.938272;
    double pipluspid = 211;
    double piminuspid =- 211;
    double uppid = 2;
    double downpid = 1;
    double antiuppid = -2;
    double antidownpid = -1;
    //
    //Initialization
    //

    // Debugging counters
    int zbin0count = 0;
    int zbin1count = 0;
    int zbin2count = 0;
    int zbin3count = 0;
    int zbin4count = 0;
    int zbin5count = 0;
    int zbin6count = 0;
    
    //Kinematics class for functions
    Kinematics _k;
    //MC::Lund bank entries
    int pid;
    int id;
    double px;
    double py;
    double pz;
    int daughter;
    int parent;
    double mass;
    double vz;

    //Calculated SIDIS kinematics
    double cth;
    double Q2;
    double x;
    double pt_lab;
    double z_h;
    double zpiplus;
    double zpiminus;
    
    //Cut Kinematics
    double W;
    double xFpiplus;
    double xFpiminus;
    double nu;
    double Mx;
    
    //Affinity ratios
    double R0;
    double R1;
    double R2;
    int qparent;
    int diparent;
    double s;
    double y;

    //Dihadron kinematics
    double Mdihadron;
    double dihadronpt;
    double Pdihadron;

    
    //MCParticle variables
    double P;
    double E;
    
    //Initializing particle vectors
    TLorentzVector q;
    TLorentzVector init_electron;
    TLorentzVector init_target;
    TLorentzVector dihadron;
    
    TLorentzVector k;
    TLorentzVector kf;
    TLorentzVector ki;
    TVector2 deltak; // Transverse light cone vector - (V_x,V_y)
    double deltakx;
    double deltaky;
    double deltak2;
    
    //gNframe
    TLorentzVector lv_p1_gN;
    TLorentzVector lv_p2_gN;
    TLorentzVector dihadron_gN;
    
    TLorentzVector lv_q_gN;
    
    TLorentzVector target_gN;
    
    TLorentzVector ki_gN;
    TLorentzVector kf_gN;
    TLorentzVector k_gN;
    
    double pt_gN;
    
    TLorentzVector gN;
    TVector3 gNBoost;
    TVector3 gNBoostNeg;
    //For checking momentum conservation / if the lundstring contains momentum for all hadrons
//    TLorentzVector diquark;
    TLorentzVector lundstring;
//    TLorentzVector photon;
//    TLorentzVector proton;
    
    //Breit frame variables
    TLorentzVector Breit;
    TLorentzVector Breit_target;
    TVector3 BreitBoost;
    TLorentzVector kfBreit;
    TVector2 kfBreitTran;
    TLorentzVector dihadronBreit;
    TVector2 dihadronBreitTran;
    
    //Photon Frame variables
    TLorentzVector PFFrame;
    TVector3 PFBoost;
    TLorentzVector qPF;
    TVector3 qPFVect;
    TVector3 qPFVectUnit;
    TVector3 zAxis(0,0,1);
    double PFAngle;
    TVector3 PFAxis;
    
    double qPFMinus;
    TLorentzVector dihadronPF;
    double dihadronPFMinus;
    
    double z_N;
    TVector2 q_T;
    
    // Bin objects for collecting kinematic variables
    
    BinVariable zbin0 = BinVariable();
    BinVariable zbin1 = BinVariable();
    BinVariable zbin2 = BinVariable();
    BinVariable zbin3 = BinVariable();
    BinVariable zbin4 = BinVariable();
    BinVariable zbin5 = BinVariable();
    BinVariable zbin6 = BinVariable();
    
    BinVariable xbin0 = BinVariable();
    BinVariable xbin1 = BinVariable();
    BinVariable xbin2 = BinVariable();
    BinVariable xbin3 = BinVariable();
    BinVariable xbin4 = BinVariable();
    BinVariable xbin5 = BinVariable();
    BinVariable xbin6 = BinVariable();
    
    BinVariable Mhbin0 = BinVariable();
    BinVariable Mhbin1 = BinVariable();
    BinVariable Mhbin2 = BinVariable();
    BinVariable Mhbin3 = BinVariable();
    BinVariable Mhbin4 = BinVariable();
    BinVariable Mhbin5 = BinVariable();
    BinVariable Mhbin6 = BinVariable();
    
    vector<double> xbins{0.1,0.13,0.16,0.19,0.235,0.3,0.5};
    vector<double> zbins{0.35,0.43,0.49,0.55,0.62,0.7,0.83};
    vector<double> Mhbins;
    for(int i = 0;i < 7; i++) {
        Mhbins.push_back(0.3 + i / 6.);
    }
    //Vectors for calculating means
    vector<BinVariable> zbinv = {zbin0, zbin1, zbin2, zbin3, zbin4, zbin5, zbin6};
    vector<BinVariable> xbinv = {zbin0, zbin1, zbin2, zbin3, zbin4, zbin5, zbin6};       
    vector<BinVariable> Mhbinv = {zbin0, zbin1, zbin2, zbin3, zbin4, zbin5, zbin6};
    vector<string> vinfoString = {"0th bin", "1st bin", "2nd bin", "3rd bin", "4th bin", "5th bin", "6th bin"};
    
    //Add MC::Lund bank for taking Lund data
    auto idx_MCLund= config_c12->addBank("MC::Lund");
    //Add a few items
    auto iPid=config_c12->getBankOrder(idx_MCLund,"pid");
    auto ipx=config_c12->getBankOrder(idx_MCLund,"px"); 
    auto ipy=config_c12->getBankOrder(idx_MCLund,"py");
    auto ipz=config_c12->getBankOrder(idx_MCLund,"pz");
    auto idaughter=config_c12->getBankOrder(idx_MCLund,"daughter");
    auto iparent=config_c12->getBankOrder(idx_MCLund,"parent");
    auto imass=config_c12->getBankOrder(idx_MCLund,"mass");
    auto ivz = config_c12->getBankOrder(idx_MCLund,"vz");
    
    //Creating vectors to fill using push_back in loop
    
    vector<int> vdiquarklist;
    vector<int> vhadronlist;
    vector<int> vquarklist;
    vquarklist = {-8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8};
    vdiquarklist = {1103, 2101, 2103, 2203, 3101, 3103, 3201, 3203, 3303, 4101, 4103, 4201, 4203, 4301, 4303, 4403, 5101, 5103, 5201, 5203, 5301, 5303, 5401, 5403, 5503};
    vhadronlist = {-3122, -211, 111, 211, 1114, 2114, 2212, 2214, 2224, 3112, 3114, 3122, 3214, 3222, 3224, 3312, 3324, -323, -313, -213, 113, 213, 221, 223, 310, 313, 323, 331, 333};
    
    int hash_count = 0;
    
    //Making new MC tree
    TTree *t = new TTree("tree_MC","Tree with MC data");

    t->Branch("z_h",&z_h);
    t->Branch("x",&x);
    t->Branch("pt",&pt_gN);
    t->Branch("Q2",&Q2);
    t->Branch("Ph",&Pdihadron);
    t->Branch("Mdihadron",&Mdihadron); //dihadron mass
    t->Branch("R0",&R0); //initial parton momentum
    t->Branch("R1",&R1); //final parton momentum
    t->Branch("R2",&R2);
//    t->Branch("PFkix",&PFkix); //Photon frame partonic momentum for checking ki values
//    t->Branch("PFkiy",&PFkiy);
//    t->Branch("PFkiz",&PFkiz);
//    t->Branch("PFkit",&PFkit);
//    t->Branch("PFkfx",&PFkfx);
//    t->Branch("PFkfy",&PFkfy);
//    t->Branch("PFkfz",&PFkfz);
//    t->Branch("PFkft",&PFkft);

    
    //Tell the user that the loop is starting
    cout << "Start Event Loop" << endl;
    int tree_count = 0;
    //now get reference to (unique)ptr for accessing data in loop
    //this will point to the correct place when file changes
    //
    //This line comes from AnalysisWithExtraBanks.C
    auto& c12=chain.C12ref();
    int event_count = 0;
    
    //Loop over all events in the file
    while(chain.Next()==true){
        
        event_count += 1;
        if(event_count == 1) {
            cout << '\n';
            cout << "\033[96m";
            cout << "\t\t\t\t\t\t\t\t" << " ~~~~~~~~~~~~" << '\n';
            cout << "\t\t\t\t\t\t\t\t" << "|Progress Bar|" << '\n';
            cout << "\t\t\t\t\t\t\t\t" << " ~~~~~~~~~~~~" << '\n';
            cout << "\t\t[";
        }
        if(event_count % 16388 == 0) {
            
            hash_count += 1;
//            cout << "\033[A" << "\033[A";
            cout << '\r' << "\t\t[";
            for (int i = 1; i < hash_count + 1;i++) {
                cout << '#';
            }
            for (int i = 0; i < 100 - hash_count; i++) {
                cout << ' ';
            }
            cout << ']' << ' ';
            cout << event_count / 16388. << '%';
            if(event_count / 16388. == 100) {
                cout << endl;
            }
            cout << flush;
        }
        //Break at event 100 for testing with shorter run time
/*        if(event_count >= 1000) {
            cout << "Breaking at event: " << event_count << '\n';
            break;
        }*/
        if(c12->getDetParticles().empty())
            continue;
        
        //Intializing MCParticles
        MCParticle electron = MCParticle();
        MCParticle proton = MCParticle();
        MCParticle photon = MCParticle();
        MCParticle Lund = MCParticle();

        Pidi piplus = Pidi();
        Pidi piminus = Pidi();

        Quark quark = Quark();

        Pidi diquark = Pidi();
        
        MultiParticle Hadron = MultiParticle();
        //Loop over MC::Lund entries in this event using its ID = idx_MCLund
        //Get PID from its id = iPid
        for(auto imc=0;imc<c12->getBank(idx_MCLund)->getRows();imc++){
            auto mcparticles = c12->mcparts();
            
            id = mcparticles->getIndex(imc);
            pid = mcparticles->getPid(imc);
            px = mcparticles->getPx(imc);
            py = mcparticles->getPy(imc);
            pz = mcparticles->getPz(imc);
            daughter = mcparticles->getDaughter(imc);
            parent = mcparticles->getParent(imc);
            mass = mcparticles->getMass(imc);
            P = _k.Pfunc(px,py,pz);
            E = _k.Efunc(mass,P);
            vz = mcparticles->getVz(imc);
            //
            //Kinematics
            // 


            //Setting scattered electron
            if(pid==11 && parent==1){
                electron.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                electron.setVectors();
            }
            //pi+
            else if(pid==pipluspid){
                piplus.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                piplus.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
            //pi-
            else if(pid==piminuspid){
                piminus.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                piminus.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
            //all quarks
            else if(std::count(vquarklist.begin(), vquarklist.end(), pid)){
                quark.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                quark.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
            //MCParticle
            else if(pid==92 || pid == 91){
                Lund.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                Lund.setVectors();
            }
            
            else if(std::count(vdiquarklist.begin(), vdiquarklist.end(), pid)){
                diquark.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                diquark.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
            else if(pid == 22 && parent == 1){
                photon.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                photon.setVectors();
            }
            else if(id == 2){
                proton.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                proton.setVectors();
            }
            else if(std::count(vhadronlist.begin(), vhadronlist.end(), pid)) {
                Hadron.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                Hadron.update(id, pid, px, py, pz, daughter, parent, 
                              mass, vz);
            }
        }
        
        //Selecting pions that come from Lund particle
        for(int i = 0; i < piplus.v_id.size(); i++) {
            if(piplus.v_parent[i] == Lund.id) {
                piplus.select_id = i;
            }
        }
        for(int i = 0; i < piminus.v_id.size(); i++) {
            if(piminus.v_parent[i] == Lund.id) {
                piminus.select_id = i;
            }
        }
        
        //Selecting initial quark
        for(int i = 0; i < quark.v_id.size(); i++) {
            if(quark.v_parent[i] == 0) {
                quark.final_id = i;
            }
        }
        
        //Selecting diquark
        for(int i = 0; i < diquark.v_id.size(); i++) {
            if(diquark.v_parent[i] == 2) {
                diquark.select_id = i;
            }
        }
        
        //Skip non-dipion events
        if(piplus.select_id == -999) {
            continue;
        }
        else if(piminus.select_id == -999) {
            continue;
        }
        if(diquark.select_id != -999) {
            diquark.fillParticle(diquark.v_id[diquark.select_id], diquark.v_pid[diquark.select_id], diquark.v_px[diquark.select_id], diquark.v_py[diquark.select_id], 
                           diquark.v_pz[diquark.select_id], diquark.v_daughter[diquark.select_id], diquark.v_parent[diquark.select_id], diquark.v_mass[diquark.select_id], diquark.v_vz[diquark.select_id]);
            diquark.setVectors();
        }
        piplus.fillParticle(piplus.v_id[piplus.select_id], piplus.v_pid[piplus.select_id], piplus.v_px[piplus.select_id], piplus.v_py[piplus.select_id], 
                           piplus.v_pz[piplus.select_id], piplus.v_daughter[piplus.select_id], piplus.v_parent[piplus.select_id], piplus.v_mass[piplus.select_id], piplus.v_vz[piplus.select_id]);
        piplus.setVectors();
        piminus.fillParticle(piminus.v_id[piminus.select_id], piminus.v_pid[piminus.select_id], piminus.v_px[piminus.select_id], piminus.v_py[piminus.select_id], 
                           piminus.v_pz[piminus.select_id], piminus.v_daughter[piminus.select_id], piminus.v_parent[piminus.select_id], piminus.v_mass[piminus.select_id], piminus.v_vz[piminus.select_id]);
        piminus.setVectors();
        
        quark.fillParticle(quark.v_id[quark.final_id], quark.v_pid[quark.final_id], quark.v_px[quark.final_id], quark.v_py[quark.final_id], 
                           quark.v_pz[quark.final_id], quark.v_daughter[quark.final_id], quark.v_parent[quark.final_id], quark.v_mass[quark.final_id], quark.v_vz[quark.final_id]);
        quark.setVectors();
       
        //Setting inital beam and target particles
        init_electron.SetPxPyPzE(0, 0, sqrt(electron_beam_energy * electron_beam_energy - electronMass * electronMass), electron_beam_energy);
        init_target.SetPxPyPzE(0, 0, 0, proton.E);
        
        
        dihadron = piplus.lv + piminus.lv;
        Mdihadron = dihadron.M();
        Pdihadron = dihadron.P();
        q = init_electron - electron.lv; //virtual photon
        
        //Missing mass
        Mx = _k.Mxfunc(q, init_target, piplus.lv, piminus.lv);

        cth = _k.cthfunc(electron.px,electron.py,electron.pz);
        Q2 = _k.Q2func(electron_beam_energy,electron.E,cth); //Momentum transfer
        zpiplus = (init_target * piplus.lv) / (init_target * q);
        zpiminus = (init_target * piminus.lv) / (init_target * q);
        z_h = zpiplus + zpiminus;
        s = _k.sfunc(protonMass, electronMass, electron_beam_energy);
        y = _k.yfunc(electron_beam_energy,electron.E);
        x = Q2/s/y; // Bjorken x
        pt_lab = _k.Ptfunc(dihadron.Px(), dihadron.Py()); //hadron transverse momentum
        
        kf = quark.lv;
        ki = kf - q;
        //Cut Kinematics
        
        gN = q;
        gN += init_target;
        gNBoost = gN.BoostVector();
        gNBoostNeg = -gNBoost;
        
        lv_p1_gN = piplus.lv;
        lv_p2_gN = piminus.lv;
        lv_p1_gN.Boost(gNBoostNeg);
        lv_p2_gN.Boost(gNBoostNeg);
        
        lv_q_gN = q;
        lv_q_gN.Boost(gNBoostNeg);
        
        //Need dihadron in gN frame for pT
        dihadron_gN = dihadron;
        dihadron_gN.Boost(gNBoostNeg);
        pt_gN = _k.Ptfunc(dihadron_gN);
        
        //Need target in gN
        target_gN = init_target;
        target_gN.Boost(gNBoostNeg);
        
        //Need partonic in gN
        ki_gN = ki;
        ki_gN.Boost(gNBoostNeg);
        
        kf_gN = kf;
        kf_gN.Boost(gNBoostNeg);
        
        //Feynman x
        xFpiplus = _k.xFfunc(lv_p1_gN,lv_q_gN,W);
        xFpiminus = _k.xFfunc(lv_p2_gN,lv_q_gN,W);
        
        //nu and W
        nu = _k.nufunc(electron_beam_energy,electron.E);
        W = _k.Wfunc(Q2,protonMass,nu);
        
        // Breit Frame Kinematics for delta k
        Breit = q;
        Breit_target.SetPxPyPzE(0,0,0,2 * x * protonMass); // E^2 = M^2 + P^2 --> P = 0 so E = M = 2 * x * protonmass
        Breit += Breit_target;
        BreitBoost = Breit.BoostVector();
        BreitBoost = -1 * BreitBoost;
        kfBreit = kf;
        kfBreit.Boost(BreitBoost);
        kfBreitTran = _k.PtVectfunc(kfBreit); //kfbT in delta k calculation - needs to be a transverse light cone vector of form (V_x, V_y)
        
        dihadronBreit = dihadron;
        dihadronBreit.Boost(BreitBoost);
        dihadronBreitTran = _k.PtVectfunc(dihadronBreit); //PBbT in qT part of delta k calculation
        
        PFFrame = q + init_target;
        PFBoost = PFFrame.BoostVector();
        PFBoost = -1 * PFBoost;
        qPF = q;
        qPF.Boost(PFBoost);
        qPFVect = qPF.Vect();
        qPFVectUnit = qPFVect.Unit();
        PFAngle = qPFVectUnit.Angle(zAxis);
        PFAxis = qPFVectUnit.Cross(zAxis);
        //To rotate -> vector.Rotate(PFAngle,PFAxis);
        
        //
        //Photon Frame
        //
        
        //Dihadron
        dihadronPF = dihadron;
        dihadronPF.Boost(PFBoost);
        dihadronPF.Rotate(PFAngle,PFAxis);
        dihadronPFMinus = _k.LightConeMinus(dihadronPF);
        //Virtual Photon
        qPF.Rotate(PFAngle,PFAxis);
        qPFMinus = _k.LightConeMinus(qPF);
        //z_N and q_T
        z_N = dihadronPFMinus / qPFMinus;
        q_T = -1 * dihadronBreitTran / z_N;
        //ki, k, and delta k
        deltak = kfBreitTran - (-1 * z_N * q_T); 
        
        k = kf - q;
        k_gN = k;
        k_gN.Boost(gNBoostNeg);
        //These ratios are calculated in lab frame
//        R0 = R0func(ki, kf, deltak, Q2);
//        R1 = R1func(dihadron, ki, kf);
//        R2 = R2func(k, Q2);
        
        //Ratios in gN frame
        R0 = _k.R0func(ki_gN, kf_gN, deltak, Q2);
        R1 = _k.R1func(dihadron_gN, ki_gN, kf_gN);
        R2 = _k.R2func(k_gN, Q2);
        
        //CUTS:
        //Region cuts:
        
/*
        //TMD and Collinear region selection:
        if(R0 > 0.3 || R1 > 0.3) {
            continue;
        }
        //TMD:
        if(R2 > 0.3) {
            continue;
        }
        //Collinear
//        if(R2 < 0.9) {
//            continue;
//        }
         
*/
        //Missing mass
        if(Mx <= 1.5) {
            continue;
        }
        
        //Feynman x
        if(xFpiplus <= 0 || xFpiminus <= 0) {
            continue;
        }
        
        //Vertex Position
        if(abs(electron.vz - piplus.vz) >= 20) {
            continue;
        }
        if(abs(electron.vz - piminus.vz >= 20)) {
            continue;
        }
        if(electron.vz <= -8 || electron.vz >= 3) {
            continue;
        }
        
        //Channel Selection
        
        //Virtual photon mass / momentum transfer
        if(Q2 <= 1 || Q2 >= 100) {
            continue;
        }
        //Hadronic system mass
        if(W <= 2 || W >= 100) {
            continue;
        }
        //Energy fraction
        if(y <= 0 || y >= 0.8) {
            continue;
        }
        if(piplus.P <= 1.25 || piminus.P <= 1.25) {
            continue;
        }
        
        tree_count += 1;
        t->Fill();
        /*
        //Need: x, z, Q2, pT, R0, R1, R2
        //zbins:
        for(int i = 0; i < zbins.size(); i++) {
            if(z_h <= zbins[i]) {
                zbinv[i].zFillVectors(x, Q2, pt_gN, R0, R1, R2);
                break;
            }
        }
        //Mh bins
        for(int i = 0; i < Mhbins.size(); i++) {
            if(Mdihadron <= Mhbins[i]) {
                Mhbinv[i].mhFillVectors(x, z_h, Q2, pt_gN, R0, R1, R2);
                break;
            }
        }
        //x bins
        for(int i = 0; i < xbins.size(); i++) {
            if(x <= xbins[i]) {
                xbinv[i].zFillVectors(z_h, Q2, pt_gN, R0, R1, R2);
                break;
            }
        }
        //print out tree count every 100 to give update to user
        if(tree_count % 100 == 0) {
//            cout << "Tree_count: " << tree_count << '\n';
        }
        */
    }//Event loop end bracket
    cout << "\033[0m" << "\033[49m";
    /*
    cout << "Final tree_count: " << tree_count;
    //Making new Affinity trees
    TTree *t_z_h = new TTree("tree_z_h_bins","Tree with mean values binned by z_h affinity calculations");
    TTree *t_x = new TTree("tree_x_bins","Tree with mean values binned by x affinity calculations");
    TTree *t_Mh = new TTree("tree_Mh_bins","Tree with mean values binned by Mh affinity calculations");
    
    string infoString;
    Double_t z_h_t;
    Double_t x_t;
    Double_t Q2_t;
    Double_t pT_t;
    Double_t R0_t;
    Double_t R1_t;
    Double_t R2_t;
    
    t_z_h->Branch("Name",&infoString);
    t_z_h->Branch("x", &x_t);
    t_z_h->Branch("Q2", &Q2_t);
    t_z_h->Branch("pT", &pT_t);
    t_z_h->Branch("R0", &R0_t);
    t_z_h->Branch("R1", &R1_t);
    t_z_h->Branch("R2", &R2_t);
    
    t_x->Branch("Name",&infoString);
    t_x->Branch("z_h", &z_h_t);
    t_x->Branch("Q2", &Q2_t);
    t_x->Branch("pT", &pT_t);
    t_x->Branch("R0", &R0_t);
    t_x->Branch("R1", &R1_t);
    t_x->Branch("R2", &R2_t);
    
    t_Mh->Branch("Name",&infoString);
    t_Mh->Branch("x", &x_t);
    t_Mh->Branch("z_h", &z_h_t);
    t_Mh->Branch("Q2", &Q2_t);
    t_Mh->Branch("pT", &pT_t);
    t_Mh->Branch("R0", &R0_t);
    t_Mh->Branch("R1", &R1_t);
    t_Mh->Branch("R2", &R2_t);
    
    //Calculating means
    //Setting zbin means
    for(int i = 0; i < vinfoString.size(); i++) {
        zbinv[i].meanZ_h();
        infoString = vinfoString[i];
        x_t = zbinv[i].xmean;
        Q2_t = zbinv[i].Q2mean;
        pT_t = zbinv[i].pTmean;
        R0_t = zbinv[i].R0mean;
        R1_t = zbinv[i].R1mean;
        R2_t = zbinv[i].R2mean;
        t_z_h->Fill();
        }
    for(int i = 0; i < vinfoString.size(); i++) {
        xbinv[i].meanx();
        infoString = vinfoString[i];
        z_h_t = xbinv[i].z_hmean;
        Q2_t = xbinv[i].Q2mean;
        pT_t = xbinv[i].pTmean;
        R0_t = xbinv[i].R0mean;
        R1_t = xbinv[i].R1mean;
        R2_t = xbinv[i].R2mean;
        t_x->Fill();
        }
    for(int i = 0; i < vinfoString.size(); i++) {
        Mhbinv[i].meanmh();
        infoString = vinfoString[i];
        x_t = Mhbinv[i].xmean;
        z_h_t = Mhbinv[i].z_hmean;
        Q2_t = Mhbinv[i].Q2mean;
        pT_t = Mhbinv[i].pTmean;
        R0_t = Mhbinv[i].R0mean;
        R1_t = Mhbinv[i].R1mean;
        R2_t = Mhbinv[i].R2mean;
        t_Mh->Fill();
        }
    */
    f->Write();
    delete f;
    
    return 0;
}