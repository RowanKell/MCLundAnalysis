#include "src/LundAnalysis.h"
#include <random>
#include <cmath>

//Main program of MCLundAnalysis - Reads Hipo file, calculates kinematics/Ratios, bins (sometimes) and saves
// to root file
int LundAnalysis_single_pion(

                    // hipoFile is the file we read in
                   const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3052_3.hipo",
                   // rootfile is the file we save data to
                   const char * rootfile = "/work/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/August_15/file_8.root"
)
{
    //I'm not sure why this is here, but I think the vector class isn't included by default?
    gROOT->ProcessLine("#include <vector>");
    
    //Create the rootfile - if it exists, recreate it
    TFile *f = TFile::Open(rootfile,"RECREATE");
    
    //CLAS12 object to help load hipo files
    HipoChain chain;
    
    //Add file to HipoChain - could have multiple, but just run 1 per job and run many jobs to improve efficiency
    chain.Add(hipoFile);
    auto config_c12 = chain.GetC12Reader();
    
    //Pre-select only events with 1 electron and 1 proton | SHOULD CHECK IF THESE SHOULD BE ATLEAST
    config_c12->addAtLeastPid(11,1);    //exactly 1 electron
    config_c12->addAtLeastPid(2212,1);    //exactly 1 proton
    
    //Add MC::Lund bank for taking Lund data
    auto idx_MCLund= config_c12->addBank("MC::Lund");
    
    // These are used to create the loading bar, just aesthetic
    int hash_count = 0;

    int num_events = 18297;
         
    //Making new MC tree for dihadron·
    //This tree is used for BOX affinity calculations, diff trees are used for AffinityCalc
    //We do not bin in this tree, rather we bin later for tree_MC
    TTree *tree_MC = new TTree("tree_MC","Tree with MC data from dihadron");
    
    //These are the kinematics/Ratios we need to bin and calculate affinity BOX style
    tree_MC->Branch("z",&z_h);
    tree_MC->Branch("x",&x);
    tree_MC->Branch("pT_BF",&pT_BF);
    tree_MC->Branch("Q2",&Q2);
    tree_MC->Branch("R0",&R0); //initial parton momentum
    tree_MC->Branch("R1",&R1); //Measured in gN frame
    tree_MC->Branch("R2",&R2);
    tree_MC->Branch("R2_adjust",&R2_adjust);
    tree_MC->Branch("qTQ_HF",&qTQ_HF);
    tree_MC->Branch("qT_HF",&qT_HF);
    tree_MC->Branch("qT_from_pT_HF_z",&qT_from_pT_HF_z);
    tree_MC->Branch("qT_from_pT_PF_z",&qT_from_pT_PF_z);
    tree_MC->Branch("qT_from_pT_BF",&qT_from_pT_BF);
    tree_MC->Branch("qT_from_zN",&qT_from_zN);
    tree_MC->Branch("zeta", &zeta);
    tree_MC->Branch("xi", &xi);
    
    TTree *tree_maxmin = new TTree("tree_maxmin","Tree with max and min Ri values");
    
    tree_maxmin->Branch("R0_max",&R0_max); //initial parton momentum
    tree_maxmin->Branch("R1_max",&R1_max); //Measured in gN frame
    tree_maxmin->Branch("R2_max",&R2_max);
    
    tree_maxmin->Branch("R0_min",&R0_min); //initial parton momentum
    tree_maxmin->Branch("R1_min",&R1_min); //Measured in gN frame
    tree_maxmin->Branch("R2_min",&R2_min);
    
//     TTree *tree_test = new TTree("tree_test","Tree with kinematics for viewing distributions");
//     tree_test->Branch("delta_k_T",&delta_k_T);
//     tree_test->Branch("M_ki",&M_ki);
//     tree_test->Branch("M_kf",&M_kf);
    
//     tree_test->Branch("r_delta_k_T",&r_delta_k_T);
//     tree_test->Branch("r_M_ki",&r_M_ki);
//     tree_test->Branch("r_M_kf",&r_M_kf);
    
    qTQ_low_aff = 1.102616;
    qTQ_low_tol = 0.11;
    qTQ_high_aff = 0.15;
    qTQ_high_tol = 0.015;
    
    Q2_low_aff = 3;
    Q2_low_tol = 0.6;
    pT_low_aff = 0.3;
    pT_low_tol = 0.3;
    x_low_aff = 0.15;
    x_low_tol = 0.15;
    z_low_aff = 0.3;
    z_low_tol = 0.3;
    
    Q2_high_aff = 3;
    Q2_high_tol = 0.6;
    pT_high_aff = 0.1;
    pT_high_tol = 0.1;
    x_high_aff = 0.15;
    x_high_tol = 0.15;
    z_high_aff = 0.3;
    z_high_tol = 0.3;
    
    
    TTree *tree_low = new TTree("tree_low","Tree with kinematics for lower TMD affinity bin only");
    tree_low->Branch("M",&M);
    tree_low->Branch("M_h",&m_1);
    tree_low->Branch("pT_BF",&pT_BF);
    tree_low->Branch("x",&x);
    tree_low->Branch("z",&z_h);
    tree_low->Branch("Q2",&Q2);
    tree_low->Branch("R0",&R0); //initial parton momentum
    tree_low->Branch("R1",&R1); //Measured in gN frame
    tree_low->Branch("R2",&R2);
    tree_low->Branch("R2_adjust",&R2_adjust);
    tree_low->Branch("qTQ_HF",&qTQ_HF);
    tree_low->Branch("qT_HF",&qT_HF);
    tree_low->Branch("M_ki",&M_ki);
    tree_low->Branch("M_kf",&M_kf);
    tree_low->Branch("delta_k_T",&delta_k_T);
    tree_low->Branch("ki_T",&ki_T);
    tree_low->Branch("xi",&xi);
    tree_low->Branch("zeta",&zeta);
    tree_low->Branch("theta_ki",&theta_ki);
    tree_low->Branch("theta_H",&theta_H);
    tree_low->Branch("theta_deltak",&theta_deltak);

    TTree *tree_high = new TTree("tree_high","Tree with kinematics for higher TMD affinity bin only");
    tree_high->Branch("M",&M);
    tree_high->Branch("M_h",&m_1);
    tree_high->Branch("pT_BF",&pT_BF);
    tree_high->Branch("x",&x);
    tree_high->Branch("z",&z_h);
    tree_high->Branch("Q2",&Q2);
    tree_high->Branch("R0",&R0); //initial parton momentum
    tree_high->Branch("R1",&R1); //Measured in gN frame
    tree_high->Branch("R2",&R2);
    tree_high->Branch("R2_adjust",&R2_adjust);
    tree_high->Branch("qTQ_HF",&qTQ_HF);
    tree_high->Branch("qT_HF",&qT_HF);
    tree_high->Branch("M_ki",&M_ki);
    tree_high->Branch("M_kf",&M_kf);
    tree_high->Branch("delta_k_T",&delta_k_T);
    tree_high->Branch("ki_T",&ki_T);
    tree_high->Branch("xi",&xi);
    tree_high->Branch("zeta",&zeta);
    tree_high->Branch("theta_ki",&theta_ki);
    tree_high->Branch("theta_H",&theta_H);
    tree_high->Branch("theta_deltak",&theta_deltak);
    
    int low_R1_count = 0;
    
    int up_count = 0;
    int down_count = 0;
    double quark_ratio = 0;
    double qTQmax = 0;
    double qTQcounts[9] = {};
    double xcounts[7] = {};
    double zcounts[7] = {};
    
    int one_bin = 0;
    double qTQcut_min = 0.3;
    double qTQcut_max = 0.5;
    
    double xcut_min = 0.235;
    double xcut_max = 0.3;
    
    double zcut_min = 0.43;
    double zcut_max = 0.49;
    
    //My bins
//     double xb[4] = {0.3,0.35,0.25,0.3};
//     double zb[4] = {0.42,0.48,0.48,0.54};
//     double pTb[4] = {0.26,0.39,0.26,0.39};
//     double Q2b[4] = {1.8,2.6,1.8,2.6};
    
    //Tetiana Bins
    double xb[4] = {0.1,0.15,0.5,0.6};
    double zb[4] = {0.3,0.36,0.4,0.8};
    double pTb[4] = {0.13,0.26,0,0};
    double Q2b[4] = {1.0,1.8,9,12};
    
    //prints this text, but just aesthetic
    cout << "Start Event Loop" << endl;
    int tree_count = 0;
    int passed_count = 0;
    
    //now get reference to (unique)ptr for accessing data in loop
    //this will point to the correct place when file changes
    //
    //This line comes from AnalysisWithExtraBanks.C
    auto& c12=chain.C12ref();
    //This value counts the number of events in total that pass the clas12reader cuts
    int event_count = 0;
    
    //Loop over all events in the file that pass proton+electron cuts
    while(chain.Next()==true){
//         if(tree_count > 10) {break;} //Uncomment this line to stop the program after 1000 events, useful for debugging/testing
        event_count += 1;
        //Aesthetics/loading bar
        if(event_count == 1) {
            cout << '\n';
            cout << "\033[96m";
            cout << "\t\t\t\t" << " ~~~~~~~~~~~~" << '\n';
            cout << "\t\t\t\t" << "|Progress Bar|" << '\n';
            cout << "\t\t\t\t" << " ~~~~~~~~~~~~" << '\n';
            cout << "\t\t[";
        }
        //This logic increments the loading bar
        if(event_count % num_events == 0) {
            
            hash_count += 1;
//            cout << "\033[A" << "\033[A";
            cout << '\r' << "\t[";
            for (int i = 1; i < hash_count + 1;i++) {
                cout << '#';
            }
            for (int i = 0; i < 100 - hash_count; i++) {
                cout << ' ';
            }
            cout << ']' << ' ';
            cout << event_count / num_events << '%';
            if(event_count / num_events == 100) {
                cout << endl;
            }
            cout << flush;
        }
        //Skip events with no particles
        if(c12->getDetParticles().empty())
            continue;
        
        //Intializing MCParticles - Apr 2024, not sure if these can be moved to declarations.C, but don't care tbh
        MCParticle electron;
        MCParticle proton;
        MCParticle photon;
        MCParticle Lund;

        Pidi pi_v;
        
        Pidi pi1;
        Pidi pi2;

        Quark quark;

        Pidi diquark;
        
        MultiParticle Hadron;
        //Loop over MC::Lund entries (particles) in this event using its ID = idx_MCLund
        for(auto imc=0;imc<c12->getBank(idx_MCLund)->getRows();imc++){
            auto mcparticles = c12->mcparts();
            
            // grab particle variables
            id = mcparticles->getIndex(imc);
            pid = mcparticles->getPid(imc);
            px = mcparticles->getPx(imc);
            py = mcparticles->getPy(imc);
            pz = mcparticles->getPz(imc);
            daughter = mcparticles->getDaughter(imc);
            parent = mcparticles->getParent(imc);
            mass = mcparticles->getMass(imc);
            P = Pfunc(px,py,pz);
            E = Efunc(mass,P);
            vz = mcparticles->getVz(imc);
            
            //
            //Kinematics
            // 

            if(pid == 2){up_count++;}
            if(pid == 1){down_count++;}

            //Setting scattered electron
            if(pid==11 && parent==1){
                //See src/MCParticle.C for these definitions
                electron.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                electron.setVectors();
            }
            //save pion variables
            else if(pid==pipluspid || pid ==piminuspid){
                //We use a vector to make sure we catch all pions so we can look at every dihadron pair
                pi_v.fillParticle(id, pid, px, py, pz, daughter, parent, mass, vz);
                pi_v.update(id, pid, px, py, pz, daughter, parent, 
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
        
        //Selecting initial quark
        for(int i = 0; i < quark.v_id.size(); i++) {
            if(quark.v_parent[i] == 0) {
                quark.final_id = i;
            }
        }
        
        //Selecting diquark - not really relevant tho
        for(int i = 0; i < diquark.v_id.size(); i++) {
            if(diquark.v_parent[i] == 2) {
                diquark.select_id = i;
            }
        }
        
        //Skip events without pions
        if(pi_v.v_id.size() < 1) {
            continue;
        }
        //Fill diquark if it exists
        if(diquark.select_id != -999) {
            diquark.fillParticle(diquark.v_id[diquark.select_id], diquark.v_pid[diquark.select_id], diquark.v_px[diquark.select_id], diquark.v_py[diquark.select_id], 
                           diquark.v_pz[diquark.select_id], diquark.v_daughter[diquark.select_id], diquark.v_parent[diquark.select_id], diquark.v_mass[diquark.select_id], diquark.v_vz[diquark.select_id]);
            diquark.setVectors();
        }
        

        //Loop over all combinations of pion pairs
        for(int i = 0; i < pi_v.v_id.size(); i++) {
            //
            // Set vectors from MC bank and beam/target knowledge
            //
            pi1.fillParticle(pi_v.v_id[i], pi_v.v_pid[i], pi_v.v_px[i], pi_v.v_py[i],pi_v.v_pz[i], pi_v.v_daughter[i], pi_v.v_parent[i], pi_v.v_mass[i], pi_v.v_vz[i]); 
            pi1.setVectors();


            quark.fillParticle(quark.v_id[quark.final_id], quark.v_pid[quark.final_id], quark.v_px[quark.final_id], quark.v_py[quark.final_id], quark.v_pz[quark.final_id], quark.v_daughter[quark.final_id], quark.v_parent[quark.final_id], quark.v_mass[quark.final_id], quark.v_vz[quark.final_id]);
            quark.setVectors();
            
            init_electron.SetPxPyPzE(0, 0, sqrt(electron_beam_energy * electron_beam_energy - electronMass * electronMass), electron_beam_energy);
            init_target.SetPxPyPzE(0, 0, 0, proton.E);

            //
            // SIDIS Kinematics in lab frame (invariants)
            //
            
            m_1 = pi1.mass;
            
            q = photon.lv;
            q_calc = init_electron - electron.lv;

            Mx = Mxfunc(q, init_target, pi1.lv);

            cth = cthfunc(electron.px,electron.py,electron.pz);
            Q2 = -(q * q);
            Q2_calc = Q2func(electron_beam_energy,electron.E,cth); //Momentum transfer

            z_h = (init_target * pi1.lv) / (init_target * q);
            s = sfunc(protonMass, electronMass, electron_beam_energy);
            y = yfunc(electron_beam_energy,electron.E);
            x = Q2/s/y;
            
            nu = nufunc(electron_beam_energy,electron.E);
            W = Wfunc(Q2,protonMass,nu);
            
            
            kf = quark.lv;
            ki = kf - q;

            
            
            //
            // Breit Frame Kinematics
            //
            
            //Boost stuff
            BF = q;
            BF_target.SetPxPyPzE(0,0,0,2 * x *protonMass); // E^2 = M^2 + P^2 --> P = 0 so E = M = 2 * x * protonmass
            BF += BF_target;
            BFBoost = BF.BoostVector();
            BFBoost = -1 * BFBoost;
            
            //Rotation
            q_BF_no_rot = q;
            q_BF_no_rot.Boost(BFBoost);
            q_BF_no_rot_VectUnit = q_BF_no_rot.Vect().Unit();
            BFAngle = q_BF_no_rot_VectUnit.Angle(zAxis);
            BFAxis = q_BF_no_rot_VectUnit.Cross(zAxis);
            
            //virtual photon
            q_BF = q;
            q_BF.Boost(BFBoost);
            q_BF.Rotate(BFAngle,BFAxis);
            
            //proton target
            proton_BF = init_target;
            proton_BF.Boost(BFBoost);
            proton_BF.Rotate(BFAngle,BFAxis);
            proton_BF_no_rot = init_target;
            proton_BF_no_rot.Boost(BFBoost);
            
            //partons
            kf_BF = kf;
            kf_BF.Boost(BFBoost);
            kfT_BF_vect = PtVectfunc(kf_BF); //BF frame boostkfbT in delta k calculation 
            kfT_BF = Ptfunc(kf_BF);


            M = proton_BF.M();
            
            //hadron (P_h)
            pi1_BF = pi1.lv;
            pi1_BF.Boost(BFBoost);
            pT_BF_vect = PtVectfunc(pi1_BF);
            pT_BF = Ptfunc(pi1_BF);
            
            //
            //Photon Frame
            //
            
            //Boosting
            PFFrame = q + init_target;
            PFBoost = PFFrame.BoostVector();
            PFBoost = -1 * PFBoost;
            
            //Rotation
            q_PF_no_rot = q;
            q_PF_no_rot.Boost(PFBoost);
            q_PF_no_rot_VectUnit = q_PF_no_rot.Vect().Unit();
            PFAngle = q_PF_no_rot_VectUnit.Angle(zAxis);
            PFAxis = q_PF_no_rot_VectUnit.Cross(zAxis);
            
            //proton, q, pion
            P_PF = init_target;
            q_PF = q;
            pi1_PF = pi1.lv;
            P_PF.Boost(PFBoost);
            q_PF.Boost(PFBoost);
            pi1_PF.Boost(PFBoost);
            P_PF.Rotate(PFAngle,PFAxis);
            q_PF.Rotate(PFAngle,PFAxis);
            pi1_PF.Rotate(PFAngle,PFAxis);
            
            //pT
            pT_PF = Ptfunc(pi1_PF);
            
            //partons
            ki_PF = ki;
            kf_PF = kf;
            ki_PF.Boost(PFBoost);
            kf_PF.Boost(PFBoost);
            ki_PF.Rotate(PFAngle,PFAxis);
            kf_PF.Rotate(PFAngle,PFAxis);
            
            //misc:
            xFpi1 = xFfunc(pi1_PF,q_PF,W);
            
            
            //
            //Hadron frame 
            //
            
            //Boosting
            Hadron_frame = pi1.lv;
            Hadron_frame += BF_target;
            HadronBoost = Hadron_frame.BoostVector();
            HadronBoost = -1 * HadronBoost;
            
            //rotation
            pi1_HF_no_rot = pi1.lv;
            pi1_HF_no_rot.Boost(HadronBoost);
            pi1_HF_no_rot_VectUnit = pi1_HF_no_rot.Vect().Unit();
            HFAngle = pi1_HF_no_rot_VectUnit.Angle(zAxis);
            HFAxis = pi1_HF_no_rot_VectUnit.Cross(zAxis);
            
            //outgoing hadron
            pi1_HF = pi1.lv;
            pi1_HF.Boost(HadronBoost);
            pi1_HF.Rotate(HFAngle,HFAxis);
            
            pT_HF = Ptfunc(pi1_HF);
            
            //virtual photon
            q_HF = q;
            q_HF.Boost(HadronBoost);
            q_HF.Rotate(HFAngle,HFAxis);
                
            qT_HF = Ptfunc(q_HF);
            qTQ_HF = qT_HF / pow(Q2,0.5);
            
            qT_from_pT_PF_z = pT_PF / z_h;
            


            //
            // z_N x_N calcs
            //
            
            //z_N = P_B,b minus / q_b minus
            P_B_b_minus = pi1_BF.Minus();
            q_b_minus = q_BF.Minus();
            z_N = P_B_b_minus / q_b_minus;
            x_N_kin =  (2 * x) / (1 + pow(1 + (4 * pow(x,2) * pow(M,2) / Q2),0.5));
            x_N = -1 * q_PF.Plus() / P_PF.Plus();

            
            qT_from_zN_vect = -1 * pT_BF_vect / z_N;

            qT_from_pT_BF_vect = pT_BF_vect / z_h;
            qT_from_pT_BF = Ptfunc(qT_from_pT_BF_vect);
            qTQ_from_pT_BF = qT_from_pT_BF / sqrt(Q2);


            //qT / Q for plotting
            qT_from_zN = Ptfunc(qT_from_zN_vect);

            //qT from pT/z
            qT_from_pT_HF_z = pT_HF / z_h;
            
            //ki, k, and delta k
            deltak = kfT_BF_vect - (-1 * z_N * qT_from_zN_vect); 

            k = kf - q;
            k_PF = k;
            k_PF.Boost(PFBoost);

            //Ratios in PF frame
            R0 = R0func(ki_PF, kf_PF, deltak, Q2);
            
            ki2 = abs(ki_PF * ki_PF);
            kf2 = abs(kf_PF * kf_PF);
            deltak2 = abs(deltak * deltak);
            
            if(deltak2 > ki2 && deltak2 > kf2) {
                R0check = 0;//DeltaK is biggest
            }
            else if(ki2 > kf2) {
                R0check = 1;//ki is biggest
            }
            else {
                R0check = 2;//kf is biggest
            }
            
            //Checking parton distributions
            M_ki = pow(ki2,0.5);
            M_kf = pow(kf2,0.5);
            delta_k_T = pow(deltak2,0.5);
            
            ki_T = Ptfunc(ki_PF);
            
            R1 = R1func(pi1_PF, ki_PF, kf_PF);
            if(R1 < -100) {
                low_R1_count++;
            }
            ki_BF = ki;
            ki_BF.Boost(BFBoost);
            kf_BF = kf;
            kf_BF.Boost(BFBoost);
            k_BF = k;
            k_BF.Boost(BFBoost);

            R2 = R2func(k_PF, Q2);
            R2_adjust = R2func_adjust(qT_HF, Q2);
            xF = xFpi1;
            
            //np params
            zeta = (z_N * q_BF.Minus()) / kf_BF.Minus();
            xi = (-1 * x_N * ki_BF.Plus()) / q_BF.Plus();
            
            //phi values
            theta_ki = ki_BF.Phi();
            theta_H = pi1_BF.Phi();
            theta_deltak = deltak.Phi();

            //Missing mass
            if(Mx <= 1.5) {
                continue;
            }

            //Feynman x
            if(xFpi1 <= 0) {
                continue;
            }

            //Vertex Position
            if(abs(electron.vz - pi1.vz) >= 20) {
                continue;
            }
            if(abs(electron.vz - pi2.vz) >= 20) {
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
            if(pi1.P <= 1.25) {
                continue;
            }
            if(one_bin) {
                if(x < xcut_min || x > xcut_max){continue;}
                if(z_h < zcut_min || z_h > zcut_max){continue;}
                if(qTQ_HF < qTQcut_min || qTQ_HF > qTQcut_max){continue;}
            }
            if(i == 0){passed_count++;}
            
            //Calculate r_M_ki etc from normal distributions
            r_M_ki = std::abs(distribution_M_ki(generator_M_ki));
            
            r_M_kf = std::abs(distribution_M_kf(generator_M_kf));
            
            r_delta_k_T = std::abs(distribution_delta_k_T(generator_delta_k_T));
            
            //Calculating max and min ratios
            if(R0 > R0_max) {R0_max = R0;}
            if(R1 > R1_max) {R1_max = R1;}
            if(R2 > R2_max) {R2_max = R2;}
            
            if(R0 < R0_min) {R0_min = R0;}
            if(R1 < R1_min) {R1_min = R1;}
            if(R2 < R2_min) {R2_min = R2;}
            
            
            
            tree_count += 1;
            tree_MC->Fill();
            low_bool = ((x > xb[0]) && (x < xb[1]) && (z_h > zb[0]) && (z_h < zb[1]) && (Q2 > Q2b[0]) && (Q2 < Q2b[1])&& (pT_PF > pTb[0]) && (pT_PF < pTb[1]));
            high_bool = ((x > xb[2]) && (x < xb[3]) && (z_h > zb[2]) && (z_h < zb[3]) && (qTQ_HF < 0.2) && (Q2 > Q2b[2]) && (Q2 < Q2b[3]));
            if(low_bool) {
                tree_low->Fill();
            }
            if(high_bool) {
                tree_high->Fill();
            }
            
            /*
            //zbins:
            for(int i = 0; i < zbins.size(); i++) {
                if(z_h <= zbins[i]) {
                    zbinv[i].zFillVectors(x, Q2, pT_gN, R0, R1, R2_adjust);
//                     cout << "filling z_h bin #" << i << " with: x - " << x << " | pT_gN: " << pT_gN << "\n";
                    zcounts[i] += 1;
                    break;
                }
            }
            //x bins
            for(int i = 0; i < xbins.size(); i++) {
                if(x <= xbins[i]) {
                    xbinv[i].xFillVectors(z_h, Q2, pT_gN, R0, R1, R2_adjust);
                    xcounts[i] += 1;
                    break;
                }
            }
            //Q2 bins
            for(int i = 0; i < Q2bins.size(); i++) {
                if(Q2 <= Q2bins[i]) {
                    Q2binv[i].Q2FillVectors(x, z_h, pT_gN, R0, R1, R2_adjust);
                    break;
                }
            }
            if(qTQ_pion > qTQmax){qTQmax = qTQ_pion;}
            //qTQ bins
            for(int i = 0; i < qTQbins.size(); i++) {
                if(qTQ_pion <= qTQbins[i]) {
                    qTQbinv[i].qTQFillVectors(x, z_h, Q2, pT_gN, R0, R1, R2_adjust);
                    qTQcounts[i] += 1;
                    break;
                }
            }*/
        }
	

    }
    cout << "\033[0m" << "\033[49m";
    cout << "Final event_count:" << event_count << '\n';
    cout << "Final tree_count: " << tree_count << '\n';
    cout << "Final events passed: " << passed_count << '\n';
    cout << "low_R1_count " << low_R1_count << '\n';

/*
    //Making new Affinity trees
    TTree *t_z_h = new TTree("tree_z_h_bins","Tree with mean values binned by z_h affinity calculations");
    TTree *t_x = new TTree("tree_x_bins","Tree with mean values binned by x affinity calculations");
    TTree *t_Q2 = new TTree("tree_Q2_bins","Tree with mean values binned by Q2 affinity calculations");
    TTree *t_qTQ = new TTree("tree_qTQ_bins","Tree with mean values binned by Q2 affinity calculations");
    
    
    string infoString;
    Double_t z_h_t;
    Double_t x_t;
    Double_t Q2_t;
    Double_t pT_t;
    
    Double_t z_h_t_1;
    Double_t pT_t_1;
    
    Double_t z_h_t_2;
    Double_t pT_t_2;
    
    Double_t R2_t; //R2 calculated as qT^2 / Q^2
    //z branch
    t_z_h->Branch("Name",&infoString);
    t_z_h->Branch("x", &x_t);
    t_z_h->Branch("Q2", &Q2_t);
    t_z_h->Branch("pT", &pT_t);
    t_z_h->Branch("R2", &R2_t);
    
    //x branch
    t_x->Branch("Name",&infoString);
    t_x->Branch("Q2", &Q2_t);
    
    t_x->Branch("z_h", &z_h_t);
    t_x->Branch("pT", &pT_t);
    t_x->Branch("R2", &R2_t);
    
    //Q2 branch
    t_Q2->Branch("Name",&infoString);
    t_Q2->Branch("x", &x_t);

    t_Q2->Branch("z_h", &z_h_t);
    t_Q2->Branch("pT", &pT_t);
    t_Q2->Branch("R2", &R2_t);
    
    //qTQ branch
    t_qTQ->Branch("Name",&infoString);
    t_qTQ->Branch("x", &x_t);
    t_qTQ->Branch("Q2", &Q2_t);
    
    t_qTQ->Branch("z_h", &z_h_t);
    t_qTQ->Branch("pT", &pT_t);
    t_qTQ->Branch("R2", &R2_t);
    //Calculating means
    //Setting zbin means
    for(int i = 0; i < vinfoString.size() - 2; i++) { //Note: we use i < vinfoString.size() - 2 for x,z,Mh 
                                                      //bc they have 7 bins, and there are 9 bins for qTQ, so vinfoString has length 9
        zbinv[i].meanZ_h(1);
        infoString = vinfoString[i];
        Q2_t = zbinv[i].Q2mean;
        x_t = zbinv[i].xmean;
        pT_t = zbinv[i].pTmean;
        R2_t = zbinv[i].R2mean;
//         cout << "z binning #" << i << ": xmean - " << zbinv[i].xmean << ": pTmean - " << zbinv[i].pTmean << "\n";
        t_z_h->Fill();
        }
    for(int i = 0; i < vinfoString.size() - 2; i++) {
        xbinv[i].meanx(1);
        infoString = vinfoString[i];
        Q2_t = xbinv[i].Q2mean;
        z_h_t = xbinv[i].z_hmean;
        pT_t = xbinv[i].pTmean;
        R2_t = xbinv[i].R2mean;
        t_x->Fill();
        }
    for(int i = 0; i < vinfoString.size() - 1; i++) {
        Q2binv[i].meanQ2(1);
        infoString = vinfoString[i];
        x_t = Q2binv[i].xmean;
        Q2_t = Q2binv[i].Q2mean;
        z_h_t = Q2binv[i].z_hmean;
        pT_t = Q2binv[i].pTmean;
        R2_t = Q2binv[i].R2mean;
        t_Q2->Fill();
        }
    
    for(int i = 0; i < vinfoString.size(); i++) {
        qTQbinv[i].meanqTQ(1);
        infoString = vinfoString[i];
        x_t = qTQbinv[i].xmean;
        Q2_t = qTQbinv[i].Q2mean;
        
        z_h_t = qTQbinv[i].z_hmean;
        pT_t = qTQbinv[i].pTmean;
        R2_t = qTQbinv[i].R2mean;
        t_qTQ->Fill();
        }
    */
    tree_maxmin->Fill();
    f->Write();
    delete f;
    return 0;
}