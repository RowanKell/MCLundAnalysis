#include "src/LundAnalysis.h"
#include <random>
#include <cmath>

//Main program of MCLundAnalysis - Reads Hipo file, calculates kinematics/Ratios, bins (sometimes) and saves
// to root file
int LundAnalysis_single_pion(

                    // hipoFile is the file we read in
                   const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis_pass1/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3052_3.hipo",
                   // rootfile is the file we save data to
                   const char * rootfile = "/work/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/October_13/R1_adjust_flipped_Mki_sign_100k.root"
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
    tree_MC->Branch("M",&M);
    tree_MC->Branch("M_h",&m_1);
    tree_MC->Branch("pT_BF",&pT_BF);
    tree_MC->Branch("x",&x);
    tree_MC->Branch("z",&z_h);
    tree_MC->Branch("Q2",&Q2);
    tree_MC->Branch("Q",&Q);
    tree_MC->Branch("R0",&R0); //initial parton momentum
    tree_MC->Branch("R1",&R1); //Measured in gN frame
    tree_MC->Branch("R2",&R2);
    tree_MC->Branch("R2_adjust",&R2_adjust);
    tree_MC->Branch("qTQ_HF",&qTQ_HF);
    tree_MC->Branch("qT_HF",&qT_HF);
    tree_MC->Branch("qT_from_zN",&qT_from_zN);
    tree_MC->Branch("qT_from_pT_BF",&qT_from_pT_BF);
    tree_MC->Branch("qT_from_pT_PF_z",&qT_from_pT_PF_z);
    tree_MC->Branch("qT_from_pT_HF_z",&qT_from_pT_HF_z);
    tree_MC->Branch("M_ki",&M_ki);
    tree_MC->Branch("M_ki2",&M_ki2);
    tree_MC->Branch("M_kf",&M_kf);
    tree_MC->Branch("R01",&R01);
    tree_MC->Branch("R02",&R02);
    tree_MC->Branch("R03",&R03);
    tree_MC->Branch("delta_k_t",&delta_k_t);
    tree_MC->Branch("ki_t",&ki_t);
    tree_MC->Branch("xi",&xi);
    tree_MC->Branch("zeta",&zeta);
    tree_MC->Branch("z_N",&z_N);
    tree_MC->Branch("z_N_hat",&z_N_hat);
    tree_MC->Branch("x_N",&x_N);
    tree_MC->Branch("x_N_hat",&x_N_hat);
    tree_MC->Branch("theta_ki",&theta_ki);
    tree_MC->Branch("theta_H",&theta_H);
    tree_MC->Branch("theta_deltak",&theta_deltak);
    tree_MC->Branch("kf_x",&kf_x);
    tree_MC->Branch("kf_y",&kf_y);
    tree_MC->Branch("kf_z",&kf_z);
    tree_MC->Branch("kf_E",&kf_E);
    tree_MC->Branch("kf_Plus",&kf_Plus);
    tree_MC->Branch("kf_Minus",&kf_Minus);
    tree_MC->Branch("event_num",&event_num);
    tree_MC->Branch("R1p",&R1p);
    
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
    tree_low->Branch("R01",&R01);
    tree_low->Branch("R02",&R02);
    tree_low->Branch("R03",&R03);
    tree_low->Branch("delta_k_t",&delta_k_t);
    tree_low->Branch("ki_t",&ki_t);
    tree_low->Branch("xi",&xi);
    tree_low->Branch("zeta",&zeta);
    tree_low->Branch("z_N",&z_N);
    tree_low->Branch("x_N",&x_N);
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
    tree_high->Branch("delta_k_t",&delta_k_t);
    tree_high->Branch("ki_t",&ki_t);
    tree_high->Branch("xi",&xi);
    tree_high->Branch("z_N",&z_N);
    tree_high->Branch("x_N",&x_N);
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
    int diff_count = 0;
    
    //now get reference to (unique)ptr for accessing data in loop
    //this will point to the correct place when file changes
    //
    //This line comes from AnalysisWithExtraBanks.C
    auto& c12=chain.C12ref();
    //This value counts the number of events in total that pass the clas12reader cuts
    int event_count = 0;
    
    //Loop over all events in the file that pass proton+electron cuts
    while(chain.Next()==true){
        if(tree_count > 100000) {break;} //Uncomment this line to stop the program after 1000 events, useful for debugging/testing
        event_count += 1;
        event_num = event_count;
        // cout << "\nstarting event #" << event_count << "\n";
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
        
        //Selecting final quark
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
            M = init_target.M();
            
            m_1 = pi1.mass;
            
            q = photon.lv;
            q_calc = init_electron - electron.lv;

            Mx = Mxfunc(q, init_target, pi1.lv);

            cth = cthfunc(electron.px,electron.py,electron.pz);
            Q2 = -(q * q);
            Q = pow(Q2,0.5);
            Q2_calc = Q2func(electron_beam_energy,electron.E,cth); //Momentum transfer

            z_h = (init_target * pi1.lv) / (init_target * q);
            s = sfunc(protonMass, electronMass, electron_beam_energy);
            y = yfunc(electron_beam_energy,electron.E);
            x = xfunc(Q2,s,y);
            
            nu = nufunc(electron_beam_energy,electron.E);
            W = Wfunc(Q2,protonMass,nu);
            
            
            kf = quark.lv;
            ki = kf - q;
            k = kf - q;

            
            
            //
            //// Breit Frame Kinematics
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
            BFAngle = q_BF_no_rot_VectUnit.Angle(zAxis) - pi;
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
             
            //hadron (P_h)
            pi1_BF = pi1.lv;
            pi1_BF.Boost(BFBoost);
            pi1_BF.Rotate(BFAngle,BFAxis);
            
            pT_BF_vect = PtVectfunc(pi1_BF);
            pT_BF = Ptfunc(pi1_BF);
            
            //partons
            ki_BF = ki;
            kf_BF = kf;
            k_BF = k;

            
            ki_BF.Boost(BFBoost);
            kf_BF.Boost(BFBoost);
            k_BF.Boost(BFBoost);
            
            ki_BF.Rotate(BFAngle,BFAxis);
            kf_BF.Rotate(BFAngle,BFAxis);
            k_BF.Rotate(BFAngle,BFAxis);
            
            kiT_BF_vect = PtVectfunc(ki_BF); 
            kfT_BF_vect = PtVectfunc(kf_BF); //BF frame boostkfbT in delta k calculation 
            kfT_BF = Ptfunc(kf_BF);
           

            
            //tempt
            kf_x = kf_BF.X();
            kf_y = kf_BF.Y();
            kf_z = kf_BF.Z();
            kf_E = kf_BF.E();
            kf_Plus = kf_BF.Plus() / pow(2,0.5);
            kf_Minus = kf_BF.Minus() / pow(2,0.5);
            
            
            //
            //// Photon Frame
            //
            
            //Boosting
            PFFrame = q + init_target;
            PFBoost = PFFrame.BoostVector();
            PFBoost = -1 * PFBoost;
            
            //Rotation
            q_PF_no_rot = q;
            q_PF_no_rot.Boost(PFBoost);
            q_PF_no_rot_VectUnit = q_PF_no_rot.Vect().Unit();
            PFAngle = q_PF_no_rot_VectUnit.Angle(zAxis) - pi;
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
            k_PF = k;
            ki_PF.Boost(PFBoost);
            kf_PF.Boost(PFBoost);
            k_PF.Boost(PFBoost);
            ki_PF.Rotate(PFAngle,PFAxis);
            kf_PF.Rotate(PFAngle,PFAxis);
            k_PF.Rotate(PFAngle,PFAxis);
            
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
            HFAngle = pi1_HF_no_rot_VectUnit.Angle(zAxis) - pi;
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
            qT_HF_vect = PtVectfunc(q_HF);
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
            
            //np params
            zeta = pi1_BF.Minus()/ kf_BF.Minus();
            xi = ki_BF.Plus() / proton_BF.Plus();
            
            x_N_hat = x_N / xi;
            z_N_hat = z_N / zeta;
            
            qT_from_zN_vect = -1 * pT_BF_vect / z_N;

            qT_from_pT_BF_vect = pT_BF_vect / z_h;
            qT_from_pT_BF = Ptfunc(qT_from_pT_BF_vect);
            qTQ_from_pT_BF = qT_from_pT_BF / sqrt(Q2);


            //qT / Q for plotting
            qT_from_zN = Ptfunc(qT_from_zN_vect);

            //qT from pT/z
            qT_from_pT_HF_z = pT_HF / z_h;
            //ki, k, and delta k
            deltak_T_vect = kfT_BF_vect - (-1 * z_N_hat * qT_from_zN_vect); 
            

            //Ratios in PF frame
            
            ki2 = abs(ki_BF * ki_BF);
            kf2 = abs(kf_BF * kf_BF);
            deltak2 = abs(deltak_T_vect * deltak_T_vect);
            
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
            M_ki2 = ki_BF * ki_BF;
            M_kf = pow(kf2,0.5);
            delta_k_t = pow(deltak2,0.5);
            
            ki_t = Ptfunc(ki_BF);

            //cout << "ki_BF^2: " << ki_BF * ki_BF << "\n" << "kf_BF^2: " << kf_BF * kf_BF << "\n" << "deltak_T_vect^2: " << deltak_T_vect * deltak_T_vect << "\n"; 


            R0 = R0func(ki_BF, kf_BF, deltak_T_vect, Q2);
            R1 = R1func(pi1_BF, ki_BF, kf_BF);
            R2 = R2func(k_BF, Q2);
            /*
            cout << "\nR2: " << R2 << "\n";
            cout << "k^2: " << k_BF * k_BF << "\n";
            cout << "Q2: " << Q2 << "\n";
            k_BF.Print();
            */
            
            R2_adjust = R2func_adjust(qT_HF, Q2);

            R1p = (pi1_BF * proton_BF) / Q2;
            
            R01 = abs((ki_BF * ki_BF) / Q2);
            R02 = abs((kf_BF * kf_BF) / Q2);
            R03 = abs((deltak_T_vect * deltak_T_vect) / Q2);
            
            if(R1 < -100) {
                low_R1_count++;
            }
            xF = xFpi1;
            
            
            //phi values
            // theta_ki = ki_BF.Phi();
            // theta_H = pi1_BF.Phi();
            // theta_deltak = deltak_T_vect.Phi();

            //phi values diff
            // theta_ki = atan(ki_BF.Py() / ki_BF.Px());
            // theta_H = atan(pi1_BF.Py() / pi1_BF.Px());
            // theta_deltak = atan(deltak_T_vect.Py() / deltak_T_vect.Px());
            
            //phi values diff
            theta_ki = atan2(ki_BF.Py() , ki_BF.Px());
            theta_H = atan2(qT_from_zN_vect.Py() , qT_from_zN_vect.Px());
            theta_deltak = atan2(deltak_T_vect.Py() , deltak_T_vect.Px());

            double R41 = abs(ki * ki / (k * k));
            double R42 = abs(kf * kf / (k * k));
            double R43 = abs(deltak_T_vect * deltak_T_vect / (k * k));
            // cout << "\nevent_num: " << event_num << "\n";
            // cout << "\nR41: " << R41 << "\n";
            // cout << "\nR42: " << R42 << "\n";
            // cout << "\nR43: " << R43 << "\n";
            // cout << "theta_H, theta_qT: " << theta_H << ", " << theta_qT << "\n";
/*
            //thetaki
            double ki_BF_quad_1 = atan(abs(ki_BF.Py()) / abs(ki_BF.Px()));
            if(ki_BF.Px() < 0){
                if(ki_BF.Py() < 0){
                    theta_ki = ki_BF_quad_1 + pi;
                        }
                else {
                        theta_ki = pi - ki_BF_quad_1;
                }
            }
            else if (ki_BF.Py() < 0){
                theta_ki = 2 * pi - ki_BF_quad_1;
            }
            else {
                theta_ki = ki_BF_quad_1;
            }

            //theta_H
            double pi1_BF_quad_1 = atan(abs(pi1_BF.Py()) / abs(pi1_BF.Px()));
            if(pi1_BF.Px() < 0){
                if(pi1_BF.Py() < 0){
                    theta_H = pi1_BF_quad_1 + pi;
                        }
                else {
                        theta_H = pi - pi1_BF_quad_1;
                }
            }
            else if (pi1_BF.Py() < 0){
                theta_H = 2 * pi - pi1_BF_quad_1;
            }
            else {
                theta_H = pi1_BF_quad_1;
            }
            

            //theta_deltak
            double deltak_T_vect_quad_1 = atan(abs(deltak_T_vect.Py()) / abs(deltak_T_vect.Px()));
            if(deltak_T_vect.Px() < 0){
                if(deltak_T_vect.Py() < 0){
                    theta_deltak = deltak_T_vect_quad_1 + pi;
                        }
                else {
                        theta_deltak = pi - deltak_T_vect_quad_1;
                }
            }
            else if (deltak_T_vect.Py() < 0){
                theta_deltak = 2 * pi - deltak_T_vect_quad_1;
            }
            else {
                theta_deltak = deltak_T_vect_quad_1;
            }
            */
            /*
            cout << "\n\ndeltak_T_vect: "; deltak_T_vect.Print();
            cout << "(delta_k_t, theta_deltak): "<< delta_k_t << ", "<< theta_deltak << "\n";
            cout << "deltak_T_vect.Phi(): "<< deltak_T_vect.Phi() << "\n";
            cout << "Px, Py from theta: (" << delta_k_t * cos(theta_deltak) << ", " << delta_k_t * sin(theta_deltak) << ")\n";
            cout << "sin(pi),sin(180): " << sin(pi) << ", " << sin(180) << "\n";

            cout << "\n\n qT_from_zN_vect: "; qT_from_zN_vect.Print();
            cout << "(qT_from_zN,theta_H): ("<<qT_from_zN << ", "<< theta_H << ")\n";
            cout << "qT_from_zN_vect.Phi(): "<< qT_from_zN_vect.Phi() << "\n";
            cout << "Px, Py from theta: (" << qT_from_zN * cos(theta_deltak) << ", " << qT_from_zN * sin(theta_deltak) << ")\n";
            cout << "sin(pi),sin(180): " << sin(pi) << ", " << sin(180) << "\n";
            */
            
            //                                              //
            ////                                          ////
            ////// BEGIN checking light cone components //////
            ////                                          ////
            //                                              //
            // cout << "\n starting pion #" << i << "\n";
/*
            //CHECKING LIGHT CONE FOR Q
            cout << "q_BF: " << "\n";
            q_BF.Print();
            cout << "\nq in Breit Frame (not rotated) LC variables: (+,-) = (" << q_BF_no_rot.Plus()/ pow(2,0.5) << ", " << q_BF_no_rot.Minus()/ pow(2,0.5) << ")\n";
            
            cout << "q in Breit Frame, rotated, LC variables: (+,-) = (" << q_BF.Plus()/ pow(2,0.5) << ", " << q_BF.Minus()/ pow(2,0.5) << ")\n";
            cout << "Q/root(2): " << pow(Q2,0.5) / pow(2,0.5) << '\n';
*/
            //CHECKING LIGHT CONE FOR PROTON
/*
            cout << "\nproton_BF: " << "\n";
            proton_BF.Print();
            cout << "proton in Breit Frame (not rotated) LC variables: (+,-) = (" << proton_BF_no_rot.Plus()/ pow(2,0.5) << ", " << proton_BF_no_rot.Minus()/ pow(2,0.5) << ")\n";
            
            cout << "proton in Breit Frame, rotated, LC variables: (+,-) = (" << proton_BF.Plus()/ pow(2,0.5) << ", " << proton_BF.Minus()/ pow(2,0.5) << ")\n";
            cout << "Plus: Q/(x_N * root(2)): " << pow(Q2,0.5) / (x_N * pow(2,0.5)) << '\n';
            cout << "Minus: (x_N * M^2)/(root(2) * Q): " << (x_N * pow(M,2)) / (pow(Q2,0.5) * pow(2,0.5)) << '\n';
*/
/*
            //CHECKING LIGHT CONE FOR PION

            //use qT_from_zN_vect for qT when paper uses bold q_T
            double M_B2 = pow(m_1,2);
            cout << "\n z_N: " << z_N << "\n";
            double theory_pi_BF_plus = (M_B2 + pow(z_N,2) *(qT_from_zN_vect * qT_from_zN_vect)) / (pow(2,0.5) * z_N * Q);

            
            double theory_pi_BF_plus_angles = pow(2,0.5) * (M_B2 / 2 + pow(z_N,2) *(qT_from_zN_vect * qT_from_zN_vect) * sin(theta_H) * sin(theta_H) / 2 + pow(z_N,2) *(qT_from_zN_vect * qT_from_zN_vect) * cos(theta_H) * cos(theta_H) / 2) / (z_N * Q);

            
            double theory_pi_BF_minus = ((z_N * Q) / (pow(2,0.5)) );
            TVector2 theory_pi_BF_T = -z_N * qT_from_zN_vect;
            double theory_pi_BF_x_angles = -z_N * qT_from_zN * cos(theta_H);
            double theory_pi_BF_y_angles = -z_N * qT_from_zN * sin(theta_H);
            cout << "\npion in Breit Frame, rotated, LC variables: (+,-,x,y) = (" << pi1_BF.Plus()/ pow(2,0.5) << ", " << pi1_BF.Minus()/ pow(2,0.5) << ", " << pi1_BF.Px() << ", " << pi1_BF.Py() << ")\n";
            cout << "Theory (plus,minus,x,y) (" <<theory_pi_BF_plus << ", " << theory_pi_BF_minus << ", " <<  theory_pi_BF_x_angles << ", " << theory_pi_BF_y_angles << ")\n";
*/
/*
            //CHECKING LIGHT CONE FOR ki - note - angles not used for ki
            cout << "ki_T: " << pow(kiT_BF_vect * kiT_BF_vect,0.5) << "\n";
            cout << "M_ki: " << ki_BF * ki_BF << "\n";
            double theory_ki_BF_Plus = Q / (x_N_hat * pow(2,0.5));
            double theory_ki_BF_Minus = (x_N_hat * ((ki_BF * ki_BF) + (kiT_BF_vect * kiT_BF_vect))) / (pow(2,0.5) * Q);
            
            cout << "\nki in Breit Frame, rotated, LC variables: (+,-,x,y) = (" << ki_BF.Plus()/ pow(2,0.5) << ", " << ki_BF.Minus()/ pow(2,0.5) << ", "<< ki_BF.X() << ", "<< ki_BF.Y()<< ")\n";
            cout << "ki theory (plus,minus) (" <<theory_ki_BF_Plus << ", " << theory_ki_BF_Minus << ")\n";
*/
         /*
            //CHECKING LIGHT CONE FOR kf

            //cout << "z_N_hat, Q, M_kf, qT, delta_k_t, theta_H, theta_deltak: " << z_N_hat << ", " << Q << ", " << pow(kf_BF * kf_BF,0.2)<< ", " << qT_from_zN<< ", " << delta_k_t<< ", " << theta_H<< ", " << theta_deltak<< "\n";
            double theory_kf_BF_Plus = (((kf_BF * kf_BF) + (kfT_BF_vect * kfT_BF_vect))) / (z_N_hat * pow(2,0.5) * Q);
            double theory_kf_BF_Plus_w_phi = (1 / (z_N_hat * pow(2,0.5) * Q)) * ((kf_BF * kf_BF) + pow((qT_from_zN * z_N_hat * sin(theta_H) - delta_k_t * sin(theta_deltak)),2) + pow((qT_from_zN * z_N_hat * cos(theta_H) - delta_k_t * cos(theta_deltak)),2));
            double theory_kf_BF_Plus_new_phi = (1 / (z_N_hat * pow(2,0.5) * Q)) * ((kf_BF * kf_BF) + pow((qT_from_zN_vect.Py() * z_N_hat - deltak_T_vect.Py()),2) + pow((qT_from_zN_vect.Px() * z_N_hat - deltak_T_vect.Px()),2));
            double theory_kf_BF_Minus = (z_N_hat * Q) / pow(2,0.5); //same with or without angles
            // cout << "\n\n pion #" << i << "\n";
            // cout << "kf in Breit Frame, rotated, LC variables: (+,-) = (" << kf_BF.Plus()/ pow(2,0.5) << ", " << kf_BF.Minus()/ pow(2,0.5) << ")\n";
            // cout << "kf theory (plus,minus) (" <<theory_kf_BF_Plus << ", " << theory_kf_BF_Minus << ")\n";
            // cout << "kf theory w phi (plus,minus) (" <<theory_kf_BF_Plus_w_phi << ", " << theory_kf_BF_Minus << ")\n";
            
            
                cout << "\n\n pion #" << i << "\n";
                cout << "kf in Breit Frame, rotated, LC variables: (+,-,x,y) = (" << kf_BF.Plus()/ pow(2,0.5) << ", " << kf_BF.Minus()/ pow(2,0.5) <<"," <<kf_BF.X() <<kf_BF.Y() <<  ")\n";
                // cout << "kf theory (plus,minus) (" <<theory_kf_BF_Plus << ", " << theory_kf_BF_Minus << ")\n";
                cout << "kf theory w phi (plus,minus) (" <<theory_kf_BF_Plus_w_phi << ", " << theory_kf_BF_Minus << ")\n";
                // cout << "theta_H: " << theta_H << "\n";
                // cout << "theta_deltak: " << theta_deltak << "\n";
            */
            /*

            //CHECKING LIGHT CONE FOR k (k = kf - q)
            cout << "\n\n";
            double theory_k_BF_x = -qT_from_zN * z_N_hat * cos(theta_H) + delta_k_t * cos(theta_deltak);
            double theory_k_BF_y = -qT_from_zN * z_N_hat * sin(theta_H) + delta_k_t * sin(theta_deltak);
            double theory_k_BF_Plus = (pow(2,0.5) / 2) * Q + pow(2,0.5) * (pow(M_kf,2) + pow((theory_k_BF_y),2) + pow((theory_k_BF_x),2)) / (2 * Q * z_N_hat);
            double theory_k_BF_Minus = pow(2,0.5) * Q * (z_N_hat - 1) / 2;
            cout << "k in Breit Frame, rotated, LC variables: (+,-,x,y) = (" << k_BF.Plus()/ pow(2,0.5) << ", " << k_BF.Minus()/ pow(2,0.5) << ", " << k_BF.Px() << ", "<<k_BF.Py()<<  ")\n";
            cout << "k theory (plus,minus,x,y) (" <<theory_k_BF_Plus << ", " << theory_k_BF_Minus << ", " <<theory_k_BF_x << ", " << theory_k_BF_y <<")\n";
*/
            /*
            //CHECKING LIGHT CONE FOR k based on affinity paper kf (k = kf - q)
            cout << "\n\n";

            double theory_kf_BF_plus = (1 / (pow(2,0.5) * Q * z_N_hat)) * ((kf_BF * kf_BF) + pow((qT_from_zN * z_N_hat * sin(theta_H) - (delta_k_t * sin(theta_deltak))),2) + pow((qT_from_zN * z_N_hat * cos(theta_H) - delta_k_t * cos(theta_deltak)),2));
            
            double theory_k_Plus_long = pow(2,0.5) * (M_kf * M_kf + Q2 * z_N_hat + pow(qT_from_zN,2) * pow(z_N_hat,2) - 2 * qT_from_zN * z_N_hat * delta_k_t * cos(theta_H - theta_deltak) + pow(delta_k_t,2)) / (2 * pow(Q2,0.5) * z_N_hat);
            
            double theory_k_BF_x = -qT_from_zN * z_N_hat * cos(theta_H) + delta_k_t * cos(theta_deltak);
            double theory_k_BF_y = -qT_from_zN * z_N_hat * sin(theta_H) + delta_k_t * sin(theta_deltak);
            double theory_k_BF_Plus_kf_q = theory_kf_BF_plus - q_BF.Plus();
            double theory_k_BF_Plus = (pow(2,0.5) / 2) * Q + pow(2,0.5) * (pow(M_kf,2) + pow((theory_k_BF_y),2) + pow((theory_k_BF_x),2)) / (2 * Q * z_N_hat);
            double theory_k_BF_Minus = pow(2,0.5) * Q * (z_N_hat - 1) / 2;
            cout << "k in Breit Frame, rotated, LC variables: (+,-,x,y) = (" << k_BF.Plus()/ pow(2,0.5) << ", " << k_BF.Minus()/ pow(2,0.5) << ", " << k_BF.Px() << ", "<<k_BF.Py()<<  ")\n";
            cout << "k theory (plus long,minus,x,y) (" <<theory_k_Plus_long << ", " << theory_k_BF_Minus << ", " <<theory_k_BF_x << ", " << theory_k_BF_y <<")\n";
            cout << "kf_T: " << Ptfunc(kf_BF) << "\n";
            //cout << "hadron #" << i << " delta_k_t - qTz_N_hat: " <<  Ptfunc(deltak_T_vect -  qT_from_zN_vect * z_N_hat)<< "\n";
            */
            /*
            //CHECKING R2 against theory
            double R2_theory = (pow(M_kf,2) - z_N_hat * (pow(M_kf,2) + Q2 * z_N_hat - Q2 - pow(qT_from_zN,2) * z_N_hat + 2 * qT_from_zN * delta_k_t * cos(theta_H - theta_deltak)) + pow(delta_k_t,2)) / (Q2 * z_N_hat);
            cout << "R2 from dot product: " << R2 << "\nR2 from theory: " << R2_theory << "\n";
*/
            //                                            //
            ////                                        ////
            ////// END checking light cone components //////
            ////                                        ////
            //   
            
            // cout << "\nki_t: " << ki_t<<  ", M_ki: " << M_ki<<  ", theta_H: " << theta_H << ", theta_ki: " << theta_ki << "\n";
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
            
	// break;
        }

    }
    cout << "\033[0m" << "\033[49m";
    cout << "Final event_count:" << event_count << '\n';
    cout << "Final tree_count: " << tree_count << '\n';
    cout << "Final events passed: " << passed_count << '\n';
    cout << "low_R1_count " << low_R1_count << '\n';
    tree_maxmin->Fill();
    f->Write();
    delete f;
    return 0;
}
