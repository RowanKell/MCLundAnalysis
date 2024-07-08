#include "src/LundAnalysis.h"
#include <random>
#include <cmath>

//Main program of MCLundAnalysis - Reads Hipo file, calculates kinematics/Ratios, bins (sometimes) and saves
// to root file
int LundAnalysis_single_pion(

                    // hipoFile is the file we read in
                   const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3051_0.hipo",
                   // rootfile is the file we save data to
                   const char * rootfile = "/work/clas12/users/rojokell/MCLundAnalysis/OutputFiles/Files_Spring_24/July_8/slurm_test.root"
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
    tree_MC->Branch("pT",&pt_gN);
    tree_MC->Branch("Q2",&Q2);
    tree_MC->Branch("R0",&R0); //initial parton momentum
    tree_MC->Branch("R1",&R1); //Measured in gN frame
    tree_MC->Branch("R2",&R2);
    tree_MC->Branch("q_TdivQ",&qTQ);
    tree_MC->Branch("qTQ_hadron",&qTQ_hadron);
    tree_MC->Branch("R2_adjust",&R2_adjust);
    tree_MC->Branch("qT_calc", &qT_calc);
    tree_MC->Branch("qT_diff", &qT_diff);
    tree_MC->Branch("qT_hadron_mag", &qT_hadron_mag);
    tree_MC->Branch("qTQ_calc", &qTQ_calc);
    
    TTree *tree_maxmin = new TTree("tree_maxmin","Tree with max and min Ri values");
    
    tree_maxmin->Branch("R0_max",&R0_max); //initial parton momentum
    tree_maxmin->Branch("R1_max",&R1_max); //Measured in gN frame
    tree_maxmin->Branch("R2_max",&R2_max);
    
    tree_maxmin->Branch("R0_min",&R0_min); //initial parton momentum
    tree_maxmin->Branch("R1_min",&R1_min); //Measured in gN frame
    tree_maxmin->Branch("R2_min",&R2_min);
    
    TTree *tree_test = new TTree("tree_test","Tree with kinematics for viewing distributions");
    tree_test->Branch("delta_k_T",&delta_k_T);
    tree_test->Branch("M_ki",&M_ki);
    tree_test->Branch("M_kf",&M_kf);
    
    tree_test->Branch("r_delta_k_T",&r_delta_k_T);
    tree_test->Branch("r_M_ki",&r_M_ki);
    tree_test->Branch("r_M_kf",&r_M_kf);
    
//     tree_test->Branch("ki_x",&kix);
//     tree_test->Branch("ki_y",&kiy);
//     tree_test->Branch("ki_z",&kiz);
//     tree_test->Branch("kf_x",&kfx);
//     tree_test->Branch("kf_y",&kfy);
//     tree_test->Branch("kf_z",&kfz);
    
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
//         if(tree_count > 1000) {break;} //Uncomment this line to stop the program after 1000 events, useful for debugging/testing
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
            //Fill first pion particle object
            pi1.fillParticle(pi_v.v_id[i], pi_v.v_pid[i], pi_v.v_px[i], pi_v.v_py[i],pi_v.v_pz[i], pi_v.v_daughter[i], pi_v.v_parent[i], pi_v.v_mass[i], pi_v.v_vz[i]); 
            pi1.setVectors();


            quark.fillParticle(quark.v_id[quark.final_id], quark.v_pid[quark.final_id], quark.v_px[quark.final_id], quark.v_py[quark.final_id], quark.v_pz[quark.final_id], quark.v_daughter[quark.final_id], quark.v_parent[quark.final_id], quark.v_mass[quark.final_id], quark.v_vz[quark.final_id]);
            quark.setVectors();
            
            //Setting inital beam and target particles
            init_electron.SetPxPyPzE(0, 0, sqrt(electron_beam_energy * electron_beam_energy - electronMass * electronMass), electron_beam_energy);
            init_target.SetPxPyPzE(0, 0, 0, proton.E);

            m_1 = pi1.mass;
            
            q = photon.lv;
            q_calc = init_electron - electron.lv;
            qx = q_calc.Px();
            qy = q_calc.Py();
            qz = q_calc.Pz();
            qE = q_calc.E();

            photonx = photon.lv.Px();
            photony = photon.lv.Py();
            photonz = photon.lv.Pz();
            photonE = photon.lv.E();


            //Missing mass
            Mx = Mxfunc(q, init_target, pi1.lv);

            cth = cthfunc(electron.px,electron.py,electron.pz);
            Q2 = -(q * q);
            Q2_calc = Q2func(electron_beam_energy,electron.E,cth); //Momentum transfer

            z_h = (init_target * pi1.lv) / (init_target * q);
            s = sfunc(protonMass, electronMass, electron_beam_energy);
            y = yfunc(electron_beam_energy,electron.E);
            x = Q2/s/y; // Bjorken x

            kf = quark.lv;
            ki = kf - q;
            //Cut Kinematics

            gN = q;
            gN += init_target;
            gNBoost = gN.BoostVector();
            gNBoostNeg = -gNBoost;

            lv_p1_gN = pi1.lv;
            lv_p1_gN.Boost(gNBoostNeg);

            lv_q_gN = q;
            lv_q_gN.Boost(gNBoostNeg);

            //Need pion pt in gN frame
            pt_gN = Ptfunc(lv_p1_gN);
            
            target_gN = init_target;
            target_gN.Boost(gNBoostNeg);

            //Need partonic in gN
            ki_gN = ki;
            ki_gN.Boost(gNBoostNeg);

            kf_gN = kf;
            kf_gN.Boost(gNBoostNeg);

            //nu and W
            nu = nufunc(electron_beam_energy,electron.E);
            W = Wfunc(Q2,protonMass,nu);

            //Feynman x
            xFpi1 = xFfunc(lv_p1_gN,lv_q_gN,W);

            // Breit Frame Kinematics for delta k
            Breit = q;
            Breit_target.SetPxPyPzE(0,0,0,2 * x *protonMass); // E^2 = M^2 + P^2 --> P = 0 so E = M = 2 * x * protonmass
            Breit += Breit_target;
            BreitBoost = Breit.BoostVector();
            BreitBoost = -1 * BreitBoost;

            //Setting up delta k variables
            kfBreit = kf;
            kfBreit.Boost(BreitBoost);
            kfBreitTran = PtVectfunc(kfBreit); //breit frame boostkfbT in delta k calculation - needs to be a transverse light cone vector of form (V_x, V_y)

            q_Breit = q;
            q_Breit.Boost(BreitBoost);
            proton_Breit = proton.lv;
            proton_Breit.Boost(BreitBoost);

            Breit1 = pi1.lv;
            Breit1.Boost(BreitBoost);
            BreitTran1 = PtVectfunc(Breit1);
            
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

            //Hadron frame kinematics for q_T
            //Note: hadron frame here refers to frame where the target proton and outgoing hadron are back to back
            //and the target proton has the same four momentum as it has in the breit frame
            Hadron_frame = pi1.lv;
            Hadron_frame += Breit_target;
            HadronBoost = Hadron_frame.BoostVector();
            HadronBoost = -1 * HadronBoost;
            q_hadron = q;
            q_hadron.Boost(HadronBoost);
            q_T_hadron = PtVectfunc(q_hadron);

            q_T_lab = PtVectfunc(q);

            qTQ_lab = Ptfunc(q_T_lab) / sqrt(Q2);
            qTQ_hadron = Ptfunc(q_T_hadron) / sqrt(Q2);

            //Hadron frame kinematics for q_T of pi_1
            pion_frame = pi1.lv;
            pion_frame += Breit_target;
            pionBoost = pion_frame.BoostVector();
            pionBoost = -1 * pionBoost;
            q_pion = q;
            q_pion.Boost(pionBoost);
            q_T_pion = PtVectfunc(q_pion);
            qTQ_pion = Ptfunc(q_T_pion) / sqrt(Q2);


            //
            //Photon Frame
            //

            PF1 = pi1.lv;
            PF1.Boost(PFBoost);
            PF1.Rotate(PFAngle,PFAxis);
            PF1Minus = LightConeMinus(PF1);
            
            //Virtual Photon
            qPF.Rotate(PFAngle,PFAxis);
            qPFMinus = LightConeMinus(qPF);
            
            //z_N and q_T
            z_N = PF1Minus / qPFMinus;
            
            q_T = -1 * BreitTran1 / z_N;

            q_T_frac = BreitTran1 / z_h;
            qTQfrac = Ptfunc(q_T_frac) / sqrt(Q2);


            //q_T / Q for plotting
            qTQ = Ptfunc(q_T) / sqrt(Q2);

            //qT from pT/z
            qT_calc = pt_gN / z_h;
            qT_hadron_mag = Ptfunc(q_T_hadron);
            qTQ_calc = qT_calc / sqrt(Q2);
            
            // diff between qT calc and qT hadron
            qT_diff = qT_calc - qT_hadron_mag;
            
            //ki, k, and delta k
            deltak = kfBreitTran - (-1 * z_N * q_T); 

            k = kf - q;
            k_gN = k;
            k_gN.Boost(gNBoostNeg);

            //Ratios in gN frame
            R0 = R0func(ki_gN, kf_gN, deltak, Q2);
            
            ki2 = abs(ki_gN * ki_gN);
            kf2 = abs(kf_gN * kf_gN);
            deltak2 = abs(deltak * deltak);
            
            kix = ki_gN.Px();
            kiy = ki_gN.Py();
            kiz = ki_gN.Pz();
            kfx = kf_gN.Px();
            kfy = kf_gN.Py();
            kfz = kf_gN.Pz();
            
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
            
            R1 = R1func(lv_p1_gN, ki_gN, kf_gN);
            if(R1 < -100) {
                low_R1_count++;
            }
            ki_Breit = ki;
            ki_Breit.Boost(BreitBoost);
            kf_Breit = kf;
            kf_Breit.Boost(BreitBoost);
            k_Breit = k;
            k_Breit.Boost(BreitBoost);

            R2 = R2func(k_gN, Q2);
            R2_adjust = R2func_adjust(Ptfunc(q_T_hadron), Q2);
            xF = xFpi1;


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
                if(qTQ_pion < qTQcut_min || qTQ_pion > qTQcut_max){continue;}
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
            tree_test->Fill();
            //zbins:
            for(int i = 0; i < zbins.size(); i++) {
                if(z_h <= zbins[i]) {
                    zbinv[i].zFillVectors(x, Q2, pt_gN, R0, R1, R2_adjust);
//                     cout << "filling z_h bin #" << i << " with: x - " << x << " | pt_gN: " << pt_gN << "\n";
                    zcounts[i] += 1;
                    break;
                }
            }
            //x bins
            for(int i = 0; i < xbins.size(); i++) {
                if(x <= xbins[i]) {
                    xbinv[i].xFillVectors(z_h, Q2, pt_gN, R0, R1, R2_adjust);
                    xcounts[i] += 1;
                    break;
                }
            }
            //Q2 bins
            for(int i = 0; i < Q2bins.size(); i++) {
                if(Q2 <= Q2bins[i]) {
                    Q2binv[i].Q2FillVectors(x, z_h, pt_gN, R0, R1, R2_adjust);
                    break;
                }
            }
            if(qTQ_pion > qTQmax){qTQmax = qTQ_pion;}
            //qTQ bins
            for(int i = 0; i < qTQbins.size(); i++) {
                if(qTQ_pion <= qTQbins[i]) {
                    qTQbinv[i].qTQFillVectors(x, z_h, Q2, pt_gN, R0, R1, R2_adjust);
                    qTQcounts[i] += 1;
                    break;
                }
            }
        }
	

    }
    cout << "\033[0m" << "\033[49m";
    cout << "Final event_count:" << event_count << '\n';
    cout << "Final tree_count: " << tree_count << '\n';
    cout << "Final events passed: " << passed_count << '\n';
    cout << "low_R1_count " << low_R1_count << '\n';


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

    tree_maxmin->Fill();
    f->Write();
    delete f;
    return 0;
}