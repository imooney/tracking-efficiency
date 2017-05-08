//  analysis_funcs.cc
//  Created by Isaac Mooney on 5/4/17.

#include "analysis_funcs.h"

namespace analysis {
    //initializes the random seed
    std::mt19937 rands() {
        //seeding pseudorandom number generator
        auto begin = std::chrono::high_resolution_clock::now();
        
        //generating the distribution
        std::random_device rd;
        std::mt19937 g(rd());
        auto end = std::chrono::high_resolution_clock::now();
        g.seed(std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count());
        
        return g;
    }
    
    //Initializes Pythia
    void pythia_init(Pythia8::Pythia & pythia) {
        //Setting parameters for Pythia
        
        pythia.readString("Beams:eCM            = 200.");
        pythia.readString("HardQCD:all          = on");
        pythia.readString("PhaseSpace:pTHatMin  = 20.");
        pythia.readString("PhaseSpace:pTHatMax  = 30.");
        
        //turn off decays
        pythia.readString("111:mayDecay     = off");    //pi0
        pythia.readString("211:mayDecay     = off");    //pi+
        pythia.readString("-211:mayDecay    = off");    //pi-
        pythia.readString("221:mayDecay     = off");    //eta
        pythia.readString("321:mayDecay     = off");    //K+
        pythia.readString("310:mayDecay     = off");    //Kshort
        pythia.readString("130:mayDecay     = off");    //Klong
        pythia.readString("3122:mayDecay    = off");    //Lambda0
        pythia.readString("3212:mayDecay    = off");    //Sigma0
        pythia.readString("3112:mayDecay    = off");    //Sigma-
        pythia.readString("3222:mayDecay    = off");    //Sigma+
        pythia.readString("3312:mayDecay    = off");    //Xi-
        pythia.readString("3322:mayDecay    = off");    //Xi0
        pythia.readString("3334:mayDecay    = off");    //Omega-
        
        pythia.init();
        
        return;
    }
    
    //checks that constituent-level cuts are satisfied
    bool constituent_cuts(const Pythia8::Particle particle) {
        if (fabs(particle.eta()) < eta_cut && particle.pT() >= pt_cut) {
            return true;
        }
        return false;
    }
    
    //checks that pseudojet-level cuts are satisfied
    bool jet_cuts(const fastjet::PseudoJet pseudojet) {
        if (fabs(pseudojet.eta()) < jet_eta) {
            return true;
        }
        return false;
    }
    
    //checks that the tracking efficiency criteria are passed
    bool efficiency_cut(const double effic_num, const Pythia8::Particle particle) {
        if (effic_num > efficiency && particle.isCharged()) {
            //means it won't be considered (potentially slightly confusing)
            return true;
        }
        return false;
    }
    
    void add_particles(std::uniform_real_distribution<> dis, std::mt19937 g, Pythia8::Event event, containers * container) {
        for (unsigned i = 0; i < event.size(); ++i) {
            //only select final state particles within particle-level eta cut
            if (event[i].isFinal() && event[i].isVisible()) {
                //save particle charge info for later in fastjet
                fastjet::PseudoJet current = fastjet::PseudoJet(event[i].px(),event[i].py(),
                                              event[i].pz(),event[i].e());
                //gives the PseudoJet a charge if applicable
                current.set_user_index(event[i].isCharged());
                //fill vector of raw particles, regardless of cuts
                container->uncut_part.push_back(current);
                if (event[i].isCharged()) {
                    container->c_uncut_part.push_back(current);
                }
                if (constituent_cuts(event[i])) {
                    double effic_num = dis(g);
                    if (efficiency_cut(effic_num, event[i])) {
                        //fill the original vector, skip filling the tracking efficiency applied vector
                        container->cut_part.push_back(current);
                        container->cut2_part.push_back(current);
                        container->c_cut_part.push_back(current);
                        container->c_cut2_part.push_back(current);
                        continue;
                    }
                    else {
                        //fill both
                        container->cut_part.push_back(current);
                        container->cut2_part.push_back(current);
                        container->effic_part.push_back(current);
                        if (event[i].isCharged()) {
                            container->c_cut_part.push_back(current);
                            container->c_cut2_part.push_back(current);
                            container->c_effic_part.push_back(current);
                        }
                    }
                }
            }
        }
    }
    
    //void cluster(containers * container, const fastjet::JetDefinition jet_def) {}
    
    //Test to see if the pt spectrum looks reasonable without two highest-pt jets
    void without_lead2_jets(const std::vector<fastjet::PseudoJet> effic_jets, double & without_leadsublead, TTree* ttests) {
        for (unsigned i = 0; i < effic_jets.size(); ++ i) {
            if (jet_cuts(effic_jets[i])) {
                if (i != 0 && i != 1) {
                    without_leadsublead = effic_jets[i].pt();
                    ttests->Fill();
                }
            }
        }
    }
    
    //Pt difference between leading jet after efficiency correction and geometrically
    //closest (eta-phi) jet before correction
    void geometric_pt_diff(const std::vector<fastjet::PseudoJet> effic_jets, const std::vector<fastjet::PseudoJet> cut2_jets, double & ptdiff) {
        double mindist  =   99999;
        double pt_diff   =  -99999;
        if (effic_jets.size() != 0 ) {
            for (unsigned i = 0; i < cut2_jets.size(); ++ i) {
                double dist = effic_jets[0].delta_R(cut2_jets[i]);
                if (dist < mindist) {
                    mindist = dist;
                    pt_diff = effic_jets[0].pt() - cut2_jets[i].pt();
                }
            }
        }
        if (mindist != 99999) {
            ptdiff = pt_diff;
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    //IMPLEMENTATION FOR CLASS 'ANALYSIS::CONTAINERS'!
    
    containers::containers() {
        //trees
        tuncut   = new TTree("uncut","uncut");
        tcuncut  = new TTree("chargedUncut","charged uncut");
        tcut     = new TTree("cut","cut");
        tccut    = new TTree("chargedCut","charged cut");
        tcut2    = new TTree("cut2","cut2");
        tccut2   = new TTree("chargedCut2","charged cut2");
        teffic   = new TTree("efficiency","efficiency");
        tceffic  = new TTree("chargedEfficiency","charged efficiency");
        tlead    = new TTree("leadingJetpTs","leading jet pTs");
        tdiffs   = new TTree("numberAndpTDifferences","number and pT differences");
        ttests   = new TTree("tests","tests");
            
        tuncutcons    = new TTree("uncutConstituents","uncut constituents");
        tcuncutcons   = new TTree("chargedUncutConstituents","charged uncut constituents");
        
        //book space for the vectors of particles
        uncut_part.reserve(50); c_uncut_part.reserve(50);   cut_part.reserve(50);
        cut2_part.reserve(50);  c_cut_part.reserve(50);     c_cut2_part.reserve(50);
        effic_part.reserve(50); c_effic_part.reserve(50);
        
        //jets, too, if I decide to implement them
        /*
        uncut_jets.reserve(50); c_uncut_jets.reserve(50);   cut_jets.reserve(50);
        cut2_jets.reserve(50);  c_cut_jets.reserve(50);     c_cut2_jets.reserve(50);
        effic_jets.reserve(50); c_effic_jets.reserve(50);
        */
    }
        
    //getters & setters
    
    void containers::SetBranches() {
        //branches
        tuncut->Branch("uncut_px", &uncut_px);                      tuncut->Branch("uncut_py", &uncut_py);
        tuncut->Branch("uncut_pz", &uncut_pz);                      tuncut->Branch("uncut_E", &uncut_E);
        tuncut->Branch("uncut_phi", &uncut_phi);                    tuncut->Branch("uncut_eta", &uncut_eta);
        tuncut->Branch("uncut_Pt", &uncut_Pt);
            
        tcuncut->Branch("c_uncut_px", &c_uncut_px);                 tcuncut->Branch("c_uncut_py", &c_uncut_py);
        tcuncut->Branch("c_uncut_pz", &c_uncut_pz);                 tcuncut->Branch("c_uncut_E", &c_uncut_E);
        tcuncut->Branch("c_uncut_phi", &c_uncut_phi);               tcuncut->Branch("c_uncut_eta", &c_uncut_eta);
        tcuncut->Branch("c_uncut_Pt", &c_uncut_Pt);
            
        tuncutcons->Branch("uncut_cons_px", &uncut_cons_px);        tuncutcons->Branch("uncut_cons_py", &uncut_cons_py);
        tuncutcons->Branch("uncut_cons_pz", &uncut_cons_pz);        tuncutcons->Branch("uncut_cons_E", &uncut_cons_E);
        tuncutcons->Branch("uncut_cons_phi", &uncut_cons_phi);      tuncutcons->Branch("uncut_cons_eta", &uncut_cons_eta);
        tuncutcons->Branch("uncut_cons_Pt", &uncut_cons_Pt);
            
        tcuncutcons->Branch("c_uncut_cons_px", &c_uncut_cons_px);   tcuncutcons->Branch("c_uncut_cons_py", &c_uncut_cons_py);
        tcuncutcons->Branch("c_uncut_cons_pz", &c_uncut_cons_pz);   tcuncutcons->Branch("c_uncut_cons_E", &c_uncut_cons_E);
        tcuncutcons->Branch("c_uncut_cons_phi", &c_uncut_cons_phi); tcuncutcons->Branch("c_uncut_cons_eta", &c_uncut_cons_eta);
        tcuncutcons->Branch("c_uncut_cons_Pt", &c_uncut_cons_Pt);
            
        tcut->Branch("cut_px", &cut_px);                            tcut->Branch("cut_py", &cut_py);
        tcut->Branch("cut_pz", &cut_pz);                            tcut->Branch("cut_E", &cut_E);
        tcut->Branch("cut_phi", &cut_phi);                          tcut->Branch("cut_eta", &cut_eta);
        tcut->Branch("cut_Pt", &cut_Pt);
            
        tccut->Branch("c_cut_px", &c_cut_px);                       tccut->Branch("c_cut_py", &c_cut_py);
        tccut->Branch("c_cut_pz", &c_cut_pz);                       tccut->Branch("c_cut_E", &c_cut_E);
        tccut->Branch("c_cut_phi", &c_cut_phi);                     tccut->Branch("c_cut_eta", &c_cut_eta);
        tccut->Branch("c_cut_Pt", &c_cut_Pt);
            
        tcut2->Branch("cut2_px", &cut2_px);                         tcut2->Branch("cut2_py", &cut2_py);
        tcut2->Branch("cut2_pz", &cut2_pz);                         tcut2->Branch("cut2_E", &cut2_E);
        tcut2->Branch("cut2_phi", &cut2_phi);                       tcut2->Branch("cut2_eta", &cut2_eta);
        tcut2->Branch("cut2_Pt", &cut2_Pt);
            
        tccut2->Branch("c_cut2_px", &c_cut2_px);                    tccut2->Branch("c_cut2_py", &c_cut2_py);
        tccut2->Branch("c_cut2_pz", &c_cut2_pz);                    tccut2->Branch("c_cut2_E", &c_cut2_E);
        tccut2->Branch("c_cut2_phi", &c_cut2_phi);                  tccut2->Branch("c_cut2_eta", &c_cut2_eta);
        tccut2->Branch("c_cut2_Pt", &c_cut2_Pt);
            
        teffic->Branch("effic_px", &effic_px);                      teffic->Branch("effic_py", &effic_py);
        teffic->Branch("effic_pz", &effic_pz);                      teffic->Branch("effic_E", &effic_E);
        teffic->Branch("effic_phi", &effic_phi);                    teffic->Branch("effic_eta", &effic_eta);
        teffic->Branch("effic_Pt", &effic_Pt);
            
        tceffic->Branch("c_effic_px", &c_effic_px);                 tceffic->Branch("c_effic_py", &c_effic_py);
        tceffic->Branch("c_effic_pz", &c_effic_pz);                 tceffic->Branch("c_effic_E", &c_effic_E);
        tceffic->Branch("c_effic_phi", &c_effic_phi);               tceffic->Branch("c_effic_eta", &c_effic_eta);
        tceffic->Branch("c_effic_Pt", &c_effic_Pt);
            
        tdiffs->Branch("num_diff_raw", &num_diff_raw);              tdiffs->Branch("c_num_diff_raw", &c_num_diff_raw);
        tdiffs->Branch("num_diff", &num_diff);                      tdiffs->Branch("c_num_diff", &c_num_diff);
        tdiffs->Branch("ptdiff", &ptdiff);
            
        tlead->Branch("uncut_leadPt", &uncut_leadPt);               tlead->Branch("c_uncut_leadPt", &c_uncut_leadPt);
        tlead->Branch("cut_leadPt", &cut_leadPt);                   tlead->Branch("c_cut_leadPt", &c_cut_leadPt);
        tlead->Branch("cut2_leadPt", &cut2_leadPt);                 tlead->Branch("c_cut2_leadPt", &c_cut2_leadPt);
        tlead->Branch("effic_leadPt", &effic_leadPt);               tlead->Branch("c_effic_leadPt", &c_effic_leadPt);
            
        ttests->Branch("without_leadsublead", &without_leadsublead);
        
    }
    
    //void containers::SetJets();
    
    void containers::Clear() {
        uncut_part.clear();     cut_part.clear();
        cut2_part.clear();      effic_part.clear();
        c_uncut_part.clear();   c_cut_part.clear();
        c_cut2_part.clear();    c_effic_part.clear();
    }
    
    void containers::write() {
        tuncut->Write();  tcuncut->Write(); tcut->Write();    tccut->Write();
        tcut2->Write();   tccut2->Write();  teffic->Write();  tceffic->Write();
        tlead->Write();   tdiffs->Write();  ttests->Write();
        
        tuncutcons->Write(); tcuncutcons->Write();
    }
}
