#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>

#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TMath.h"

#include <TGraph.h>
#include <TAxis.h>
#include <TPaveText.h>

struct Particle {
    int pdgId;
    int status;
    int mother;
    double px, py, pz, E;
};

struct Lepton {
    TLorentzVector v;
    int charge;
    int flavor;
    int role;  // 0 = from WR, 1 = from N
};

bool comparePt(const TLorentzVector& a, const TLorentzVector& b) {
    return a.Pt() > b.Pt();
}

void read_KS_events_scan() {
    std::cout << "User Test 027:" << std::endl;
    std::cout << "Analyzing KS process: pp->W_R->lN->lljj" << std::endl;
    
    const char* absFigPath = "/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/";
    TString baseFigPath = absFigPath;
    
    gSystem->mkdir(baseFigPath + "KS_Figures_FullLRSM2", kTRUE);
    gSystem->mkdir(baseFigPath + "KS_Figures_FullLRSM2/leptons", kTRUE);
    gSystem->mkdir(baseFigPath + "KS_Figures_FullLRSM2/jets", kTRUE);
    gSystem->mkdir(baseFigPath + "KS_Figures_FullLRSM2/reconstruction", kTRUE);
    
    std::ofstream summary("/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/KS_Figures_FullLRSM2/summary.txt");
    summary << std::fixed << std::setprecision(2);
    
    // Define the Event Files
    std::vector<std::pair<std::string, std::string>> eventFiles = {
        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN100_MWR3000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=100GeV, m_{W_{R}}=3TeV"},
        // // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN100_MWR3500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=100GeV, m_{W_{R}}=3.5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN100_MWR4000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=100GeV, m_{W_{R}}=4TeV"},
        // // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN100_MWR4500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=100GeV, m_{W_{R}}=4.5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN100_MWR5000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=100GeV, m_{W_{R}}=5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN100_MWR6000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=100GeV, m_{W_{R}}=6TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN100_MWR7000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=100GeV, m_{W_{R}}=7TeV"},

        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR3000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=500GeV, m_{W_{R}}=3TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN500_MWR3500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=500GeV, m_{W_{R}}=3.5TeV"},
        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR4000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=500GeV, m_{W_{R}}=4TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN500_MWR4500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=500GeV, m_{W_{R}}=4.5TeV"},
        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR5000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=500GeV, m_{W_{R}}=5TeV"},
        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR6000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=500GeV, m_{W_{R}}=6TeV"},
        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR7000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=500GeV, m_{W_{R}}=7TeV"},

        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN1000_MWR3000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1TeV, m_{W_{R}}=3TeV"},
        // // // // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN1000_MWR3500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1TeV, m_{W_{R}}=3.5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN1000_MWR4000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1TeV, m_{W_{R}}=4TeV"},
        // // // // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN1000_MWR4500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1TeV, m_{W_{R}}=4.5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN1000_MWR5000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1TeV, m_{W_{R}}=5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN1000_MWR6000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1TeV, m_{W_{R}}=6TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN1000_MWR7000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1TeV, m_{W_{R}}=7TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN1500_MWR3000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1.5TeV, m_{W_{R}}=3TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN1500_MWR3500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1.5TeV, m_{W_{R}}=3.5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN1500_MWR4000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1.5TeV, m_{W_{R}}=4TeV"},
        // // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN1500_MWR4500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1.5TeV, m_{W_{R}}=4.5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN1500_MWR5000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1.5TeV, m_{W_{R}}=5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN1500_MWR6000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=1.5TeV, m_{W_{R}}=6TeV"},

        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN2000_MWR3000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=2TeV, m_{W_{R}}=3TeV"},
        // // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN2000_MWR3500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=2TeV, m_{W_{R}}=3.5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN2000_MWR4000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=2TeV, m_{W_{R}}=4TeV"},
        // // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN2000_MWR4500/Events/run_01/unweighted_events.lhe.gz", "m_{N}=2TeV, m_{W_{R}}=4.5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN2000_MWR5000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=2TeV, m_{W_{R}}=5TeV"},
        // {"/storage/alexanum/WORKBENCH/runs/KS_WR_MN2000_MWR6000/Events/run_01/unweighted_events.lhe.gz", "m_{N}=2TeV, m_{W_{R}}=6TeV"}
    };
    
    // Create ONE set of histograms that will be filled from all files
    TH1F *h_lep_invmass = new TH1F("h_lep_invmass", ";M_{ll} [GeV];Events", 100, 0, 3000);
    TH1F *h_lep1_pt = new TH1F("h_lep1_pt", ";p_{T}(l1) [GeV]; Events", 100, 0, 4000);
    TH1F *h_lep2_pt = new TH1F("h_lep2_pt", ";p_{T}(l2) [GeV]; Events", 100, 0, 2000);
    TH1F *h_lep_pt = new TH1F("h_lep_pt", ";p_{T}(l) [GeV]; Events", 100, 0, 4000);
    TH1F *h_lep_deltaR = new TH1F("h_lep_deltaR", ";#DeltaR(ll); Events", 100, 0, 5);
    TH1F *h_jet_pt = new TH1F("h_jet_pt", ";p_{T}(j) [GeV]; Events", 100, 0, 2000);
    TH1F *h_jj_mass = new TH1F("h_jj_mass", ";M(jj) [GeV]; Events", 100, 0, 2000);
    TH1F *h_jj_deltaR = new TH1F("h_jj_deltaR", ";#DeltaR(jj); Events", 100, 0, 5);
    TH1F *h_l2jj_deltaR = new TH1F("h_l2jj_deltaR", ";#DeltaR(l_{2}jj); Events", 100, 0, 5);
    TH1F *h_WR_mass = new TH1F("h_WR_mass", ";m_{W_{R}} [GeV]; Events", 100, 0, 8000);
    TH1F *h_N1_mass = new TH1F("h_N1_mass", ";m_{N} [GeV]; Events", 100, 0, 3000);
    TH1F *h_ll_charge = new TH1F("h_ll_charge", ";Dilepton charge; Events", 2, 0, 2);
    TH1F *h_ll_flavor = new TH1F("h_ll_flavor", ";Dilepton flavor; Events", 6, 0, 6);
    
    h_ll_charge->GetXaxis()->SetBinLabel(1, "OS");
    h_ll_charge->GetXaxis()->SetBinLabel(2, "SS");
    h_ll_flavor->GetXaxis()->SetBinLabel(1, "ee");
    h_ll_flavor->GetXaxis()->SetBinLabel(2, "#mu#mu");
    h_ll_flavor->GetXaxis()->SetBinLabel(3, "#tau#tau");
    h_ll_flavor->GetXaxis()->SetBinLabel(4, "e#mu");
    h_ll_flavor->GetXaxis()->SetBinLabel(5, "e#tau");
    h_ll_flavor->GetXaxis()->SetBinLabel(6, "#mu#tau");
    
    // For overlaid plots
    std::vector<TH1F*> h_lep_invmass_vec;
    std::vector<TH1F*> h_lep1_pt_vec;
    std::vector<TH1F*> h_lep2_pt_vec;
    std::vector<TH1F*> h_lep_pt_vec;
    std::vector<TH1F*> h_lep_deltaR_vec;
    std::vector<TH1F*> h_jet_pt_vec;
    std::vector<TH1F*> h_jj_mass_vec;
    std::vector<TH1F*> h_jj_deltaR_vec;
    std::vector<TH1F*> h_l2jj_deltaR_vec;
    std::vector<TH1F*> h_WR_mass_vec;
    std::vector<TH1F*> h_N1_mass_vec;
    std::vector<TH1F*> h_ll_charge_vec;
    std::vector<TH1F*> h_ll_flavor_vec;

    std::vector<int> colors = {kRed, kBlue, kGreen, kBlack, kViolet, kBlue, kGreen, kBlack};
    std::vector<int> wrStyles = {1, 1, 1, 1, 1, 3, 3, 3, 3};
    // std::vector<int> colors = {kRed, kRed+1, kRed+2, kRed+3, kGreen, kGreen+1, kGreen+2, kGreen+3, kBlue, kBlue+1, kBlue+2, kBlue+3};
    
    // Create histograms for each dataset
    for (size_t i = 0; i < eventFiles.size(); i++) {
        h_lep_invmass_vec.push_back(new TH1F(Form("h_lep_invmass_%zu", i), ";m(#it{ll}) [GeV]; Events", 100, 0, 7000));
        h_lep1_pt_vec.push_back(new TH1F(Form("h_lep1_pt_%zu", i), ";p_{T}(#it{l_{1}}) [GeV]; Events", 100, 0, 4000));
        h_lep2_pt_vec.push_back(new TH1F(Form("h_lep2_pt_%zu", i), ";p_{T}(#it{l_{2}}) [GeV]; Events", 100, 0, 2000));
        h_lep_pt_vec.push_back(new TH1F(Form("h_lep_pt_%zu", i), ";p_{T}(#it{l_{i}}) [GeV]; Events", 100, 0, 4000));
        h_lep_deltaR_vec.push_back(new TH1F(Form("h_lep_deltaR_%zu", i), ";#DeltaR(#it{ll}); Events", 100, 0, 5));
        h_jet_pt_vec.push_back(new TH1F(Form("h_jet_pt_%zu", i), ";p_{T}(#it{j}) [GeV]; Events", 100, 0, 2000));
        h_jj_mass_vec.push_back(new TH1F(Form("h_jj_mass_%zu", i), ";m(#it{jj}) [GeV]; Events", 100, 0, 2000));
        h_jj_deltaR_vec.push_back(new TH1F(Form("h_jj_deltaR_%zu", i), ";#DeltaR(#it{jj}); Events", 100, 0, 5));
        h_l2jj_deltaR_vec.push_back(new TH1F(Form("h_l2jj_deltaR_%zu", i), ";#DeltaR(#it{l_{2}jj}); Events", 100, 0, 5));
        h_WR_mass_vec.push_back(new TH1F(Form("h_WR_mass_%zu", i), ";m_{W_{R}} [GeV]; Events", 100, 0, 8000));
        h_N1_mass_vec.push_back(new TH1F(Form("h_N1_mass_%zu", i), ";m_{N} [GeV]; Events", 100, 0, 3000));
        h_ll_charge_vec.push_back(new TH1F(Form("h_ll_charge_%zu", i), ";Dilepton charge; Events", 2, 0, 2));
        h_ll_flavor_vec.push_back(new TH1F(Form("h_ll_flavor_%zu", i), ";Dilepton flavor; Events", 6, 0, 6));
        
        // Set colors
        h_lep_invmass_vec[i]->SetLineColor(colors[i]);
        h_lep_invmass_vec[i]->SetLineStyle(wrStyles[i]);
        h_lep1_pt_vec[i]->SetLineColor(colors[i]);
        h_lep1_pt_vec[i]->SetLineStyle(wrStyles[i]);
        // h_lep1_pt_vec[i]->SetLogy();
        h_lep2_pt_vec[i]->SetLineColor(colors[i]);
        h_lep2_pt_vec[i]->SetLineStyle(wrStyles[i]);
        h_lep_pt_vec[i]->SetLineColor(colors[i]);
        h_lep_pt_vec[i]->SetLineStyle(wrStyles[i]);
        h_lep_deltaR_vec[i]->SetLineColor(colors[i]);
        h_lep_deltaR_vec[i]->SetLineStyle(wrStyles[i]);
        h_jet_pt_vec[i]->SetLineColor(colors[i]);
        h_jet_pt_vec[i]->SetLineStyle(wrStyles[i]);
        h_jj_mass_vec[i]->SetLineColor(colors[i]);
        h_jj_mass_vec[i]->SetLineStyle(wrStyles[i]);
        h_jj_deltaR_vec[i]->SetLineColor(colors[i]);
        h_jj_deltaR_vec[i]->SetLineStyle(wrStyles[i]);
        h_l2jj_deltaR_vec[i]->SetLineColor(colors[i]);
        h_l2jj_deltaR_vec[i]->SetLineStyle(wrStyles[i]);
        h_WR_mass_vec[i]->SetLineColor(colors[i]);
        h_WR_mass_vec[i]->SetLineStyle(wrStyles[i]);
        h_N1_mass_vec[i]->SetLineColor(colors[i]);
        h_N1_mass_vec[i]->SetLineStyle(wrStyles[i]);
        h_ll_charge_vec[i]->SetLineColor(colors[i]);
        h_ll_charge_vec[i]->SetLineStyle(wrStyles[i]);
        h_ll_flavor_vec[i]->SetLineColor(colors[i]);
        h_ll_flavor_vec[i]->SetLineStyle(wrStyles[i]);
        
        // Set bin labels for categorical histograms
        h_ll_charge_vec[i]->GetXaxis()->SetBinLabel(1, "OS");
        // h_ll_charge_vec[i]->SetMaximum(5000);
        h_ll_charge_vec[i]->GetXaxis()->SetBinLabel(2, "SS");
        h_ll_flavor_vec[i]->GetXaxis()->SetBinLabel(1, "ee");
        h_ll_flavor_vec[i]->GetXaxis()->SetBinLabel(2, "#mu#mu");
        h_ll_flavor_vec[i]->GetXaxis()->SetBinLabel(3, "#tau#tau");
        h_ll_flavor_vec[i]->GetXaxis()->SetBinLabel(4, "e#mu");
        h_ll_flavor_vec[i]->GetXaxis()->SetBinLabel(5, "e#tau");
        h_ll_flavor_vec[i]->GetXaxis()->SetBinLabel(6, "#mu#tau");
    }
    
    // Process each event file
    for (size_t fileIdx = 0; fileIdx < eventFiles.size(); fileIdx++) {
        int oddEvents=0;
        std::cout << "\nProcessing: " << eventFiles[fileIdx].second << std::endl;
        std::cout << "File: " << eventFiles[fileIdx].first << std::endl;
        
        const char* filename = eventFiles[fileIdx].first.c_str();
        TString cmd = TString::Format("gunzip -c %s", filename);
        FILE* pipe = gSystem->OpenPipe(cmd, "r");
        
        if (!pipe) {
            std::cerr << "ERROR: Cannot open LHE file: " << filename << std::endl;
            continue;
        }
        
        char line[4096];
        bool inEvent = false, firstLine = false;
        std::vector<Particle> particles;
        int totalEvents = 0, successfulReco = 0, WL_instances = 0;
        
        while (fgets(line, sizeof(line), pipe)) {
            std::string s(line);
            s.erase(0, s.find_first_not_of(" \t\n\r"));
            s.erase(s.find_last_not_of(" \t\n\r") + 1);
            if (s.empty() || s[0] == '#') continue;
            
            if (s.find("<event>") != std::string::npos) {
                inEvent = true;
                firstLine = true;
                particles.clear();
                continue;
            }
            
            if (s.find("</event>") != std::string::npos) {
                totalEvents++;
                
                std::vector<Lepton> leptons;
                std::vector<TLorentzVector> jets;
                
                for (auto &p : particles) {
                    if (p.status != 1) {
                        if (abs(p.pdgId)==24) {WL_instances++;}
                        continue;
                    }
                    TLorentzVector v(p.px, p.py, p.pz, p.E);
                    int apid = abs(p.pdgId);
                    
                    if (apid == 11 || apid == 13 || apid == 15) {
                        int role = -1;
                        if (p.mother == 3) role = 0; // l1 (from WR)
                        if (p.mother == 4) role = 1; // l2 (from N)

                        leptons.push_back({v, p.pdgId > 0 ? 1 : -1, apid, role}); //pdgID=11->e^- pdgID=-11->e^+  
                        h_lep_pt_vec[fileIdx]->Fill(v.Pt());
                        h_lep_pt->Fill(v.Pt());
                    }
                    
                    if (apid >= 1 && apid <= 6) {   //included 6
                        jets.push_back(v);
                        h_jet_pt_vec[fileIdx]->Fill(v.Pt());
                        h_jet_pt->Fill(v.Pt());
                    }
                }

                // Dilepton properties
                const Lepton* l1 = nullptr; // from WR
                const Lepton* l2 = nullptr; // from N
                if (leptons.size() >= 2) {
                    // std::sort(leptons.begin(), leptons.end(),
                    //     [](const Lepton& a, const Lepton& b) { return a.v.Pt() > b.v.Pt(); });
        
                    for (const auto& l : leptons) {
                        if (l.role == 0) l1 = &l;
                        if (l.role == 1) l2 = &l;
                    }
                    if (!l1 || !l2) continue; // Safety: skip event if assignment failed
                    h_lep1_pt_vec[fileIdx]->Fill(l1->v.Pt());
                    h_lep2_pt_vec[fileIdx]->Fill(l2->v.Pt());
                    h_lep1_pt->Fill(l1->v.Pt());
                    h_lep2_pt->Fill(l2->v.Pt());
                    
                    TLorentzVector ll = l1->v + l2->v;
                    h_lep_invmass_vec[fileIdx]->Fill(ll.M());
                    h_lep_deltaR_vec[fileIdx]->Fill(l1->v.DeltaR(l2->v));
                    h_lep_invmass->Fill(ll.M());
                    h_lep_deltaR->Fill(l1->v.DeltaR(l2->v));
                    
                    bool ss = l1->charge * l2->charge > 0;
                    h_ll_charge_vec[fileIdx]->Fill(ss ? 1 : 0);
                    h_ll_charge->Fill(ss ? 1 : 0);
                    
                    int f1 = l1->flavor;
                    int f2 = l2->flavor;
                    if ((f1 == 11 && f2 == 11)) h_ll_flavor_vec[fileIdx]->Fill(0), h_ll_flavor->Fill(0);
                    else if ((f1 == 13 && f2 == 13)) h_ll_flavor_vec[fileIdx]->Fill(1), h_ll_flavor->Fill(1);
                    else if ((f1 == 15 && f2 == 15)) h_ll_flavor_vec[fileIdx]->Fill(2), h_ll_flavor->Fill(2);
                    else if ((f1 == 11 && f2 == 13) || (f1 == 13 && f2 == 11)) h_ll_flavor_vec[fileIdx]->Fill(3), h_ll_flavor->Fill(3);
                    else if ((f1 == 11 && f2 == 15) || (f1 == 15 && f2 == 11)) h_ll_flavor_vec[fileIdx]->Fill(4), h_ll_flavor->Fill(4);
                    else if ((f1 == 13 && f2 == 15) || (f1 == 15 && f2 == 13)) h_ll_flavor_vec[fileIdx]->Fill(5), h_ll_flavor->Fill(5);
                }
                
                // Dijet properties
                if (jets.size() >= 2) {
                    std::sort(jets.begin(), jets.end(), comparePt);
                    h_jj_mass_vec[fileIdx]->Fill((jets[0] + jets[1]).M());
                    h_jj_deltaR_vec[fileIdx]->Fill(jets[0].DeltaR(jets[1]));
                    h_jj_mass->Fill((jets[0] + jets[1]).M());
                    h_jj_deltaR->Fill(jets[0].DeltaR(jets[1]));
                }
                
                // Reconstruction and DR2(l2jj)
                if (leptons.size() >= 2 && jets.size() >= 2 && l1 && l2) {
                    successfulReco++;
                    TLorentzVector N1 = l2->v + jets[0] + jets[1];
                    TLorentzVector WR = l1->v + N1;
                    h_N1_mass_vec[fileIdx]->Fill(N1.M());
                    h_WR_mass_vec[fileIdx]->Fill(WR.M());
                    h_N1_mass->Fill(N1.M());
                    h_WR_mass->Fill(WR.M());

                    TLorentzVector dijet = jets[0] + jets[1];
                    h_l2jj_deltaR_vec[fileIdx]->Fill(l2->v.DeltaR(dijet));
                    h_l2jj_deltaR->Fill(l2->v.DeltaR(dijet));
                }
                
                inEvent = false;
                continue;
            }
            
            if (inEvent) {
                if (firstLine) { 
                    firstLine = false; 
                    int nParticles = 0;
                    std::istringstream iss(s);
                    iss >> nParticles;  // first number is NUP
                    if (nParticles != 8) {  // 8 = normal, 7 or 9 = unusual
                        oddEvents++;
                        // summary << "Event " << totalEvents+1 << ": unusual NUP = " << nParticles << "\n";
                    }
                    continue;
                }
                if (s[0] == '<') continue;
                
                Particle p;
                int m1, m2, c1, c2;
                double m, l, sp;
                std::istringstream iss(s);
                
                if (iss >> p.pdgId >> p.status >> p.mother >> m2 >> c1 >> c2
                    >> p.px >> p.py >> p.pz >> p.E >> m >> l >> sp) {
                    particles.push_back(p);
                }
            }
        }
        
        gSystem->ClosePipe(pipe);
        
        // Write summary for this dataset
        summary << "\nDataset: " << eventFiles[fileIdx].second << "\n";
        summary << "Total events: " << totalEvents << "\n";
        summary << "Events with #particles!=8: " << oddEvents << "\n";
        summary << "Instances of WL (pdg=24, state=2): " << WL_instances << "\n";
        summary << "Successful reco: " << successfulReco << " ("
                << (totalEvents ? 100.0 * successfulReco / totalEvents : 0) << "%)\n";
        summary << "Mean M_ll = " << h_lep_invmass_vec[fileIdx]->GetMean() << " GeV\n";
        summary << "Mean M_jj = " << h_jj_mass_vec[fileIdx]->GetMean() << " GeV\n";
        summary << "Mean M_WR = " << h_WR_mass_vec[fileIdx]->GetMean() << " GeV\n";
        summary << "Mean M_N1 = " << h_N1_mass_vec[fileIdx]->GetMean() << " GeV\n";
        summary << "--------------------------------------------------\n";
        
        std::cout << "  Processed " << totalEvents << " events" << std::endl;
    }
    
    std::cout << "\nAnalysis complete! Creating overlaid plots with legends..." << std::endl;
    
    TCanvas c("c", "c", 800, 600);
    gStyle->SetOptStat(0);

    auto saveOverlaidPlot = [&](const char* subfolder, const char* plotname, 
                                const std::vector<TH1F*>& hists, 
                                const std::string& drawOption = "HIST") {
        if (hists.empty()) return;
        
        // Check if this is charge plot
        TString plotNameStr(plotname);
        bool isChargePlot = plotNameStr.Contains("charge");
        bool isLep_pt = plotNameStr.Contains("lep1_pt");
        // Find maximum for scaling
        double maxY = 0;
        if (isChargePlot) {
            // Force maximum for charge plot
            maxY = 5100;
        } else {
            // Auto-scale for other plots
            for (auto hist : hists) {
                double histMax = hist->GetMaximum();
                if (histMax > maxY) maxY = histMax;
            }
            maxY *= 1.3;  // Add margin
        }
        
        // Draw first histogram
        hists[0]->SetMaximum(maxY);
        hists[0]->Draw(drawOption.c_str());
        
        // Draw remaining histograms
        for (size_t i = 1; i < hists.size(); i++) {
            hists[i]->Draw((drawOption + " SAME").c_str());
        }
        

        // Create Textboxes for peaks
        // so that text apears above loval minima need a silly check with a hardcoded x-axis threshold
        // note: hardcoded m_N = 500
        
        TLatex lat;
        lat.SetTextSize(0.025);
        lat.SetTextAlign(13);

        for (size_t i = 0; i < h_lep_pt_vec.size(); i++) {
            if (isLep_pt && eventFiles[i].second.find("m_{N}") != std::string::npos){
                TH1* h = h_lep1_pt_vec[i];
                // int bin = h_lep_pt_vec[i]->GetMaximumBin();
                double xmin_cut  = h->GetXaxis()->GetXmax()/5;
                int binMin = h->GetXaxis()->FindBin(xmin_cut);
                int binMax = h->GetNbinsX();
                int bestBin = -1;
                double bestVal = -1.0;

                for (int b = binMin; b <= binMax; ++b) {
                    double val = h->GetBinContent(b);
                    if (val > bestVal) {
                        bestVal = val;
                        bestBin = b;
                    }
                }
                // Safety check: no valid maximum found
                if (bestBin < 0 || bestVal <= 0)
                    continue;

                double x = h->GetBinCenter(bestBin);
                double y = bestVal;

                lat.SetTextColor(h_lep_pt_vec[i]->GetLineColor());
                lat.DrawLatex(
                    x + 60,
                    y + 28,
                    Form("m_{W_{R}}=%zuTeV", (i % 4) + 3)
                );
            }
        }

        // Create and draw legend
        TLegend* leg = new TLegend(0.55, 0.65, 0.85, 0.9);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        
        for (size_t i = 0; i < hists.size(); i++) {
            leg->AddEntry(hists[i], eventFiles[i].second.c_str(), "l");
            leg->SetTextFont(32);
        }
        
        leg->Draw();
        
        // Save with proper TString concatenation
        TString fullPath = TString::Format("%sKS_Figures_FullLRSM2/%s/%s", baseFigPath.Data(), subfolder, plotname);
        c.SaveAs(fullPath);
        c.Clear();
        delete leg;
    };


    // Lepton plots
    saveOverlaidPlot("leptons", "lep_invmass_overlaid.png", h_lep_invmass_vec);
    saveOverlaidPlot("leptons", "lep1_pt_overlaid.png", h_lep1_pt_vec);
    saveOverlaidPlot("leptons", "lep2_pt_overlaid.png", h_lep2_pt_vec);
    saveOverlaidPlot("leptons", "lep_pt_overlaid.png", h_lep_pt_vec);
    saveOverlaidPlot("leptons", "lep_deltaR_overlaid.png", h_lep_deltaR_vec);
    saveOverlaidPlot("leptons", "lep_charge_overlaid.png", h_ll_charge_vec, "HIST");
    saveOverlaidPlot("leptons", "dilepton_flavor_overlaid.png", h_ll_flavor_vec, "HIST");

    // Jet plots
    saveOverlaidPlot("jets", "jet_pt_overlaid.png", h_jet_pt_vec);
    saveOverlaidPlot("jets", "jj_mass_overlaid.png", h_jj_mass_vec);
    saveOverlaidPlot("jets", "jj_deltaR_overlaid.png", h_jj_deltaR_vec);

    // Reconstruction plots
    saveOverlaidPlot("reconstruction", "WR_mass_overlaid.png", h_WR_mass_vec);
    saveOverlaidPlot("reconstruction", "N1_mass_overlaid.png", h_N1_mass_vec);
    saveOverlaidPlot("reconstruction", "l2jj_deltaR_overlaid.png", h_l2jj_deltaR_vec);
    
    summary << "\n=== ANALYSIS COMPLETE ===" << std::endl;
    summary << "Overlaid plots created with 3 datasets shown together with legend." << std::endl;
    summary.close();
    
    std::cout << "\nAll overlaid plots saved to: " << baseFigPath << "KS_Figures_FullLRSM2/" << std::endl;
    std::cout << "Summary saved to: " << baseFigPath << "KS_Figures_FullLRSM2/summary.txt" << std::endl;
}