#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <optional>
#include <random>
#include <iomanip>

#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"

struct Particle {
    int pdgId;
    int status;
    int mother;
    double px, py, pz, E;
};

bool sortByPt(const TLorentzVector& a, const TLorentzVector& b) {
    return a.Pt() > b.Pt();
}

void KS_ApplyKinematicCuts() {
    std::cout << "User Test 004:" << std::endl;
    std::cout << "Applying Kinematic Cuts to (lljj)" << std::endl;
    gStyle->SetOptStat(0);

    const char* absFigPath = "/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/KS_KinematicCuts/Figures/";
    TString baseFigPath = absFigPath;
    gSystem->mkdir(baseFigPath, kTRUE);

    std::ofstream summary("/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/KS_KinematicCuts/CutEfficiencySummary.txt");
    summary << std::fixed << std::setprecision(2);

    // ================================================================================================================
    // Defining Kinematic Cuts (comment out = no cut)
    // ================================================================================================================
    // use: [ = std::nullopt; ] to pass all, note: efficiecy log = 100%
    std::optional<double> pt_cut_l1 = 100.0;   // leading lepton
    std::optional<double> pt_cut_l2 = 20.0;    // subleading lepton
    std::optional<double> jet_pt_cut = 30.0;   // pT > pT_cut to pass
    std::optional<double> eta_cut = 2.5;       // |η| < eta_cut to pass
       
    // ================================================================================================================
    // Defining Event Files
    // ================================================================================================================
    std::vector<std::pair<std::string,std::string>> eventFiles = {
        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR3000/Events/run_01/unweighted_events.lhe.gz","m_{N}=500GeV, m_{W_{R}}=3TeV"},
        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR4000/Events/run_01/unweighted_events.lhe.gz","m_{N}=500GeV, m_{W_{R}}=4TeV"},
        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR5000/Events/run_01/unweighted_events.lhe.gz","m_{N}=500GeV, m_{W_{R}}=5TeV"},
        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR6000/Events/run_01/unweighted_events.lhe.gz","m_{N}=500GeV, m_{W_{R}}=6TeV"},
        {"/storage/alexanum/WORKBENCH/runs/KS_offWR_MN500_MWR7000/Events/run_01/unweighted_events.lhe.gz","m_{N}=500GeV, m_{W_{R}}=7TeV"},
    };


    // Create Histograms per dataset
    std::vector<TH1F*> h_l1_pt, h_l2_pt, h_N_mass, h_WR_mass;
    std::vector<int> colors = {kRed, kBlue, kGreen, kBlack, kViolet};

    for (size_t i=0;i<eventFiles.size();++i) {
        h_l1_pt.push_back(new TH1F(Form("h_l1_pt_%zu",i),";p_{T}(l_{1}) [GeV];Events",100,0,4000));
        h_l2_pt.push_back(new TH1F(Form("h_l2_pt_%zu",i),";p_{T}(l_{2}) [GeV];Events",100,0,2000));
        h_N_mass.push_back(new TH1F(Form("h_N_%zu",i),";m_{N} [GeV];Events",100,0,3000));
        h_WR_mass.push_back(new TH1F(Form("h_WR_%zu",i),";m_{W_{R}} [GeV];Events",100,0,8000));

        h_l1_pt[i]->SetLineColor(colors[i]);
        h_l2_pt[i]->SetLineColor(colors[i]);
        h_N_mass[i]->SetLineColor(colors[i]);
        h_WR_mass[i]->SetLineColor(colors[i]);
    }

    // Loop over each Event File
    for (size_t fileIdx=0; fileIdx<eventFiles.size(); ++fileIdx) {
        std::cout << "\nProcessing: " << eventFiles[fileIdx].second << std::endl;
        std::cout << "File: " << eventFiles[fileIdx].first << std::endl;
        const char* filename = eventFiles[fileIdx].first.c_str();
        TString cmd = TString::Format("gunzip -c %s", filename);
        FILE* pipe = gSystem->OpenPipe(cmd, "r");
        if (!pipe) {
            std::cerr << "ERROR: Cannot open LHE file: " << filename << std::endl;
            continue;
        }
        // Initiate pool of parant-less jets and leptons
        // jets are not primary/secondary so they come in pair-vector
        std::vector<TLorentzVector> l1_pool, l2_pool;
        // std::vector<std::pair<TLorentzVector,std::vector<TLorentzVector>>> jet_pool; 
        std::vector<std::vector<TLorentzVector>> jet_pool;

        // initiate filter efficiency count
        size_t n_lep_eta_pass = 0;
        size_t n_lep_loPT_pass = 0;
        size_t n_lep_hiPT_pass = 0;
        size_t n_jet_eta_pass = 0;
        size_t n_jet_pt_pass = 0;

        char line[4096];
        bool inEvent=false, firstLine=false;
        std::vector<Particle> particles;

        while (fgets(line,sizeof(line),pipe)) {
            std::string s(line);
            s.erase(0,s.find_first_not_of(" \t\n\r"));
            s.erase(s.find_last_not_of(" \t\n\r") + 1);
            if (s.empty() || s[0] == '#') continue;

            if (s.find("<event>")!=std::string::npos) {
                inEvent=true; 
                firstLine=true; 
                particles.clear(); 
                continue;
            }

            if (s.find("</event>")!=std::string::npos) {
                std::vector<TLorentzVector> leptons;
                std::vector<TLorentzVector> jets;

                for (auto &p: particles) {
                    if (p.status!=1) continue; // not final state particle

                    TLorentzVector v(p.px,p.py,p.pz,p.E);
                    int apid = abs(p.pdgId);

                    if (apid==11 || apid==13 || apid==15) {
                        // Apply Eta cut to leptons
                        // if eta_cut undefined: pass all events
                        if (eta_cut && fabs(v.Eta())>*eta_cut) continue;
                        n_lep_eta_pass++;
                        leptons.push_back(v);
                    }

                    if (apid>=1 && apid<=6) {
                        // Apply Eta and pT cut to jets
                        //if cut undefined: pass all events
                        if (eta_cut && fabs(v.Eta())>*eta_cut) continue;
                        n_jet_eta_pass++;
                        if (jet_pt_cut && v.Pt()<*jet_pt_cut) continue;
                        n_jet_pt_pass++;
                        jets.push_back(v);
                    }
                }

                if (leptons.size()<2 || jets.size()<2) { inEvent=false; continue; }

                std::sort(leptons.begin(),leptons.end(),sortByPt);
                std::sort(jets.begin(),jets.end(),sortByPt);

                TLorentzVector l1 = leptons[0];
                TLorentzVector l2 = leptons[1];

                if (pt_cut_l1 && l1.Pt()<*pt_cut_l1) { inEvent=false; continue; }
                n_lep_hiPT_pass++;
                if (pt_cut_l2 && l2.Pt()<*pt_cut_l2) { inEvent=false; continue; }
                n_lep_loPT_pass++;
                l1_pool.push_back(l1);
                l2_pool.push_back(l2);
                // jet_pool.emplace_back(l2,std::vector<TLorentzVector>{jets[0],jets[1]});
                jet_pool.push_back({jets[0], jets[1]});

                h_l1_pt[fileIdx]->Fill(l1.Pt());
                h_l2_pt[fileIdx]->Fill(l2.Pt());

                inEvent=false;
                continue;
            }

            if (inEvent) {
                if (firstLine) { firstLine=false; continue; }
                if (s[0]=='<') continue;

                Particle p;
                int m2, c1, c2; 
                double m, l, sp;
                std::istringstream iss(s);
                if (iss>>p.pdgId>>p.status>>p.mother>>m2>>c1>>c2
                    >>p.px>>p.py>>p.pz>>p.E>>m>>l>>sp) {
                    particles.push_back(p);
                }
            }
        } // Pipe loop closing bracket
        gSystem->ClosePipe(pipe);


        // Shuffle Pool Vectors and Trim if Uneven Sizes
        // auto trimToMin = [&](auto& a, auto& b) {
        //     if (a.size()>b.size()) {
        //         std::sort(a.begin(),a.end(),
        //             [](auto& x, auto& y){return x.Pt()<y.Pt();});
        //         a.resize(b.size());
        //     }
        // };

        // trimToMin(l1_pool,l2_pool);
        // trimToMin(l2_pool,l1_pool);
        // // Jets are allowed to be recorded together (come in pairs). But the l2 is shuffled
        // std::shuffle(l1_pool.begin(),l1_pool.end(),std::mt19937{std::random_device{}()});
        // std::shuffle(l2_pool.begin(),l2_pool.end(),std::mt19937{std::random_device{}()});
        // std::shuffle(jet_pool.begin(),jet_pool.end(),std::mt19937{std::random_device{}()});

        // size_t nReco = std::min({l1_pool.size(),l2_pool.size(),jet_pool.size()});

        // for (size_t i=0;i<nReco;++i) {
        //     TLorentzVector N = l2_pool[i] + jet_pool[i].second[0] + jet_pool[i].second[1];
        //     TLorentzVector WR = l1_pool[i] + N;
        //     h_N_mass[fileIdx]->Fill(N.M());
        //     h_WR_mass[fileIdx]->Fill(WR.M());
        // }

        /*
        The Great Equilization: Compare zie of pool l1 and l2 and jets
        --------------------------------------------------------------
        From big to small lepton pool:
            l2 majority:  highest l2 promote to l1
            l1 majority: lowest l1 fall to l2
        Compare jet group and leptons:
            jet majority: purge jets (smallest summed pT suffer)
            jet minority purge equally l1 and l2 (smallest pT suffer)
        */

        if (l1_pool.size() > l2_pool.size()) {
            std::sort(l1_pool.begin(), l1_pool.end(),
                    [](auto& a, auto& b){ return a.Pt() < b.Pt(); });

            size_t diff = l1_pool.size() - l2_pool.size();
            for (size_t i = 0; i < diff; ++i) {
                l2_pool.push_back(l1_pool.front());
                l1_pool.erase(l1_pool.begin());
            }
        }

        if (l2_pool.size() > l1_pool.size()) {
            std::sort(l2_pool.begin(), l2_pool.end(),
                    [](auto& a, auto& b){ return a.Pt() > b.Pt(); });

            size_t diff = l2_pool.size() - l1_pool.size();
            for (size_t i = 0; i < diff; ++i) {
                l1_pool.push_back(l2_pool.front());
                l2_pool.erase(l2_pool.begin());
            }
        }

        size_t nMin = std::min({l1_pool.size(), l2_pool.size(), jet_pool.size()});

        // Trim l1
        std::sort(l1_pool.begin(), l1_pool.end(),
                [](auto& a, auto& b){ return a.Pt() < b.Pt(); });
        l1_pool.resize(nMin);

        // Trim l2
        std::sort(l2_pool.begin(), l2_pool.end(),
                [](auto& a, auto& b){ return a.Pt() < b.Pt(); });
        l2_pool.resize(nMin);

        // Trim jets by summed pT
        std::sort(jet_pool.begin(), jet_pool.end(),
                [](auto& a, auto& b){
                    return (a[0].Pt()+a[1].Pt()) < (b[0].Pt()+b[1].Pt());
                });
        jet_pool.resize(nMin);

        // Shuffle
        std::mt19937 rng{std::random_device{}()};
        std::shuffle(l1_pool.begin(), l1_pool.end(), rng);
        std::shuffle(l2_pool.begin(), l2_pool.end(), rng);  
        std::shuffle(jet_pool.begin(), jet_pool.end(), rng);

        // Reconstruction
        for (size_t i=0; i<nMin; ++i) {
            TLorentzVector N  = l2_pool[i] + jet_pool[i][0] + jet_pool[i][1];
            TLorentzVector WR = l1_pool[i] + N;

            h_N_mass[fileIdx]->Fill(N.M());
            h_WR_mass[fileIdx]->Fill(WR.M());
        }

        // Log Cut Efficiencies per File
        constexpr double N_EVENTS = 10000.0;
        summary << "=============================================\n";
        summary << "Dataset: " << eventFiles[fileIdx].second << "\n";
        summary << "File:    " << eventFiles[fileIdx].first << "\n\n";

        summary << "Lepton efficiencies:\n";
        summary << "  Passed eta cut        : " << 100.0 * n_lep_eta_pass / N_EVENTS / 2 << "\n";
        summary << "  Passed l1 pT cut      : " << 100.0 * n_lep_hiPT_pass / N_EVENTS << "\n";
        summary << "  Passed l2 pT cut      : " << 100.0 * n_lep_loPT_pass / N_EVENTS << "\n\n";

        summary << "Jet efficiencies:\n";
        summary << "  Passed eta cut        : " << 100.0 * n_jet_eta_pass / N_EVENTS / 2 << "\n";
        summary << "  Passed pT cut         : " << 100.0 * n_jet_pt_pass / N_EVENTS / 2<< "\n\n";

        summary << "Final pool sizes after equalization:\n";
        summary << "  l1 pool size          : " << l1_pool.size() << "\n";
        summary << "  l2 pool size          : " << l2_pool.size() << "\n";
        summary << "  jet-pair pool size    : " << jet_pool.size() << "\n\n";

    } // close of File loop

    // Plot Histograms
    TCanvas c("c","c",800,600);
    TLegend leg(0.55,0.65,0.85,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.SetTextFont(32);

    auto draw = [&](std::vector<TH1F*>& v,const char* name){
        double max=0;
        for(auto h:v) max=std::max(max,h->GetMaximum());
        for(size_t i=0;i<v.size();++i){
            v[i]->SetMaximum(1.3*max);
            v[i]->Draw(i?"HIST SAME":"HIST");
            leg.AddEntry(v[i],eventFiles[i].second.c_str(),"l");
        }
        leg.Draw();
        TString fullPath = TString::Format("%s%s", baseFigPath.Data(), name);
        c.SaveAs(fullPath);
        c.Clear(); 
        leg.Clear();
    };

    draw(h_l1_pt,"l1_pt.png");
    draw(h_l2_pt,"l2_pt.png");
    draw(h_N_mass,"N_mass.png");
    draw(h_WR_mass,"WR_mass.png");

    std::cout << "\nAll plots saved to: " << baseFigPath << std::endl;
    std::cout << "Summary saved to: " << ".../KS_KinematicCuts/CutEfficiencySummary.txt" << std::endl;
}