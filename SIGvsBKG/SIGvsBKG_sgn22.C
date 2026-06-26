#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"
#include "TLorentzVector.h"
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

/*
(==================================-RUN ME-==================================)
Do setup in in direcoty above /ROOT_files, /scripts
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase  
    alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh' 
    setupATLAS -c centos7 
    asetup Athena,24.0.11 
Compile/Run:
    root -l -q -b SIG_vs_BKG2.C
*/

struct Sample {
    std::string name;
    std::string type;  // f.ex: "SIG", "BKG_Z", "BKG_T"
    double sigma_pb;
    double lumi_pb;
    double filter_eff;
    double k_factor;
    std::vector<std::string> files;
};

// Global selection cuts
const float lep_pt_min  = 10.0;  // GeV
const float lep_eta_max = 2.5;
const float jet_pt_min  = 10.0;  // GeV
const float jet_eta_max = 2.5;
const float lep_rel_ptvarcone20_max  = 0.15;
const float lep_rel_topoetcone20_max = 0.15;
const float jvt_min = 0.5;  
const float jWidth_cut = 0.05;
const float jNumTrkPt500_cut = 10;

struct MyLepton {
    float pt;
    float eta;
    float phi;
    float charge;
    char type;  // 'e' or 'm'
};

struct MyJet {
    float pt;
    float eta;
    float phi;
    float m;
    float Jvt;
    float Width;
    float NumTrkPt500;
};

struct MyLepTrack {
    float pt;
    float x;
    float y;
    float z;
    char type;
};

// sort by pT (descending)
bool sortByPtLep(const MyLepton &a, const MyLepton &b) { 
    return a.pt > b.pt;
}

bool sortByPtJet(const MyJet &a, const MyJet &b) { 
    return a.pt > b.pt; 
}

std::string basename(const std::string &path) {
    auto pos = path.find_last_of("/\\");
    return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

std::string extract_meta_topkey(const std::string &path) {
    // Extract the top-level JSON key from the file path.
    // BKG: ...mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu...               -> "601190" (DSID)
    // SIG: ...mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_70_2L2J_MadSpin... -> "MGPy8EG_LRSM_WR_HNL_5000_70_2L2J_MadSpin" (generator tag)
    const std::string tag = "TeV.";
    auto pos = path.find(tag);
    if (pos == std::string::npos) return "";
    pos += tag.size();
    auto pos2 = path.find('.', pos);
    if (pos2 == std::string::npos) return "";
    return path.substr(pos, pos2 - pos); 
}

std::string extract_meta_filekey(const std::string &path, const json &topBlock) {
    // Decide which file key to use inside meta[topkey]
    std::string base = basename(path);
    if (topBlock.contains(base)) {
        return base;
    }
    if (topBlock.is_object() && topBlock.size() == 1) {
        return topBlock.begin().key();
    }
    return base;
}

double get_total_sumw_from_metadata(const std::vector<std::string> &file_paths, const json &meta) {
    /*
    Sum AllExecutedEvents.W over all files in this sample.
    Inputs:
        topkey = token after 'TeV.' in the filename (the BKG DSID or SIG generator tag)
        filekey = basename or (if only one key) the single inner key
    */
    double totalW = 0.0;
    for (const auto &fullPath : file_paths) {
        std::string topkey = extract_meta_topkey(fullPath);
        if (topkey.empty()) {
            std::cerr << "WARNING: could not extract metadata topkey from path: " << fullPath << "\n";
            continue;
        }
        auto itTop = meta.find(topkey);
        if (itTop == meta.end()) {
            std::cerr << "WARNING: topkey " << topkey << " not found in metadata\n";
            continue;
        }
        const json &topBlock = *itTop;
        std::string fileKey = extract_meta_filekey(fullPath, topBlock);
        auto itFile = topBlock.find(fileKey);
        if (itFile == topBlock.end()) {
            std::cerr << "WARNING: file key " << fileKey << " not found in metadata for topkey " << topkey << "\n";
            continue;
        }
        const json &fileBlock = *itFile;
        auto itAll = fileBlock.find("AllExecutedEvents");
        if (itAll == fileBlock.end()) {
            std::cerr << "WARNING: AllExecutedEvents not found for " << fileKey << " (topkey " << topkey << ")\n";
            continue;
        }
        const json &allExec = *itAll;
        if (!allExec.contains("W")) {
            std::cerr << "WARNING: AllExecutedEvents has no W for " << fileKey << " (topkey " << topkey << ")\n";
            continue;
        }
        double W = allExec["W"].get<double>();
        totalW += W;
    }
    return totalW;
}



// DATA Quality 
std::map< TString, std::vector<UInt_t> > GRL;
std::vector<UInt_t > GRL_vec;
std::map< UInt_t, std::vector<UInt_t> > GRL_lb;

bool isCleanEvent(UInt_t larFlag, UInt_t tileFlag, UInt_t sctFlag, UInt_t coreFlag) {
    // Called like this for cleaning (every event)
    //isCleanEvent(EventInfoAuxDyn.larFlags, EventInfoAuxDyn.tileFlags, EventInfoAuxDyn.sctFlags, EventInfoAuxDyn.coreFlags)
    bool isCleanEvent = !( ((( larFlag >> 28 ) & 0xF) == 2) || ((( tileFlag >> 28 ) & 0xF) == 2) || ((( sctFlag >> 28 ) & 0xF) == 2) || ( ( coreFlag >> 18 ) & 0x1 ));
    return isCleanEvent;
}


void setGRL(int vb) {
    // Called in initialize (before starting loop over events)
    std::ifstream input("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data22_13p6TeV/20250321/data22_13p6TeV.periodAllYear_DetStatus-v134-pro28-09_MERGED_PHYS_StandardGRL_All_Good_25ns.xml");
    bool startLB = false;
    bool firstRun = true;
    std::vector<UInt_t> thislbvec;
    UInt_t thisrun;
    for( std::string line; getline( input, line );) {
        TString l = line;
        // <Run>430536</Run>
        if(l.Contains("<Run>") && l.Contains("</Run>")){
            if(!firstRun) {
                //GRL_lb[thisrun] = thislbvec;
                GRL_lb.insert( std::make_pair(thisrun, thislbvec) );
                thislbvec.clear();
                //firstRun = false;
            }
            int st = l.First('>');
            int en = 6;//l.Last("Run");
            if(vb)printf("Index start = %i, stop = %i\n",st,en);
            TString runnum( l(st+1,en) );
            thisrun = (UInt_t)runnum.Atoi();
            firstRun = false;
        } else if(l.Contains("<LBRange")) {
            TObjArray *tx = l.Tokenize("\"");
            UInt_t start = (UInt_t)(((TObjString *)tx->At(1))->String()).Atoi();
            UInt_t end   = (UInt_t)(((TObjString *)tx->At(3))->String()).Atoi();
            if(vb)printf("Run %i has LB start = %i and end = %i\n",thisrun,start,end);
            for(UInt_t i = start; i<=end; i++) {
                //printf("Adding %i to run %i\n",i,thisrun);
                thislbvec.push_back(i);
            } 
        }
    }
  // Need to get the last run too!
  GRL_lb.insert( std::make_pair(thisrun, thislbvec) );
  thislbvec.clear();
}


Int_t checkLB(UInt_t rnum, UInt_t lb, int vb) {
    // Called for every event like this:
    // checkLB(EventInfoAuxDyn.runNumber,EventInfoAuxDyn.lumiBlock)
    if(vb) {
        printf("Checking LB %i\n",lb);
        for (auto x = GRL_lb.begin(); x != GRL_lb.end(); x++) {
            if(x->first != rnum)continue;
            printf("Run %i has following LBs:\n",x->first);
            for(UInt_t j= 0; j< (x->second).size(); j++) {
                printf("%i  ",(x->second).at(j));
            }
        }
    }
    int isGRL = 0;//GRLlist.HasRunLumiBlock(rnum,lb);
    if(GRL_lb.find(rnum) != GRL_lb.end()) {
        if ( std::find(GRL_lb[rnum].begin(), GRL_lb[rnum].end(), lb) != GRL_lb[rnum].end() ) {
            isGRL = 1;
        }
    }
    return isGRL;
}

// build leptons from electron + muon collections
void build_leptons(const std::vector<float> &e_pt,
                   const std::vector<float> &e_eta,
                   const std::vector<float> &e_phi,
                   const std::vector<float> &e_charge,    
                   const std::vector<float> &e_ptvarcone20,
                   const std::vector<float> &e_topoetcone20,
                   const std::vector<char>  &e_DFquality,
                   const std::vector<float> &m_pt,
                   const std::vector<float> &m_eta,
                   const std::vector<float> &m_phi,
                   const std::vector<float> &m_charge,   
                   const std::vector<float> &m_ptvarcone20,
                   const std::vector<float> &m_topoetcone20,
                   std::vector<MyLepton> &out_Leptons)
{
    out_Leptons.clear();

    // electrons
    for (size_t i = 0; i < e_pt.size(); ++i) {
        float pt_mev = e_pt[i];
        if (pt_mev <= 0) continue;
        float pt  = pt_mev * 0.001f; // MeV -> GeV
        float eta = e_eta[i];
        float phi = e_phi[i];
        if (pt <= lep_pt_min) continue;
        if (std::fabs(eta) > lep_eta_max) continue;
        if (!e_DFquality[i]) continue;
        float rel_ptvar = e_ptvarcone20[i] / pt_mev;
        float rel_topo = e_topoetcone20[i] / pt_mev;
        if (rel_ptvar > lep_rel_ptvarcone20_max) continue;
        if (rel_topo > lep_rel_topoetcone20_max) continue;

        MyLepton lep;
        lep.pt = pt;
        lep.eta = eta;
        lep.phi = phi;
        lep.charge = e_charge[i];
        lep.type = 'e';
        out_Leptons.push_back(lep);
    }

    // muons
    for (size_t i = 0; i < m_pt.size(); ++i) {
        float pt_mev = m_pt[i];
        if (pt_mev <= 0) continue;
        float pt = pt_mev * 0.001f; // MeV -> GeV
        float eta = m_eta[i];
        float phi = m_phi[i];
        if (pt <= lep_pt_min) continue;
        if (std::fabs(eta) > lep_eta_max) continue;
        float rel_ptvar = m_ptvarcone20[i] / pt_mev;
        float rel_topo = m_topoetcone20[i] / pt_mev;
        if (rel_ptvar > lep_rel_ptvarcone20_max) continue;
        if (rel_topo > lep_rel_topoetcone20_max) continue;

        MyLepton lep;
        lep.pt = pt;
        lep.eta = eta;
        lep.phi = phi;
        lep.charge = m_charge[i];
        lep.type = 'm';
        out_Leptons.push_back(lep);
    }
    std::sort(out_Leptons.begin(), out_Leptons.end(), sortByPtLep);
}

// build jets
void build_jets(const std::vector<float> &j_pt,
                const std::vector<float> &j_eta,
                const std::vector<float> &j_phi,
                const std::vector<float> &j_m,
                const std::vector<float> &j_Jvt,
                const std::vector<float> &j_Width,
                const std::vector<std::vector<int>> &j_NumTrkPt500,
                std::vector<MyJet> &out_Jets)
{
    out_Jets.clear();
    size_t n = j_pt.size();
    out_Jets.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        if (j_pt[i] * 0.001f <= jet_pt_min) continue;
        if (std::fabs(j_eta[i]) > jet_eta_max) continue;
        if (j_Jvt[i] < jvt_min) continue;
        if (j_Width[i] > jWidth_cut) continue;
        int nTrk = 0;
        if (i < j_NumTrkPt500.size() && !j_NumTrkPt500[i].empty()) {
            nTrk = j_NumTrkPt500[i][0];  // 0th entry for PV
        }
        // if (nTrk > jNumTrkPt500_cut) continue;
        MyJet jet;
        jet.pt = j_pt[i] * 0.001f;
        jet.eta = j_eta[i];
        jet.phi = j_phi[i];
        jet.m = j_m[i] * 0.001f;
        jet.Jvt = j_Jvt[i];
        jet.Width = j_Width[i];
        jet.NumTrkPt500 = nTrk;
        out_Jets.push_back(jet);
    }
    std::sort(out_Jets.begin(), out_Jets.end(), sortByPtJet);
}

void fill_reco_variables(const std::string &sname,
                         const std::vector<MyLepton> &leptons,
                         const std::vector<MyJet> &jets_sel,
                         double w_event,
                         bool buff_weight,
                         std::map<std::string, std::map<std::string, TH1F*>> &histos)
{
    if (jets_sel.size() < 2 || leptons.size() < 2) return;
    
    // Build jet system jj
    TLorentzVector j1, j2, jj;
    j1.SetPtEtaPhiM(jets_sel[0].pt, jets_sel[0].eta, jets_sel[0].phi, jets_sel[0].m);
    j2.SetPtEtaPhiM(jets_sel[1].pt, jets_sel[1].eta, jets_sel[1].phi, jets_sel[1].m);
    jj = j1 + j2;
    double dR_jj = j1.DeltaR(j2);
    // Validate dijet system. If wrong any reconstruction is futile in this event
    // In this version is categorical. We do not make a sliding DR(j[i]j[i+1]) window. 
    if (dR_jj > 2.0) {
        return;
    }
    // Leptons arrive sorted by pT: leptons[0] = l1, leptons[1] = l2
    // Select same sign signature
    float qprod = leptons[0].charge*leptons[1].charge;   //double
    if (qprod<0) return;

    double pt1 = leptons[0].pt;
    double pt2 = leptons[1].pt;
    if (pt1 <= 0 || pt2 <= 0) return; // Safety: avoid log of non‑positive pT

    // Decision logic
    bool use_event = false;
    int  idx_lep_for_m = -1;  // index of lepton to be used in m(l2jj)

    // 1) pT-based criterion
    double logdiff = std::log(pt1) - std::log(pt2); 
    if (logdiff > 1.5) {
        // leptons are discriminated by pT, keep input ordering and use [1] for mN
        idx_lep_for_m = 1;  
        use_event = true;
    } else {
        // 2) DR-based criterion
        TLorentzVector l1, l2;
        l1.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, 0.0);
        l2.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, 0.0);
        double dR1 = l1.DeltaR(jj);
        double dR2 = l2.DeltaR(jj);
        // If BOTH DR(l1,jj) and DR(l2,jj) are between 1 and 2:
        // if (dR1 > 0.0 && dR1 < 2.0 &&dR2 > 0.0 && dR2 < 2.0) {
            if (fabs(dR1 - dR2) >= 1.0) {
            // assumption "l2 has lowest DR" is valid 
            idx_lep_for_m = (dR1 < dR2 ? 0 : 1);
            use_event = true;
        } else {
            // 3) Neither pT nor DR criteria satisfied -> kill this event for m(l2jj)
            return;
        }
    }

    // Get histogram map once
    auto itSample = histos.find(sname);
    if (itSample == histos.end()) return;
    auto &varMap = itSample->second;

    // Fill: m(lljj) reco-logic invariant
    TLorentzVector l1_vec, l2_vec;
    l1_vec.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, 0.0);
    l2_vec.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, 0.0);
    TLorentzVector lljj = l1_vec + l2_vec + j1 + j2;
    double m_lljj = lljj.M();
    auto itHist_lljj = varMap.find("m_lljj");
    if (itHist_lljj != varMap.end()) {
        if (buff_weight) {
            itHist_lljj->second->Fill(m_lljj, w_event*5e4);  //1e2  5e1   //c
        } else {
            itHist_lljj->second->Fill(m_lljj, w_event);    //w_event*2e5
        }
    }
    
    if (!use_event || idx_lep_for_m < 0) return;

    // Fill: m(l2jj) reco-logic dependent
    TLorentzVector l2jj_lep;
    l2jj_lep.SetPtEtaPhiM(leptons[idx_lep_for_m].pt, leptons[idx_lep_for_m].eta, leptons[idx_lep_for_m].phi, 0.0);
    TLorentzVector l2jj = jj + l2jj_lep;
    double m_l2jj = l2jj.M();
    auto itHist_l2jj = varMap.find("m_l2jj");
    if (itHist_l2jj != varMap.end()) {
        if (buff_weight) {
            itHist_l2jj->second->Fill(m_l2jj, w_event*5e4);   //1e3  5e3
        } else {
            itHist_l2jj->second->Fill(m_l2jj, w_event);   //w_event*2e5
        }
    }

    // Fill Again Selected LLJJ if l1 and l2 are distinguished
    TLorentzVector l1_vec_sel, l2_vec_sel;
    l1_vec_sel.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, 0.0);
    l2_vec_sel.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, 0.0);
    TLorentzVector lljj_sel = l1_vec + l2_vec + j1 + j2;
    double m_lljj_sel = lljj_sel.M();
    auto itHist_lljj_sel = varMap.find("m_lljj_sel");
    if (itHist_lljj_sel != varMap.end()) {
        if (buff_weight) {
            itHist_lljj_sel->second->Fill(m_lljj_sel, w_event*5e4);  //1e2  5e1   //c
        } 
    }
}


void spike_killer(TH1F *h, double min_y, double r_factor=1.5) {
    // Depricated Function, used to remove isolated histogram bins (likely background activity in the BKG files)
    int nb = h->GetNbinsX();
    for (int ib = 2; ib < nb; ++ib) { // skip first & last
        double y  = h->GetBinContent(ib);
        double ym = h->GetBinContent(ib - 1);
        double yp = h->GetBinContent(ib + 1);
        double y_next = 0.5 * (ym + yp);
        if (ym < min_y/r_factor || yp < min_y/r_factor) {
            h->SetBinContent(ib, 0);
            h->SetBinError(ib, 0);
        }
    }
}

std::vector<double> make_log_bins(double xmin, double xmax, int nbins) {
    // Build progressively wider bins between xmin and xmax.
    std::vector<double> edges;
    edges.reserve(nbins + 1);
    // Require positive xmin for log.
    if (xmin <= 0) xmin = 1e-3;
    double log_min = std::log10(xmin);
    double log_max = std::log10(xmax);
    double dlog    = (log_max - log_min) / nbins;
    for (int i = 0; i <= nbins; ++i) {
        double log_x = log_min + i * dlog;
        edges.push_back(std::pow(10.0, log_x));
    }
    return edges;
}

// Helper to build a TChain from sample.files
TChain* build_chain(const Sample &s, const std::string &treename) {
    TChain *chain = new TChain(treename.c_str());
    for (const auto &fname : s.files) {
        chain->Add(fname.c_str());
    }
    return chain;
}

// Process one sample: fill histos[sname][var] and compute super_factor
void process_sample(const Sample &sample,
                    std::map<std::string, std::map<std::string, TH1F*>> &histos,
                    std::map<std::string, double> &super_factors,
                    const json *metadata = nullptr,
                    bool use_metadata_totalW = true)
{
    const std::string &sname = sample.name;
    const std::string &stype = sample.type;
    std::cout << "Processing sample: " << sname << " (" << stype << ")\n";
    TChain *chain = build_chain(sample, "CollectionTree");
    // get the list of appended files in the chain for metadata
    std::vector<std::string> chain_files;
    {
        TObjArray *fileElements = chain->GetListOfFiles();
        for (int i = 0; i < fileElements->GetEntries(); ++i) {
            auto *el = dynamic_cast<TChainElement*>(fileElements->At(i));
            if (!el) continue;
            chain_files.emplace_back(el->GetTitle()); // full path of file
        }
    }
    // Slim branches
    chain->SetBranchStatus("*", 0);

    chain->SetBranchStatus("ElectronsAuxDyn.pt",                       1);
    chain->SetBranchStatus("ElectronsAuxDyn.eta",                      1);
    chain->SetBranchStatus("ElectronsAuxDyn.phi",                      1);
    chain->SetBranchStatus("ElectronsAuxDyn.charge",                   1);
    chain->SetBranchStatus("ElectronsAuxDyn.ptvarcone20",              1);
    chain->SetBranchStatus("ElectronsAuxDyn.topoetcone20",             1);
    chain->SetBranchStatus("ElectronsAuxDyn.DFCommonElectronsLHLoose", 1);  //Loose

    chain->SetBranchStatus("MuonsAuxDyn.pt",           1);
    chain->SetBranchStatus("MuonsAuxDyn.eta",          1);
    chain->SetBranchStatus("MuonsAuxDyn.phi",          1);
    chain->SetBranchStatus("MuonsAuxDyn.charge",       1);
    chain->SetBranchStatus("MuonsAuxDyn.ptvarcone20",  1);
    chain->SetBranchStatus("MuonsAuxDyn.topoetcone20", 1);

    chain->SetBranchStatus("AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_pt",  1);
    chain->SetBranchStatus("AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_eta", 1);
    chain->SetBranchStatus("AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_phi", 1);
    chain->SetBranchStatus("AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_m",   1);
    chain->SetBranchStatus("AntiKt4EMPFlowJetsAuxDyn.Jvt",                         1);
    chain->SetBranchStatus("AntiKt4EMPFlowJetsAuxDyn.Width",                       1);
    chain->SetBranchStatus("AntiKt4EMPFlowJetsAuxDyn.NumTrkPt500",                 1);

    chain->SetBranchStatus("PrimaryVerticesAuxDyn.vertexType", 1);
    chain->SetBranchStatus("EventInfoAuxDyn.DFCommonJets_eventClean_LooseBad",     1);
    chain->SetBranchStatus("EventInfoAuxDyn.mcEventWeights", 1);

    // Enable TTree cache
    chain->SetCacheSize(50 * 1024 * 1024); // 50 MB
    chain->AddBranchToCache("*", kTRUE);
    Long64_t nentries = chain->GetEntries();
    std::cout << "  Entries: " << nentries << std::endl;

    // Branch pointers
    std::vector<float> *e_pt           = nullptr;
    std::vector<float> *e_eta          = nullptr;
    std::vector<float> *e_phi          = nullptr;
    std::vector<float> *e_charge       = nullptr;   //float
    std::vector<float> *e_ptvarcone20  = nullptr;
    std::vector<float> *e_topoetcone20 = nullptr;
    std::vector<char> *e_DFquality     = nullptr;

    std::vector<float> *m_pt           = nullptr;
    std::vector<float> *m_eta          = nullptr;
    std::vector<float> *m_phi          = nullptr;
    std::vector<float> *m_charge       = nullptr;   //float
    std::vector<float> *m_ptvarcone20  = nullptr;
    std::vector<float> *m_topoetcone20 = nullptr;

    std::vector<float> *j_pt           = nullptr;
    std::vector<float> *j_eta          = nullptr;
    std::vector<float> *j_phi          = nullptr;
    std::vector<float> *j_m            = nullptr;
    std::vector<float> *j_Jvt          = nullptr;
    std::vector<float> *j_Width        = nullptr;
    std::vector<std::vector<int>> *j_NumTrkPt500 = nullptr;

    std::vector<short> *PV_types = nullptr;
    int eventClean = 0;
    std::vector<float> *mcEventWeights = nullptr;
    
    chain->SetBranchAddress("ElectronsAuxDyn.pt",                       &e_pt);
    chain->SetBranchAddress("ElectronsAuxDyn.eta",                      &e_eta);
    chain->SetBranchAddress("ElectronsAuxDyn.phi",                      &e_phi);
    chain->SetBranchAddress("ElectronsAuxDyn.charge",                   &e_charge);
    chain->SetBranchAddress("ElectronsAuxDyn.ptvarcone20",              &e_ptvarcone20);
    chain->SetBranchAddress("ElectronsAuxDyn.topoetcone20",             &e_topoetcone20);
    chain->SetBranchAddress("ElectronsAuxDyn.DFCommonElectronsLHLoose", &e_DFquality);   //Loose

    chain->SetBranchAddress("MuonsAuxDyn.pt",           &m_pt);
    chain->SetBranchAddress("MuonsAuxDyn.eta",          &m_eta);
    chain->SetBranchAddress("MuonsAuxDyn.phi",          &m_phi);
    chain->SetBranchAddress("MuonsAuxDyn.charge",       &m_charge);
    chain->SetBranchAddress("MuonsAuxDyn.ptvarcone20",  &m_ptvarcone20);
    chain->SetBranchAddress("MuonsAuxDyn.topoetcone20", &m_topoetcone20);

    chain->SetBranchAddress("AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_pt",  &j_pt);
    chain->SetBranchAddress("AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_eta", &j_eta);
    chain->SetBranchAddress("AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_phi", &j_phi);
    chain->SetBranchAddress("AntiKt4EMPFlowJetsAuxDyn.JetConstitScaleMomentum_m",   &j_m);
    chain->SetBranchAddress("AntiKt4EMPFlowJetsAuxDyn.Jvt",                         &j_Jvt);
    chain->SetBranchAddress("AntiKt4EMPFlowJetsAuxDyn.Width",                       &j_Width);
    chain->SetBranchAddress("AntiKt4EMPFlowJetsAuxDyn.NumTrkPt500",                 &j_NumTrkPt500);

    chain->SetBranchAddress("PrimaryVerticesAuxDyn.vertexType", &PV_types);
    chain->SetBranchAddress("EventInfoAuxDyn.DFCommonJets_eventClean_LooseBad", &eventClean);
    chain->SetBranchAddress("EventInfoAuxDyn.mcEventWeights", &mcEventWeights);

    UInt_t larFlag   = 0;
    UInt_t tileFlag  = 0;
    UInt_t sctFlag   = 0;
    UInt_t coreFlag  = 0;
    UInt_t runNumber = 0;
    UInt_t lumiBlock = 0;
    if (stype == "DATA") {
        chain->SetBranchStatus("EventInfoAuxDyn.larFlags",  1);
        chain->SetBranchStatus("EventInfoAuxDyn.tileFlags", 1);
        chain->SetBranchStatus("EventInfoAuxDyn.sctFlags",  1);
        chain->SetBranchStatus("EventInfoAuxDyn.coreFlags", 1);
        chain->SetBranchStatus("EventInfoAuxDyn.runNumber", 1);
        chain->SetBranchStatus("EventInfoAuxDyn.lumiBlock", 1);
        chain->SetBranchAddress("EventInfoAuxDyn.larFlags",  &larFlag);
        chain->SetBranchAddress("EventInfoAuxDyn.tileFlags", &tileFlag);
        chain->SetBranchAddress("EventInfoAuxDyn.sctFlags",  &sctFlag);
        chain->SetBranchAddress("EventInfoAuxDyn.coreFlags", &coreFlag);
        chain->SetBranchAddress("EventInfoAuxDyn.runNumber", &runNumber);
        chain->SetBranchAddress("EventInfoAuxDyn.lumiBlock", &lumiBlock);
    }

    // Single loop: accumulate total_sumw and fill histos with w_nom
    //-----------------------------------------------------------------------------------------------
    double total_sumw_loop = 0.0;
    std::vector<MyLepton> leptons;
    std::vector<MyJet>    jets;
    leptons.reserve(10);
    jets.reserve(10);
    int PVrejectedevents = 0;
    int EventUncleanevents = 0;
    int RejectedData = 0;
    int vb = 0;
    if (stype=="DATA") {
        setGRL(vb); 
    }
    for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
        chain->GetEntry(ievt);
        if (ievt % 10000000 == 0 && ievt !=0) {
            std::cout << "  Yes, I'm aline, @ " << ievt << std::endl;
        }

        if (stype=="DATA") {
            bool datacheckclean = isCleanEvent(larFlag, tileFlag, sctFlag, coreFlag);
            Int_t datacheckLB = checkLB(runNumber, lumiBlock, vb);
            if (!datacheckclean || datacheckLB == 0) {
                RejectedData++;
                continue;
            }
        }

        // Reject events with false PV or eventClean
        if (eventClean != 1) { 
            EventUncleanevents++;
            continue; 
        }

        bool hasGoodPV = false;
        if (PV_types) {
            for (short vtype : *PV_types) {
                if (vtype == 1) {        
                    hasGoodPV = true;
                    break;
                }
            }
        }
        if (!hasGoodPV) {
            PVrejectedevents++;
            continue;
        }   
        double w_nom = 1.0;
        if (stype != "DATA") {
            if (!mcEventWeights || mcEventWeights->empty()) continue;
            w_nom = (*mcEventWeights)[0];
        }   
        total_sumw_loop += w_nom;
        double w_event = w_nom;  // scale later by super_factor
        build_leptons(*e_pt, *e_eta, *e_phi, *e_charge, *e_ptvarcone20, *e_topoetcone20, *e_DFquality,
                      *m_pt, *m_eta, *m_phi, *m_charge, *m_ptvarcone20, *m_topoetcone20,
                      leptons);

        build_jets(*j_pt, *j_eta, *j_phi, *j_m, *j_Jvt, *j_Width, *j_NumTrkPt500, jets);
        std::vector<MyJet> jets_sel = jets;

        // Fill jet pT histos
        if (jets_sel.size() > 1) {
            TLorentzVector j1, j2, jj;
            j1.SetPtEtaPhiM(jets_sel[0].pt, jets_sel[0].eta, jets_sel[0].phi, jets_sel[0].m);
            j2.SetPtEtaPhiM(jets_sel[1].pt, jets_sel[1].eta, jets_sel[1].phi, jets_sel[1].m);
            double dR_jj = j1.DeltaR(j2);
            if (dR_jj <= 2.0) {
                histos[sname]["jet_pt"]->Fill(jets_sel[0].pt, w_event);
                histos[sname]["jet_pt"]->Fill(jets_sel[1].pt, w_event);
                jj = j1 + j2;
                double m_jj = jj.M();
                histos[sname]["m_jj"]->Fill(m_jj, w_event);
            } 
        }

        // leptons
        if (leptons.size() > 1) {
            histos[sname]["lep1_pt"]->Fill(leptons[0].pt, w_event);
            histos[sname]["lep2_pt"]->Fill(leptons[1].pt, w_event);

            TH1F *hq = histos[sname]["dilep_q"];
            hq->GetXaxis()->SetBinLabel(1, "OS");
            hq->GetXaxis()->SetBinLabel(2, "");
            hq->GetXaxis()->SetBinLabel(3, "SS");
            float qprod = leptons[0].charge * leptons[1].charge;   //ouble
            if (qprod < 0) hq->Fill(-0.5, w_event);
            else if (qprod > 0) hq->Fill(+0.5, w_event);
        } else if (leptons.size() == 1) {
            histos[sname]["lep1_pt"]->Fill(leptons[0].pt, w_event);
        }
        // DR(l, jj) using the two leading jets as the dijet system
        if (jets_sel.size() >= 2 && !leptons.empty()) {
            TLorentzVector j1, j2, jj;
            j1.SetPtEtaPhiM(jets_sel[0].pt, jets_sel[0].eta, jets_sel[0].phi, jets_sel[0].m);
            j2.SetPtEtaPhiM(jets_sel[1].pt, jets_sel[1].eta, jets_sel[1].phi, jets_sel[1].m);
            double dR_jj = j1.DeltaR(j2);
            if (dR_jj <= 2.0) {
                jj = j1 + j2;
                // Leading lepton
                TLorentzVector l1;
                l1.SetPtEtaPhiM(leptons[0].pt, leptons[0].eta, leptons[0].phi, 0.0);
                double dR1 = l1.DeltaR(jj);
                histos[sname]["dr_lep1_jj"]->Fill(dR1, w_event);
                // Subleading lepton if it exists
                if (leptons.size() > 1) {
                    TLorentzVector l2;
                    l2.SetPtEtaPhiM(leptons[1].pt, leptons[1].eta, leptons[1].phi, 0.0);
                    double dR2 = l2.DeltaR(jj);
                    histos[sname]["dr_lep2_jj"]->Fill(dR2, w_event);
                }
            }
        }
        if (stype != "SIG") {
            fill_reco_variables(sname, leptons, jets_sel, w_event, false, histos);
        } else {
            fill_reco_variables(sname, leptons, jets_sel, w_event, true, histos);
        } 
        
    } // <-end of event loop
    std::cout<< "  #PV-rejected events: " << PVrejectedevents << " #False-eventClean events: " << EventUncleanevents << std::endl;
    // Dtermine the Weights from metadata (enabled by default) otherwise computed as sum of all event weights (also fallback if metadata fails)
    double total_sumw = 0.0;
    if (stype == "DATA") {
        std::cout<< "  #Rejected data: " << RejectedData << std::endl;
        // data: no metadata, no cross-section scaling
        total_sumw = total_sumw_loop;
    } else if (use_metadata_totalW && metadata) {
        double total_sumw_meta = get_total_sumw_from_metadata(chain_files, *metadata);
        if (total_sumw_meta > 0.0) {
            total_sumw = total_sumw_meta;
            std::cout << "  Using metadata total_sumw (sum AllExecutedEvents.W) = " << total_sumw << std::endl;
        } else {
            std::cerr << "  WARNING: metadata total_sumw <= 0 for " << sample.name << ", falling back to loop total_sumw.\n";
            total_sumw = total_sumw_loop;
        }
    } else {
        total_sumw = total_sumw_loop; 
    }
    double sigma_pb  = sample.sigma_pb;
    double lumi_pb   = sample.lumi_pb;
    double filt_eff  = sample.filter_eff;
    double k_factor  = sample.k_factor;
    double super_factor = 0.0;
    if (stype == "DATA") {
        super_factor = 1.0;    // do not scale data
    } else if (total_sumw != 0.0) {
        super_factor = sigma_pb * lumi_pb * filt_eff * k_factor / total_sumw;
    } else {
        std::cerr << "  WARNING: total_sumw = 0 for " << sname
                  << " (no valid events). Setting super_factor = 0.\n";
    }
    std::cout << "  total_sumw_loop = " << total_sumw_loop
              << ", total_sumw_used = " << total_sumw
              << ", super_factor = " << super_factor << std::endl;
    super_factors[sname] = super_factor;

    // Rescale all histograms of this sample 
    for (auto &kv : histos[sname]) {
        const std::string &hname = kv.first;
        if (hname == "dr_lep2_jj_unw") continue;
        TH1F *h = kv.second;
        if (stype != "DATA") h->Scale(super_factor);
    }
    delete chain;
}

// ============================================================================
// NEW: helper to compute and plot significance vs m(l2jj) for multiple SIG models
// ============================================================================

void make_significance_m_l2jj_canvas(
    const std::string &output_dir,
    const std::vector<Sample> &samples,
    const std::map<std::string, std::map<std::string, TH1F*>> &histos
) {
    // Collect signal, background, and data histograms for m_l2jj
    std::vector<std::string> sig_names;
    for (const auto &s : samples) {
        if (s.type == "SIG") {
            // only use if histogram exists
            auto itS = histos.find(s.name);
            if (itS != histos.end() && itS->second.count("m_l2jj")) {
                sig_names.push_back(s.name);
            }
        }
    }

    if (sig_names.empty()) {
        std::cerr << "No signal samples with m_l2jj histograms found. Nothing to plot.\n";
        return;
    }

    // Build stacked background histogram B(m_l2jj)
    TH1F *hB_template = nullptr;
    for (const auto &s : samples) {
        if (s.type.substr(0,3) == "BKG") {
            auto itB = histos.find(s.name);
            if (itB != histos.end() && itB->second.count("m_l2jj")) {
                hB_template = itB->second.at("m_l2jj");
                break;
            }
        }
    }
    if (!hB_template) {
        std::cerr << "No background m_l2jj histogram found. Cannot compute significance vs B.\n";
        return;
    }

    TH1F *hB = (TH1F*)hB_template->Clone("hB_m_l2jj");
    hB->Reset("ICES");
    for (const auto &s : samples) {
        if (s.type.substr(0,3) == "BKG") {
            auto itB = histos.find(s.name);
            if (itB == histos.end()) continue;
            auto itH = itB->second.find("m_l2jj");
            if (itH == itB->second.end()) continue;
            hB->Add(itH->second);
        }
    }

    // Build data histogram D(m_l2jj) if DATA is present
    TH1F *hD = nullptr;
    for (const auto &s : samples) {
        if (s.type == "DATA") {
            auto itD = histos.find(s.name);
            if (itD != histos.end() && itD->second.count("m_l2jj")) {
                if (!hD) {
                    hD = (TH1F*)itD->second.at("m_l2jj")->Clone("hD_m_l2jj");
                    hD->Reset("ICES");
                }
                hD->Add(itD->second.at("m_l2jj"));
            }
        }
    }
    bool haveData = (hD != nullptr);

    // Prepare Z histograms for each signal: vs B and vs DATA
    struct SigPlotInfo {
        std::string name;
        TH1F *Z_vsB;
        TH1F *Z_vsD;  // may be null if no data
        int color;
    };
    std::vector<SigPlotInfo> sig_plots;

    // Simple colour palette for up to a few models
    std::vector<int> colors = {kRed+1, kBlue+1, kMagenta+2, kGreen+2};
    int colorIndex = 0;

    for (const auto &sname : sig_names) {
        auto itSamp = histos.find(sname);
        if (itSamp == histos.end()) continue;
        auto itH = itSamp->second.find("m_l2jj");
        if (itH == itSamp->second.end()) continue;

        TH1F *hS_orig = itH->second;
        TH1F *hS = (TH1F*)hS_orig->Clone(("hS_"+sname+"_m_l2jj").c_str());

        TH1F *hZ_B = (TH1F*)hS->Clone(("hZ_"+sname+"_vsB").c_str());
        hZ_B->Reset("ICES");
        hZ_B->SetTitle(";m_{inv}(#it{l_{2}jj}) / GeV;asivmov sign. / #sigma");

        TH1F *hZ_D = nullptr;
        if (haveData) {
            hZ_D = (TH1F*)hS->Clone(("hZ_"+sname+"_vsD").c_str());
            hZ_D->Reset("ICES");
            hZ_D->SetTitle(";m_{inv}(#it{l_{2}jj}) / GeV;asivmov sign. / #sigma");
        }

        int nbins = hS->GetNbinsX();
        double b_min = 1e-3;

        for (int ibin = 1; ibin <= nbins; ++ibin) {
            double s = hS->GetBinContent(ibin);

            // Undo plotting buff used only for SIG in fill_reco_variables:
            s /= 5e4;

            double b = hB->GetBinContent(ibin);
            double zB = 0.0;
            if (b >= b_min && s > 0) {
                double arg = 1.0 + s / b;
                if (arg > 0) {
                    double val = 2.0 * ((s+b) * std::log(arg) - s);
                    if (val > 0) zB = std::sqrt(val);
                }
            }
            hZ_B->SetBinContent(ibin, zB);

            if (haveData && hZ_D) {
                double d = hD->GetBinContent(ibin);
                double zD = 0.0;
                if (d >= b_min && s > 0) {
                    double argD = 1.0 + s / d;
                    if (argD > 0) {
                        double valD = 2.0 * ((s+d) * std::log(argD) - s);
                        if (valD > 0) zD = std::sqrt(valD);
                    }
                }
                hZ_D->SetBinContent(ibin, zD);
            }
        }

        SigPlotInfo info;
        info.name  = sname;
        info.Z_vsB = hZ_B;
        info.Z_vsD = hZ_D;
        info.color = colors[colorIndex % colors.size()];
        ++colorIndex;

        sig_plots.push_back(info);
    }

    if (sig_plots.empty()) {
        std::cerr << "No valid signal Z(m_l2jj) histograms to plot.\n";
        return;
    }

    // Create canvas and draw all lines
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c_m_l2jj_significance",
                             "Significance vs m(l_{2}jj)",
                             1000, 800);
    c->cd();
    gPad->SetLogx(kFALSE);
    gPad->SetLogy(kTRUE);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetGrid(0,0);

    // Use first signal vsB hist as axis template
    TH1F *hAxis = sig_plots.front().Z_vsB;
    hAxis->GetYaxis()->SetTitle("asimov sign. / #sigma");
    hAxis->SetLineColor(kWhite);
    hAxis->SetMarkerSize(0);
    // hAxis->SetMaximum(0.0);
    // for (const auto &info : sig_plots) {
    //     if (!info.Z_vsB) continue;
    //     double localMax = info.Z_vsB->GetMaximum();
    //     if (info.Z_vsD) {
    //         localMax = std::max(localMax, info.Z_vsD->GetMaximum());
    //     }
    //     if (localMax > hAxis->GetMaximum()) hAxis->SetMaximum(1.2 * localMax);
    // }
    // hAxis->SetMinimum(0.0);
    hAxis->GetXaxis()->SetRangeUser(20.0, 600.0);  
    hAxis->SetMinimum(1.5e-6);                     
    hAxis->SetMaximum(1.5e0);  


    hAxis->GetXaxis()->SetNdivisions(510);
    hAxis->GetXaxis()->SetMoreLogLabels(kTRUE);
    hAxis->GetXaxis()->SetNoExponent(kTRUE);
    hAxis->GetXaxis()->SetTitleOffset(1.0);
    hAxis->GetYaxis()->SetTitleOffset(1.1);
    hAxis->GetXaxis()->SetTitleSize(0.045);
    hAxis->GetYaxis()->SetTitleSize(0.045);
    hAxis->GetXaxis()->SetLabelSize(0.04);
    hAxis->GetYaxis()->SetLabelSize(0.04);
    hAxis->GetYaxis()->SetTitleFont(32);
    hAxis->GetXaxis()->SetTitleFont(32);

    hAxis->Draw("HIST");

    // Legend
    TLegend *leg = new TLegend(0.15, 0.60, 0.5, 0.68);
    leg->SetBorderSize(0);
    leg->SetNColumns(2); 
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    // Draw each signal with two line styles
    for (const auto &info : sig_plots) {
        if (info.Z_vsB) {
            info.Z_vsB->SetLineColor(info.color);
            info.Z_vsB->SetLineWidth(3);
            info.Z_vsB->SetLineStyle(1); // solid: vs background
            info.Z_vsB->SetMarkerStyle(0);
            info.Z_vsB->Draw("HIST SAME");
            if (info.name == "WR5000_N500_2L2J") {
                leg->AddEntry(info.Z_vsB, "model#it{B}","l");
            } else {
                leg->AddEntry(info.Z_vsB, "model#it{A*}","l");
            }
            
        }
        // if (info.Z_vsD && haveData) {
        //     info.Z_vsD->SetLineColor(info.color);
        //     info.Z_vsD->SetLineWidth(3);
        //     info.Z_vsD->SetLineStyle(2); // dashed: vs data
        //     info.Z_vsD->SetMarkerStyle(0);
        //     info.Z_vsD->Draw("HIST SAME");
        //     leg->AddEntry(info.Z_vsD,
        //                   (info.name + " vs data").c_str(),
        //                   "l");
        // }
    }

    // ATLAS-style label
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(72);
    lat.SetTextSize(0.042);
    lat.DrawLatex(0.15, 0.82, "ATLAS");  // 0.85->2
    lat.SetTextFont(42);
    lat.DrawLatex(0.27, 0.82, "Internal");
    lat.SetTextFont(32);
    lat.DrawLatex(0.15, 0.76, "pp #rightarrow W_{R} #rightarrow l (N #rightarrow l j j )");
    lat.DrawLatex(0.15, 0.70, "#sqrt{s} = 13.6 TeV, #scale[0.7]{#int} L = 500 fb^{-1}");

    leg->Draw();
    gPad->RedrawAxis();

    std::string outname = output_dir + "/significance_m_l2jj_multiSIG_vsB_vsDATA.png";
    c->SaveAs(outname.c_str());

    std::cout << "Saved significance plot to " << outname << std::endl;
}

// ============================================================================
// main macro
// ============================================================================

void SIGvsBKG_sgn22() { 
    // silence ROOT noise: keep only serious errors
    gErrorIgnoreLevel = kFatal; //kError
    bool use_metadata_totalW = true;
    std::cout << "User Test 006:" << std::endl;
    std::string output_dir = "/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/SIGvsBKG_RECO/Figures/SIGvsBKG_SGN";
    std::string metadata_json_path = "/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/SIGvsBKG_RECO/metadataLLP1_backgrounds.json";
    json metadata;
    if (use_metadata_totalW) {
        std::ifstream jf(metadata_json_path);
        if (!jf) {
            std::cerr << "ERROR: cannot open metadata JSON: "
                      << metadata_json_path << "\n";
            use_metadata_totalW = false; // fall back if file missing
        } else {
            jf >> metadata;
        }
    }
    // Define samples
    /*
    MCprofile: CrossSection[pb], Int.Lum.[pb^-], FilterEff, k-factor
    Signal: CS(50GeV, 5TeV, 0.1)=4.072e-4pb 
            Int. Lum: 29.0493 fb^- = 29049.3 pb^-
            Filter Eff: MG (#tau events: 6%) 
    Signal: CS(70GeV, 5TeV, 0.1)=3.968-4pb  
    Signal: CS(500GeV, 5TeV, 0.1)=2.32479e-4pb  
    Signal: CS(500GeV, 10TeV, 0.1)=4.672e-6pb     (e-15mb)
    Signal: CS(50GeV, 10TeV, 0.1)=1.303-5pb  
    */
    std::vector<Sample> samples;
    {
        Sample s;
        s.name = "WR5000_N500_2L2J";
        s.type = "SIG";
        s.sigma_pb    = 2.32479e-4;  //3.968e-4 (70, 5000) 4.072e-4 (50,5000), 2.32479e-4(500, 5000)
        s.lumi_pb     = 500000;    //29049.3   500000  3549
        s.filter_eff  = 0.94;
        s.k_factor    = 1.0;
        s.files = {
            "/srv/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_500_2L2J_MadSpin.DAOD_LLP1.pool.root"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "WR10000_N50_2L2J";
        s.type = "SIG";
        s.sigma_pb    = 1.303e-5;  //3.968e-4 (70, 5000) 4.072e-4 (50,5000), 2.32479e-4(500, 5000)
        s.lumi_pb     = 500000;    //29049.3   500000  3549
        s.filter_eff  = 0.94;
        s.k_factor    = 1.0;
        s.files = {
            "/srv/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_10000_50_2L2J_MadSpin.DAOD_LLP1.pool.root"
        };
        samples.push_back(s);
    }
    // {
    //     Sample s;
    //     s.name = "Data_2022";
    //     s.type = "DATA";
    //     s.sigma_pb    = 1.0;       
    //     s.lumi_pb     = 1.0;       
    //     s.filter_eff  = 1.0;
    //     s.k_factor    = 1.0;
    //     s.files = {
    //         "/srv/data22_13p6TeV.periodJ/ALL/data22_13p6TeV.periodJ.physics_Main.PhysCont.DAOD_LLP1.grp22_v01_p7079/DAOD_LLP1.*.pool.root.1"
    //     };
    //     samples.push_back(s);
    // }
    {
        Sample s;
        s.name = "Zmumu";
        s.type = "BKG_Z";
        s.sigma_pb    = 1998.8;
        s.lumi_pb     = 500000;
        s.filter_eff  = 1.0;
        s.k_factor    = 1.0;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zmumu_CVBV";
        s.type = "BKG_Z";
        s.sigma_pb    = 2335.9;
        s.lumi_pb     = 500000;
        s.filter_eff  = 8.444979e-1;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zmumu_BF";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 500000;
        s.filter_eff  = 2.502305e-2;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zmumu_CFBV";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 500000;
        s.filter_eff  = 1.304099e-1;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zee_CVBV";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 500000;
        s.filter_eff  = 8.442598e-1;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zee_BF";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 500000;
        s.filter_eff  = 2.558475e-2;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zee_CFBV";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 500000;
        s.filter_eff  = 1.301249e-1;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "ttbar_llqq";
        s.type = "BKG_T";
        s.sigma_pb    = 811.29;
        s.lumi_pb     = 500000;
        s.filter_eff  = 4.384566e-1;
        s.k_factor    = 1.138433852;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "ttbar_llll";
        s.type = "BKG_T";
        s.sigma_pb    = 85.482;
        s.lumi_pb     = 500000;
        s.filter_eff  = 1;
        s.k_factor    = 1.138433852;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "ttbar_qqqq";
        s.type = "BKG_T";
        s.sigma_pb    = 811.29;
        s.lumi_pb     = 500000;
        s.filter_eff  = 4.561727e-1;
        s.k_factor    = 1.1384339;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
        };
        samples.push_back(s);
    }

    // Histograms per sample
    std::map<std::string, std::map<std::string, TH1F*>> histos;
    std::map<std::string, double> super_factors;

    // calculate Progressive bins
    static std::vector<double> lep1_bins = make_log_bins(lep_pt_min, 3950.0, 20); //3950
    static std::vector<double> lep2_bins = make_log_bins(lep_pt_min, 3050.0, 18); //2950
    static std::vector<double> jet_bins  = make_log_bins(jet_pt_min, 3350.0, 20); //3350
    static std::vector<double> disp_bins = make_log_bins(0.01, 100.0, 25);   //0.5, 250
    static std::vector<double> m_lljj_bins = make_log_bins(95, 6950, 28);  //28
    // static std::vector<double> m_l2jj_bins = make_log_bins(20, 690, 15);  //28
    
    for (auto &s : samples) {
        const std::string &sname = s.name;
        const char *y_title_latex = "event weight: #sigma #scale[0.7]{#int} #font[12]{L} #varepsilon k w^{i} / W_{tot}";
        const char *y_title_latex_ljj = "events / 35 GeV";
        const char *y_title_latex_lljj = "events / 343 GeV";
        int n_lep1 = lep1_bins.size() - 1;
        histos[sname]["lep1_pt"] = new TH1F(
            ("h_"+sname+"_lep1_pt").c_str(),
            (std::string(" ;pT(#it{l_{1}} ) / GeV;") + y_title_latex).c_str(), 
            n_lep1, lep1_bins.data()
        );
        int n_lep2 = lep2_bins.size() - 1;
        histos[sname]["lep2_pt"] = new TH1F(
            ("h_"+sname+"_lep2_pt").c_str(),
            (std::string(" ;pT(#it{l_{2}} ) / GeV;") + y_title_latex).c_str(),
            n_lep2, lep2_bins.data()
        );
        int n_jet = jet_bins.size() - 1;
        histos[sname]["jet_pt"] = new TH1F(
            ("h_"+sname+"_jet_pt").c_str(),
            (std::string(" ;pT(#it{j} ) / GeV;") + y_title_latex).c_str(),
            n_jet, jet_bins.data()
        );
        int n_disp = disp_bins.size() - 1;
        histos[sname]["disp"] = new TH1F(
            ("h_"+sname+"_disp").c_str(),
            (std::string(" ;Lxyz / mm;") + y_title_latex).c_str(),
            n_disp, disp_bins.data()
        );
        double edges[4] = {-0.75, -0.25, 0.25, 0.75};
        histos[sname]["dilep_q"] = new TH1F(
            ("h_"+sname+"_dilep_q").c_str(),
            (std::string(" ;;") + y_title_latex).c_str(),
            3, edges
        );
        histos[sname]["dr_lep2_jj"] = new TH1F(
            ("h_"+sname+"_dr_lep2_jj").c_str(),
            (std::string(" ;#DeltaR(#it{l_{2}, jj} );") + y_title_latex).c_str(),
            20, 0.0, 4.5
        );
        histos[sname]["dr_lep1_jj"] = new TH1F(
            ("h_"+sname+"_dr_lep1_jj").c_str(),
            (std::string(" ;#DeltaR(#it{l_{1}, jj} );") + y_title_latex).c_str(),
            20, 0.0, 5.5
        );
        // int n_m_l2jj = m_l2jj_bins.size() - 1;
        histos[sname]["m_l2jj"] = new TH1F(
            ("h_"+sname+"_m_l2jj").c_str(),
            (std::string(" ;m_{inv.}(#it{l_{2}jj}) / GeV;") + y_title_latex).c_str(),
            24, 20.0, 600.0   //25!!!!
            // n_m_l2jj, m_l2jj_bins.data()
        );
        int n_m_lljj = m_lljj_bins.size() - 1;
        histos[sname]["m_lljj"] = new TH1F(
            ("h_"+sname+"_m_lljj").c_str(),
            (std::string(" ;m_{inv.}(#it{lljj}) / GeV;") + y_title_latex).c_str(),
            n_m_lljj, m_lljj_bins.data()
        );
        histos[sname]["m_lljj_sel"] = new TH1F(
            ("h_"+sname+"_m_lljj_sel").c_str(),
            (std::string(" ;m_{inv.}(#it{lljj}) / GeV;") + y_title_latex).c_str(),
            n_m_lljj, m_lljj_bins.data()
        );
        histos[sname]["m_jj"] = new TH1F(
            ("h_"+sname+"_m_jj").c_str(),
            (std::string(" ;m_{inv.}(#it{jj}) / GeV;") + y_title_latex).c_str(),
            16, 50, 790
        );
        if (s.type == "SIG") {
            int n_disp_sig = disp_bins.size() - 1;
            histos[sname]["disp_reco_d0"] = new TH1F(
                ("h_"+sname+"_disp_reco_d0").c_str(),
                (std::string(" ;L_{reco}^{d0} / mm;") + y_title_latex).c_str(),
                n_disp_sig, disp_bins.data()
            );
            histos[sname]["disp_truth_Lxyz"] = new TH1F(
                ("h_"+sname+"_disp_truth_Lxyz").c_str(),
                (std::string(" ;L_{truth}^{xyz} / mm;") + y_title_latex).c_str(),
                n_disp_sig, disp_bins.data()
            );
        }
        TH1::SetDefaultSumw2(true);
    }

    for (auto &s : samples) {
        process_sample(s, histos, super_factors,
                       use_metadata_totalW ? &metadata : nullptr,
                       use_metadata_totalW);
    }

    make_significance_m_l2jj_canvas(output_dir, samples, histos);
    std::cout << "\nAnalysis successfully completed. Results saved to " << output_dir << std::endl;
}
