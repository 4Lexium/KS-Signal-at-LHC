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
    int charge;
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
                   const std::vector<int>   &e_charge,
                   const std::vector<float> &e_ptvarcone20,
                   const std::vector<float> &e_topoetcone20,
                   const std::vector<float> &m_pt,
                   const std::vector<float> &m_eta,
                   const std::vector<float> &m_phi,
                   const std::vector<int>   &m_charge,
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
                const std::vector<float> &j_NumTrkPt500,
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
        if (j_NumTrkPt500[i] > jNumTrkPt500_cut) continue;
        MyJet jet;
        jet.pt = j_pt[i] * 0.001f;
        jet.eta = j_eta[i];
        jet.phi = j_phi[i];
        jet.m = j_m[i] * 0.001f;
        jet.Jvt = j_Jvt[i];
        jet.Width = j_Width[i];
        jet.NumTrkPt500 = j_NumTrkPt500[i];
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
    int qprod = leptons[0].charge*leptons[1].charge;
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
            itHist_lljj->second->Fill(m_lljj, 50);  //1e2  5e1   //2e3
        } else {
            itHist_lljj->second->Fill(m_lljj, 1);    //w_event*2e5
        }
    }
    double m_jj = jj.M();
    auto itHist_jj = varMap.find("m_jj");
    if (itHist_jj != varMap.end()) {
        if (buff_weight) {
            itHist_jj->second->Fill(m_jj, w_event*1e-2);  //1e3
        } else {
            itHist_jj->second->Fill(m_jj, w_event);
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
            itHist_l2jj->second->Fill(m_l2jj, 50);   //1e3  5e3
        } else {
            itHist_l2jj->second->Fill(m_l2jj, 1);   //w_event*2e5
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

    chain->SetBranchStatus("ElectronsAuxDyn.pt",           1);
    chain->SetBranchStatus("ElectronsAuxDyn.eta",          1);
    chain->SetBranchStatus("ElectronsAuxDyn.phi",          1);
    chain->SetBranchStatus("ElectronsAuxDyn.charge",       1);
    chain->SetBranchStatus("ElectronsAuxDyn.ptvarcone20",  1);
    chain->SetBranchStatus("ElectronsAuxDyn.topoetcone20", 1);

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
    std::vector<int>   *e_charge       = nullptr;
    std::vector<float> *e_ptvarcone20  = nullptr;
    std::vector<float> *e_topoetcone20 = nullptr;

    std::vector<float> *m_pt           = nullptr;
    std::vector<float> *m_eta          = nullptr;
    std::vector<float> *m_phi          = nullptr;
    std::vector<int>   *m_charge       = nullptr;
    std::vector<float> *m_ptvarcone20  = nullptr;
    std::vector<float> *m_topoetcone20 = nullptr;

    std::vector<float> *j_pt           = nullptr;
    std::vector<float> *j_eta          = nullptr;
    std::vector<float> *j_phi          = nullptr;
    std::vector<float> *j_m            = nullptr;
    std::vector<float> *j_Jvt          = nullptr;
    std::vector<float> *j_Width        = nullptr;
    std::vector<float> *j_NumTrkPt500  = nullptr;

    std::vector<short> *PV_types         = nullptr;
    int eventClean = 0;
    std::vector<float> *mcEventWeights = nullptr;
    
    chain->SetBranchAddress("ElectronsAuxDyn.pt",           &e_pt);
    chain->SetBranchAddress("ElectronsAuxDyn.eta",          &e_eta);
    chain->SetBranchAddress("ElectronsAuxDyn.phi",          &e_phi);
    chain->SetBranchAddress("ElectronsAuxDyn.charge",       &e_charge);
    chain->SetBranchAddress("ElectronsAuxDyn.ptvarcone20",  &e_ptvarcone20);
    chain->SetBranchAddress("ElectronsAuxDyn.topoetcone20", &e_topoetcone20);

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

    // std::vector<UInt_t> *larFlags  = nullptr;
    // std::vector<UInt_t> *tileFlags = nullptr;
    // std::vector<UInt_t> *sctFlags  = nullptr;
    // std::vector<UInt_t> *coreFlags = nullptr;
    // std::vector<UInt_t> *rnums     = nullptr;
    // std::vector<UInt_t> *lbs       = nullptr;
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

        // if (ievt > 5000 && ievt < 5010) {
        //     std::cout << "Event " << ievt << "  ";
        //     if (!PV_types) {
        //         std::cout << "PV_types is NULL\n";
        //     } else {
        //         std::cout << "vertexType size = " << PV_types->size() << " values:";
        //         for (size_t i = 0; i < PV_types->size(); ++i) {
        //             std::cout << " " << (*PV_types)[i];
        //         }
        //         std::cout << "\n";
        //     }
        // }

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
        build_leptons(*e_pt, *e_eta, *e_phi, *e_charge, *e_ptvarcone20, *e_topoetcone20,
                      *m_pt, *m_eta, *m_phi, *m_charge, *m_ptvarcone20, *m_topoetcone20,
                      leptons);
        build_jets(*j_pt, *j_eta, *j_phi, *j_m, *j_Jvt, *j_Width, *j_NumTrkPt500, jets);
        std::vector<MyJet> jets_sel = jets;
        // selection done innside build function
        // jets_sel.reserve(jets.size());
        // for (auto &j : jets) {
        //     if (j.pt > jet_pt_min) jets_sel.push_back(j);
        // }
        // 
        // Fill jet pT histos
        if (jets_sel.size() > 1) {
            TLorentzVector j1, j2, jj;
            j1.SetPtEtaPhiM(jets_sel[0].pt, jets_sel[0].eta, jets_sel[0].phi, jets_sel[0].m);
            j2.SetPtEtaPhiM(jets_sel[1].pt, jets_sel[1].eta, jets_sel[1].phi, jets_sel[1].m);
            double dR_jj = j1.DeltaR(j2);
            if (dR_jj <= 2.0) {
                histos[sname]["jet_pt"]->Fill(jets_sel[0].pt, w_event);
                histos[sname]["jet_pt"]->Fill(jets_sel[1].pt, w_event);
            }
        }
        // } else if (jets_sel.size() == 1) {
        //     // histos[sname]["jet_pt"]->Fill(jets_sel[0].pt, w_event);
        // }

        // leptons
        if (leptons.size() > 1) {
            histos[sname]["lep1_pt"]->Fill(leptons[0].pt, w_event);
            histos[sname]["lep2_pt"]->Fill(leptons[1].pt, w_event);

            TH1F *hq = histos[sname]["dilep_q"];
            hq->GetXaxis()->SetBinLabel(1, "OS");
            hq->GetXaxis()->SetBinLabel(2, "");
            hq->GetXaxis()->SetBinLabel(3, "SS");
            int qprod = leptons[0].charge * leptons[1].charge;
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
            fill_reco_variables(sname, leptons, jets_sel, w_event*1e2, true, histos);
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
        if (hname == "dr_lep2_jj_unw" || hname == "m_l2jj" || hname == "m_lljj") continue;   //hname == "dr_lep2_jj_unw" debug to 
        TH1F *h = kv.second;
        if (stype != "DATA") h->Scale(super_factor);
    }
    delete chain;
}

// Combine samples into groups, compute significance, and plot
void make_sig_canvas(const std::string &var,
                     const std::string &title_top,
                     const std::string &outname,
                     bool logy_top,
                     bool draw_text_bottom,
                     const std::map<std::string, TH1F*> &group_histo_var,
                     TH1F *hZ)
{
    gStyle->SetOptStat(0);
    TGaxis::SetMaxDigits(3);   //4
    TCanvas *c = new TCanvas(("c_"+var).c_str(), title_top.c_str(), 1000, 800);
    c->Divide(1, 2);
    bool logx = (var != "dilep_q" && var != "dr_lep1_jj" && var != "dr_lep2_jj" && var != "m_l2jj" && var != "m_jj"); //&& var != "m_lljj"

    // ----------------- TOP PAD -----------------
    c->cd(1);
    gPad->SetPad(0.0, 0.30, 1.0, 1.0);
    gPad->SetBottomMargin(0.02);
    if (logy_top) gPad->SetLogy();
    if (logx)     gPad->SetLogx();
    gPad->SetTicky(1);
    // Colours
    std::map<std::string, int> group_colors;
    group_colors["Zmumu"] = kRed+2;
    group_colors["Zee"]   = kBlue+2;
    group_colors["ttbar"] = kGreen+3;
    group_colors["SIG"]   = kOrange+1;

    // legend position
    double x1, y1, x2, y2;
    if (var != "disp") {
        x1 = 0.76; y1 = 0.65;
        x2 = 0.94; y2 = 0.90;
    } else {
        x1 = 0.37; y1 = 0.68;
        x2 = 0.52; y2 = 0.90;
    }
    TLegend *leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);

    // define drawing order: backgrounds (stacked) + signal (overlaid) 
    std::vector<std::string> bkg_order = {"ttbar", "Zmumu", "Zee"};
    std::string sig_name = "SIG";

    // compute histogram limits
    double HISTMAX = 0.0;
    std::vector<std::string> group_order = {"ttbar", "Zmumu", "Zee", "SIG"};
    for (auto &gname : group_order) {
        auto it = group_histo_var.find(gname);
        if (it == group_histo_var.end()) continue;
        HISTMAX = std::max(HISTMAX, (double)it->second->GetMaximum());
    }
    double xmin_vis = 0.0;
    double xmax_vis = 0.0;
    if (var == "lep1_pt" || var == "lep2_pt" || var == "jet_pt") {
        xmin_vis = lep_pt_min;
        xmax_vis = group_histo_var.begin()->second->GetXaxis()->GetXmax();
    } else if (var == "disp") {
        xmin_vis = 0.3;
        xmax_vis = group_histo_var.begin()->second->GetXaxis()->GetXmax();
    } else if (var == "dilep_q") {
        xmin_vis = -1.0;
        xmax_vis =  1.0;
    } else {
        xmin_vis = group_histo_var.begin()->second->GetXaxis()->GetXmin();
        xmax_vis = group_histo_var.begin()->second->GetXaxis()->GetXmax();
    }

    // build THStack for backgrounds 
    THStack *hs_bkg = new THStack(("hs_"+var).c_str(), "");
    double max_bkg = 0.0;
    double max_sig = 0.0;
    for (auto &gname : bkg_order) {
        auto it = group_histo_var.find(gname);
        if (it == group_histo_var.end()) continue;
        TH1F *h = it->second;
        int color = group_colors[gname];
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
        h->SetFillColor(color);
        h->SetFillStyle(1001);
        h->SetMarkerStyle(0);
        max_bkg = std::max(max_bkg, (double)h->GetMaximum());
        hs_bkg->Add(h, "HIST");
        leg->AddEntry(h, gname.c_str(), "f");
    }

    // signal histogram on top
    TH1F *h_sig = nullptr;
    {
        auto it_sig = group_histo_var.find(sig_name);
        if (it_sig != group_histo_var.end()) {
            h_sig = it_sig->second;
            int color = group_colors[sig_name];
            h_sig->SetLineColor(color);
            h_sig->SetLineWidth(4);
            h_sig->SetLineStyle(1); 
            h_sig->SetFillColor(color);
            h_sig->SetFillStyle(3004);
            h_sig->SetMarkerStyle(20);
            h_sig->SetMarkerSize(0.8);
            h_sig->SetMarkerColor(color);
            max_sig = h_sig->GetMaximum();
            leg->AddEntry(h_sig, "model#it{B}", "l");
            if (var == "m_l2jj" || var == "m_l2jj" ) {
                leg->AddEntry((TObject*)0, " ", " ");
            }
        }
    }

    TH1F *h_axis = nullptr;
    if (!bkg_order.empty()) {
        auto it = group_histo_var.find(bkg_order.front());
        if (it != group_histo_var.end()) h_axis = it->second;
    }
    if (!h_axis && h_sig) h_axis = h_sig;

    if (h_axis) {
        h_axis->GetXaxis()->SetLabelSize(0.0);
        h_axis->GetXaxis()->SetTitleSize(0.0);
        h_axis->GetYaxis()->SetTitleOffset(0.95);   //0.8
        h_axis->GetYaxis()->SetTitleSize(0.042);
        h_axis->GetYaxis()->SetLabelSize(0.040);
        h_axis->GetYaxis()->SetNoExponent(kFALSE);
        h_axis->GetYaxis()->SetTitleFont(32);
        h_axis->GetXaxis()->SetTitleFont(32);
        h_axis->GetXaxis()->SetNdivisions(510);
        h_axis->GetXaxis()->SetMoreLogLabels(kTRUE);
    }
    hs_bkg->Draw("HIST");
    // Hardcoded limits, chosen for aestetic purpouses 
    if (var == "m_l2jj") {
        hs_bkg->SetMaximum(4.99e3);   //1.19e4 for 500fb: 1.99e5  for 3fb:1.39e3  // FNAL: 2.75e3
    } else if (var == "m_jj") {
        hs_bkg->SetMaximum(1.5e4);    //1.5e4   for 500fb:1.5e5   for 3fb:1.5e3
        hs_bkg->SetMinimum(1.5e-4);   //1.5e-3  for 500fb:1.5e-2   for 3fb:1.5e-4
    } else if (var == "m_lljj" ) {
        hs_bkg->SetMaximum(4.99e3);   //1.19e4 for 500fb:1.99e5    for 3fb:1.39e3   // FINAL: 2.75e3
    } else if (var == "dilep_q") {
        hs_bkg->SetMaximum(1.5e10);
        hs_bkg->SetMinimum(1.5e-4);   //1.5e-3  for 500fb:1.5e-2   for 3fb:1.5e-4
    } else if (var == "jet_pt") {
        hs_bkg->SetMaximum(1.5e7);    //1.5e7  for 500fb:1.5e8    for 3fb:1.5e6
        hs_bkg->SetMinimum(1.5e-4);
    } else if (var == "lep1_pt") {
        hs_bkg->SetMaximum(1.5e8);    //1.5e7  for 500fb:1.5e8    for 3fb:1.5e6
        hs_bkg->SetMinimum(1.5e-4);
    } else {
        hs_bkg->SetMaximum(1.5e6);    //1.5e7  for 500fb:1.5e8    for 3fb:1.5e6
        hs_bkg->SetMinimum(1.5e-4);   //1.5e-3  for 500fb:1.5e-2   for 3fb:1.5e-4
    }

    hs_bkg->GetXaxis()->SetRangeUser(xmin_vis, xmax_vis);
    // re‑apply axis titles and fonts 
    if (h_axis) {
        hs_bkg->GetXaxis()->SetTitle(h_axis->GetXaxis()->GetTitle());
        hs_bkg->GetYaxis()->SetTitle(h_axis->GetYaxis()->GetTitle());
        hs_bkg->GetXaxis()->SetLabelSize(0.0);
        hs_bkg->GetXaxis()->SetTitleSize(0.0);
        hs_bkg->GetXaxis()->SetTickLength(0.05);
        hs_bkg->GetYaxis()->SetTitleOffset(0.95);   //0.8
        hs_bkg->GetYaxis()->SetTitleSize(0.042);
        hs_bkg->GetYaxis()->SetLabelSize(0.04);
        hs_bkg->GetYaxis()->SetTitleFont(32);
        hs_bkg->GetXaxis()->SetTitleFont(32);
    }

    // draw signal on top
    if (h_sig) h_sig->Draw("HIST SAME");
    // draw data as error bars
    TH1F *h_data = nullptr;
    auto it_data = group_histo_var.find("DATA");
    if (it_data != group_histo_var.end()) {
        h_data = it_data->second;
        h_data->SetLineColor(kBlack);
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(1.0);
        h_data->SetFillStyle(0);
        h_data->Draw("PE0 SAME");
        leg->AddEntry(h_data, "data22", "pe");
    }
    leg->Draw();
    // ATLAS label
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.SetTextFont(72);
    if (var != "disp") {
        lat.DrawLatex(0.55, 0.85, "ATLAS");
        lat.SetTextFont(42);
        lat.DrawLatex(0.63, 0.85, "Internal");
        lat.SetTextFont(32);
        lat.DrawLatex(0.52, 0.8, "pp #rightarrow W_{R} #rightarrow l (N #rightarrow l j j )");
        lat.SetTextFont(32);
        lat.DrawLatex(0.51, 0.73, "#sqrt{s} = 13.6 TeV, #scale[0.7]{#int} L = 3.55 fb^{-1}");   //29.0  3.55
    } else {
        lat.DrawLatex(0.16, 0.85, "ATLAS");
        lat.SetTextFont(42);
        lat.DrawLatex(0.24, 0.85, "Internal");
        lat.SetTextFont(32);
        lat.DrawLatex(0.13, 0.8, "pp #rightarrow W_{R} #rightarrow l (N #rightarrow l j j )");
        lat.SetTextFont(32);
        lat.DrawLatex(0.12, 0.73, "#sqrt{s} = 13.6 TeV, #scale[0.7]{#int} L = 3.55 fb^{-1}");
    }
    gPad->RedrawAxis();

    // Draw Bottom PAD (significance) 
    c->cd(2);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetPad(0.0, 0.0, 1.0, 0.30);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.35);
    gPad->SetGridy(kTRUE);
    gPad->SetLogy(kFALSE);
    if (logx) gPad->SetLogx();
    hZ->SetStats(0);
    hZ->SetMinimum(0.0);
    hZ->SetMaximum(1.2 * hZ->GetMaximum());
    hZ->GetXaxis()->SetRangeUser(xmin_vis, xmax_vis);
    hZ->SetLineColor(kBlack);
    hZ->SetLineWidth(2);
    hZ->SetFillColor(0);
    hZ->SetMarkerStyle(20);
    hZ->SetMarkerSize(0.8);
    auto xaxis = hZ->GetXaxis();
    auto yaxis = hZ->GetYaxis();
    xaxis->SetTitleSize(0.11);
    xaxis->SetTitleOffset(1.2);
    xaxis->SetNdivisions(510);
    xaxis->SetMoreLogLabels(kTRUE);
    xaxis->SetNoExponent(kTRUE);
    xaxis->SetLabelSize(0.10);
    TGaxis::SetExponentOffset(-0.045, -0.03, "y");   //0.03
    yaxis->SetTitle("azimov sign. /#sigma");
    yaxis->SetTitleFont(32);
    xaxis->SetTitleFont(32);
    yaxis->SetTitleSize(0.10);
    yaxis->SetMaxDigits(3);
    yaxis->SetTitleOffset(0.4);  //0.4
    yaxis->SetLabelSize(0.08);
    yaxis->SetNdivisions(505);  //505
    hZ->Draw("PE1");
    gPad->Modified();
    gPad->Update();
    gPad->RedrawAxis();
    gPad->Update();
    c->SaveAs(outname.c_str());
}

// main macro
void SIGvsBKG_evt22() { 
    // silence ROOT noise: keep only serious errors
    gErrorIgnoreLevel = kFatal; //kError
    bool use_metadata_totalW = true;
    std::cout << "User Test 069:" << std::endl;
    std::string output_dir = "/home/alexanum/WORKSPACE/testing/MadGraph5/analysis/SIGvsBKG_RECO/Figures/SIGvsBKG_reco2_3fb_nodat";
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
        s.name = "WR5000_N50_2L2J";
        s.type = "SIG";
        s.sigma_pb    = 2.32479e-4;  //3.968e-4 (70, 5000) 4.072e-4 (50,5000), 2.32479e-4(500, 5000)
        s.lumi_pb     = 3549;    //29049.3   500000  3549
        s.filter_eff  = 0.94;
        s.k_factor    = 1.0;
        s.files = {
            "/srv/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_500_2L2J_MadSpin.DAOD_LLP1.pool.root"
            // "/srv/ROOT_files/DAOD_LLP1.DAOD_LLP1.mc23_13p6TeV.MGPy8EG_LRSM_WR_HNL_5000_500_2L2J_MadSpin.DAOD_LLP1.pool.root"
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
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000001.pool.root.1",
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000002.pool.root.1",
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000003.pool.root.1",
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000004.pool.root.1",
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000005.pool.root.1",
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000006.pool.root.1",
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000007.pool.root.1",
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000008.pool.root.1",
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000009.pool.root.1",
    //         // "/srv/data22_13p6TeV.periodJ/DAOD_LLP1.47616571._000010.pool.root.1",
    //     };
    //     samples.push_back(s);
    // }
    {
        Sample s;
        s.name = "Zmumu";
        s.type = "BKG_Z";
        s.sigma_pb    = 1998.8;
        s.lumi_pb     = 3549;
        s.filter_eff  = 1.0;
        s.k_factor    = 1.0;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.601190.PhPy8EG_AZNLO_Zmumu.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762661._000102.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zmumu_CVBV";
        s.type = "BKG_Z";
        s.sigma_pb    = 2335.9;
        s.lumi_pb     = 3549;
        s.filter_eff  = 8.444979e-1;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.700791.Sh_2214_Zmumu_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762549._000584.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zmumu_BF";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 3549;
        s.filter_eff  = 2.502305e-2;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.700789.Sh_2214_Zmumu_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762363._000036.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zmumu_CFBV";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 3549;
        s.filter_eff  = 1.304099e-1;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.700790.Sh_2214_Zmumu_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762751._000027.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zee_CVBV";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 3549;
        s.filter_eff  = 8.442598e-1;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.700788.Sh_2214_Zee_maxHTpTV2_CVetoBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762860._000141.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zee_BF";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 3549;
        s.filter_eff  = 2.558475e-2;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.700786.Sh_2214_Zee_maxHTpTV2_BFilter.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47656889._000048.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "Zee_CFBV";
        s.type = "BKG_Z";
        s.sigma_pb    = 2336.1;
        s.lumi_pb     = 3549;
        s.filter_eff  = 1.301249e-1;
        s.k_factor    = 0.9332;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.700787.Sh_2214_Zee_maxHTpTV2_CFilterBVeto.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47658627._000207.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "ttbar_llqq";
        s.type = "BKG_T";
        s.sigma_pb    = 811.29;
        s.lumi_pb     = 3549;
        s.filter_eff  = 4.384566e-1;
        s.k_factor    = 1.138433852;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47659319._000078.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "ttbar_llll";
        s.type = "BKG_T";
        s.sigma_pb    = 85.482;
        s.lumi_pb     = 3549;
        s.filter_eff  = 1;
        s.k_factor    = 1.138433852;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47762465._000196.pool.root.1"
        };
        samples.push_back(s);
    }
    {
        Sample s;
        s.name = "ttbar_qqqq";
        s.type = "BKG_T";
        s.sigma_pb    = 811.29;
        s.lumi_pb     = 3549;
        s.filter_eff  = 4.561727e-1;
        s.k_factor    = 1.1384339;
        s.files = {
            "/srv/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.*.pool.root.1"
            // "/srv/SMbackgrounds/mc23_13p6TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.deriv.DAOD_LLP1.e8514_s4162_r15540_p7077/DAOD_LLP1.47657677._000301.pool.root.1"
        };
        samples.push_back(s);
    }

    // Histograms per sample
    std::map<std::string, std::map<std::string, TH1F*>> histos;
    std::map<std::string, double> super_factors;
    // calculate Progressive bins
    static std::vector<double> lep1_bins = make_log_bins(lep_pt_min, 3950.0, 20); //3950
    static std::vector<double> lep2_bins = make_log_bins(lep_pt_min, 2950.0, 20); //2950
    static std::vector<double> jet_bins  = make_log_bins(jet_pt_min, 3350.0, 20); //3350
    static std::vector<double> disp_bins = make_log_bins(0.01, 100.0, 20);   //0.5, 250
    static std::vector<double> m_lljj_bins = make_log_bins(95, 6950, 30);  
    
    for (auto &s : samples) {
        const std::string &sname = s.name;
        const char *y_title_latex = "event weight: #sigma #scale[0.7]{#int} #font[12]{L} #varepsilon k w^{i} / W_{tot}";   //#scale[0.7]{#int #font[12]{L}}
        const char *y_title_latex_ljj = "events / 23 GeV";
        const char *y_title_latex_lljj = "events / 230 GeV";
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
            20, 0.0, 5.5
        );
        histos[sname]["dr_lep1_jj"] = new TH1F(
            ("h_"+sname+"_dr_lep1_jj").c_str(),
            (std::string(" ;#DeltaR(#it{l_{1}, jj} );") + y_title_latex).c_str(),
            20, 0.0, 5.5
        );
        // int n_m_l2jj = m_l2jj_bins.size() - 1;
        histos[sname]["m_l2jj"] = new TH1F(
            ("h_"+sname+"_m_l2jj").c_str(),
            (std::string(" ;m(#it{l_{2}jj}) / GeV;") + y_title_latex_ljj).c_str(),
            30, 50, 750   //50, 50, 750
        );
        int n_m_lljj = m_lljj_bins.size() - 1;
        histos[sname]["m_lljj"] = new TH1F(
            ("h_"+sname+"_m_lljj").c_str(),
            (std::string(" ;m(#it{lljj}) / GeV;") + y_title_latex_lljj).c_str(),
            n_m_lljj, m_lljj_bins.data()
            // 50, 1050, 6550  // adjust range/bins to your taste /5950
        );
        histos[sname]["m_jj"] = new TH1F(
            ("h_"+sname+"_m_jj").c_str(),
            (std::string(" ;m(#it{jj}) / GeV;") + y_title_latex).c_str(),
            20, 50, 525  // adjust range/bins to your taste
        );
        if (s.type == "SIG") {
            int n_disp = disp_bins.size() - 1;
            histos[sname]["disp_reco_d0"] = new TH1F(
                ("h_"+sname+"_disp_reco_d0").c_str(),
                (std::string(" ;L_{reco}^{d0} / mm;") + y_title_latex).c_str(),
                n_disp, disp_bins.data()
            );
            histos[sname]["disp_truth_Lxyz"] = new TH1F(
                ("h_"+sname+"_disp_truth_Lxyz").c_str(),
                (std::string(" ;L_{truth}^{xyz} / mm;") + y_title_latex).c_str(),
                n_disp, disp_bins.data()
            );
        }
        TH1::SetDefaultSumw2(true);
    }

    // Process each sample
    for (auto &s : samples) {
        process_sample(s, histos, super_factors,
                       use_metadata_totalW ? &metadata : nullptr,
                       use_metadata_totalW);
    }

    // Compute pooled S, B, and Z
    std::vector<std::string> variables = {"lep1_pt", "lep2_pt", "disp", "jet_pt", "dilep_q", "dr_lep1_jj", "dr_lep2_jj", "m_l2jj", "m_lljj", "m_jj"};
    std::map<std::string, TH1F*> pooled_S;
    std::map<std::string, TH1F*> pooled_B;
    std::map<std::string, TH1F*> pooled_Z;
    for (auto &var : variables) {
        TH1F *sig_template = nullptr;
        TH1F *bkg_template = nullptr;
        for (auto &s : samples) {
            TH1F *h = histos[s.name][var];
            if (s.type.substr(0,3) == "SIG" && !sig_template) sig_template = h;
            if (s.type.substr(0,3) == "BKG" && !bkg_template) bkg_template = h;
        }
        if (!sig_template || !bkg_template) continue;
        TH1F *hS = (TH1F*)sig_template->Clone(("h_pooled_S_"+var).c_str());
        hS->Reset("ICES");
        TH1F *hB = (TH1F*)bkg_template->Clone(("h_pooled_B_"+var).c_str());
        hB->Reset("ICES");
        for (auto &s : samples) {
            TH1F *h = histos[s.name][var];
            if (s.type == "SIG") hS->Add(h);
            else if (s.type.substr(0,3) == "BKG") hB->Add(h);
            // data type is ignored
        }
        // Do not use rebin if progressive binning has been applied
        // hS->Rebin(2);
        // hB->Rebin(2);
        pooled_S[var] = hS;
        pooled_B[var] = hB;
        TH1F *hZ = (TH1F*)hS->Clone(("h_significance_"+var).c_str());
        hZ->Reset("ICES");
        hZ->SetTitle((std::string(";")+hS->GetXaxis()->GetTitle()+";Z [#sigma]").c_str());
        double b_min = 1e-3;  // spike killer sets bin to 0, but significance plot doesn't know this yet!
        int nbins = hS->GetNbinsX();
        for (int ibin = 1; ibin <= nbins; ++ibin) {
            double s = hS->GetBinContent(ibin);
            // Undo plotting buffs for significance calculation
            if (var == "m_lljj") {
                s /= 50; //2e5
            } else if (var == "m_l2jj") {
                s /= 50;  //2e5
            } else if (var == "m_jj") {
                s /= 1.0;
            }
            double b = hB->GetBinContent(ibin);
            double Z = 0.0;
            if (b >= b_min && s > 0) {
                double arg = 1.0 + s / b;
                if (arg > 0) {
                    double val = 2.0 * ((s+b)*std::log(arg) - s);
                    if (val > 0) Z = std::sqrt(val);
                }
            }
            hZ->SetBinContent(ibin, Z);
        }
        pooled_Z[var] = hZ;
    }

    // Build group histos (SIG, Zmumu, ttbar) from the individual samples
    std::map<std::string, std::vector<std::string>> group_definitions;
    group_definitions["SIG"] = {"WR5000_N50_2L2J"};
    group_definitions["Zmumu"] = {"Zmumu", "Zmumu_CVBV", "Zmumu_BF", "Zmumu_CFBV"};
    group_definitions["Zee"] = {"Zee_CVBV", "Zee_BF", "Zee_CFBV"};
    group_definitions["ttbar"] = {"ttbar_llqq", "ttbar_llll", "ttbar_qqqq"};
    group_definitions["DATA"] = {"Data_2022"}; 
    std::map<std::string, std::map<std::string, TH1F*>> group_histos;

    for (auto &var : variables) {
        for (auto &gd : group_definitions) {
            const std::string &gname = gd.first;
            const auto &members = gd.second;
            TH1F *template_h = nullptr;
            for (auto &s : samples) {
                if (std::find(members.begin(), members.end(), s.name) != members.end()) {
                    template_h = histos[s.name][var];
                    break;
                }
            }
            if (!template_h) continue;
            TH1F *h_group = (TH1F*)template_h->Clone(("h_"+gname+"_"+var).c_str());
            h_group->Reset("ICES");
            for (auto &s : samples) {
                if (std::find(members.begin(), members.end(), s.name) != members.end()) {
                    h_group->Add(histos[s.name][var]);
                }
            }
            group_histos[var][gname] = h_group;
        }
    }

    // Apply Histogram smooting ("spike killer" depricated) for BKG
    std::map<std::string, double> n_reference_events;
    n_reference_events["SIG"]   = 0.0;
    n_reference_events["Zmumu"] = 1e3;
    n_reference_events["Zee"] = 1e3;
    n_reference_events["ttbar"] = 1e3;
    for (auto &var : variables) {
        for (auto &gh : group_histos[var]) {
            const std::string &gname = gh.first;
            TH1F *h = gh.second;
            // if (gname != "SIG" && var != "disp") {  // && var != "disp"
            //     spike_killer(h, n_reference_events[gname], 2.0);
            // }
            // if (gname != "DATA") {   //&& gname != "SIG"
            //     if (h->GetNbinsX() > 10) h->Smooth(1);
            // }
        }
    }

    // Make canvases
    std::cout << "Plotting results to " << output_dir << std::endl;
    make_sig_canvas("lep1_pt", "Leading lepton p_{T}",
                    output_dir + "/leading_lepton_pt_with_significance.png",
                    true, false, group_histos["lep1_pt"], pooled_Z["lep1_pt"]);
    make_sig_canvas("lep2_pt", "Subleading lepton p_{T}",
                    output_dir + "/subleading_lepton_pt_with_significance.png",
                    true, false, group_histos["lep2_pt"], pooled_Z["lep2_pt"]);
    make_sig_canvas("disp", "SV / HNL displacement",
                    output_dir + "/sv_displacement_with_significance.png",
                    true, false, group_histos["disp"], pooled_Z["disp"]);
    make_sig_canvas("jet_pt", "Di-jet p_{T}",
                    output_dir + "/dijet_pt_with_significance.png",
                    true, false, group_histos["jet_pt"], pooled_Z["jet_pt"]);
    make_sig_canvas("dilep_q", "Dilepton charge (OS / SS)",
                    output_dir + "/dilepton_charge_with_significance.png",
                    true, true, group_histos["dilep_q"], pooled_Z["dilep_q"]);
    make_sig_canvas("dr_lep1_jj", "#DeltaR(l_{1}, jj)",
                    output_dir + "/dr_lep1_jj_with_significance.png",
                    true, false, group_histos["dr_lep1_jj"], pooled_Z["dr_lep1_jj"]);
    make_sig_canvas("dr_lep2_jj", "#DeltaR(l_{2}, jj)",
                    output_dir + "/dr_lep2_jj_with_significance.png",
                    true, false, group_histos["dr_lep2_jj"], pooled_Z["dr_lep2_jj"]);
    make_sig_canvas("m_l2jj", "Invariant mass m(l_{2}jj)",
                    output_dir + "/m_l2jj_with_significance.png",
                    false, false, group_histos["m_l2jj"], pooled_Z["m_l2jj"]);
    make_sig_canvas("m_lljj", "Invariant mass m(l_{1}l_{2}jj)",
                    output_dir + "/m_lljj_with_significance.png",
                    false, false, group_histos["m_lljj"], pooled_Z["m_lljj"]);
    make_sig_canvas("m_jj", "Invariant mass m(jj)",
                    output_dir + "/m_jj_with_significance.png",
                    true, false, group_histos["m_jj"], pooled_Z["m_jj"]);

    // DEBUG: plot only DR(l2, jj) for signal, linear scale ===
    {
        for (const auto& gname : {
            std::string("Zmumu"),
            std::string("Zee"),
            std::string("ttbar"),
            std::string("SIG")
        })
        {
            auto& varMap = group_histos["m_lljj"]; // or the map appropriate for var
            if (varMap.count(gname) == 0) continue;

            auto h_lljj = varMap[gname];
            auto h_l2jj = group_histos["m_l2jj"][gname];

            std::cout
                << gname
                << "\n  m_lljj : " << h_lljj->GetEntries()
                << " entries"
                << "\n  m_l2jj : " << h_l2jj->GetEntries()
                << " entries"
                << "\n";
        }
        auto itS1 = pooled_S.find("dr_lep1_jj");
        auto itS2 = pooled_S.find("dr_lep2_jj");
        if (itS1 != pooled_S.end() && itS2 != pooled_S.end()) {
            TH1F *hL1 = (TH1F*)itS1->second->Clone("h_dr_lep1_jj_debug");
            TH1F *hL2 = (TH1F*)itS2->second->Clone("h_dr_lep2_jj_debug");
            hL1->SetTitle(";#DeltaR(#it{l, jj});event weight"); // common title
            // Style L1
            hL1->SetLineColor(kRed+1);
            hL1->SetLineWidth(2);
            hL1->SetFillStyle(0);
            hL1->SetMarkerStyle(20);
            hL1->SetMarkerColor(kRed+1);
            hL1->SetMarkerSize(0.8);
            // Style L2
            hL2->SetLineColor(kBlue+1);
            hL2->SetLineWidth(2);
            hL2->SetFillStyle(0);
            hL2->SetMarkerStyle(24);
            hL2->SetMarkerColor(kBlue+1);
            hL2->SetMarkerSize(0.8);
            // Linear y, linear x
            TCanvas *cDebug = new TCanvas("c_dr_lep1_lep2_jj_debug","Debug: DR(l1,jj) and DR(l2,jj), signal only", 800, 600);
            cDebug->cd();
            gPad->SetLogy(kFALSE);
            gPad->SetLogx(kFALSE);
            // Axis cosmetics
            hL1->GetXaxis()->SetTitleSize(0.045);
            hL1->GetXaxis()->SetLabelSize(0.04);
            hL1->GetYaxis()->SetTitleSize(0.045);
            hL1->GetYaxis()->SetLabelSize(0.04);
            // Set y‑range from both histos
            double maxY = std::max(hL1->GetMaximum(), hL2->GetMaximum());
            hL1->SetMinimum(0.0);
            hL1->SetMaximum(1.2 * maxY);
            // Draw
            hL1->Draw("HIST E1");
            hL2->Draw("HIST E1 SAME");
            // Legend
            TLegend *leg = new TLegend(0.60, 0.70, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->AddEntry(hL1, "#DeltaR(#it{l_{1}, jj})", "l");
            leg->AddEntry(hL2, "#DeltaR(#it{l_{2}, jj})", "l");
            leg->Draw();
            cDebug->SaveAs((output_dir + "/debug_dr_lep1_lep2_jj_signal_linear.png").c_str());
        } else {
            std::cerr << "DEBUG: pooled_S[\"dr_lep1_jj\"] or pooled_S[\"dr_lep2_jj\"] not found.\n";
        }
    }

    std::cout << "\nAnalysis successfully completed. Results saved to " << output_dir << std::endl;
}
