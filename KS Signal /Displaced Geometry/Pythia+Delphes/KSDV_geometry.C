#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h" 
#include "TCanvas.h"
#include "TSystem.h"
#include "TEllipse.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

struct Sample {
    TString path;
    TString label;
};

/*
Before run in terminal:
source /storage/user/madgraph/MG5_aMC_v2_9_25/Delphes/DelphesEnv.sh
root -l
gSystem->Load("libDelphes");
.L KS_Detector_analysis.C or .L KS_Detector_analysis.C+
KS_Detector_analysis()
*/

void KS_Detector_scan_Geometry() {
    cout << "User Test 031: Vertex Displacement Scan" << endl;
    
    // Load Delphes library
    gSystem->Load("libDelphes");
    
    // Prepare files
    //--------------------------------------------------------------------------------------------------
    TString baseFigPath = "/home/user/WORKSPACE/";
    TString outDir = "/home/user/WORKSPACE/KS_Delphes_analysis/";
    gSystem->mkdir(baseFigPath + outDir, kTRUE);

    vector<Sample> samples = {
        {"/storage/user/WORKBENCH/runs/KS_delphes_decay_PythGEOGUSSE_MN50_MWR5000/Events/run_01/tag_1_delphes_events.root", "MN50_MWR5000"},
    };
    vector<TH2F*> hists;
    vector<TH1F*> histsLxy;
    vector<TH1F*> histsLz;
    // File Loop
    //--------------------------------------------------------------------------------------------------
    for (auto &s : samples) {
        cout << "Processing: " << s.label << endl;

        TFile *f = TFile::Open(s.path);
        TTree *t = (TTree*)f->Get("Delphes");
        if (!f || f->IsZombie()) {cout << "Cannot open " << s.path << endl; continue;}
        if (!t) {cout << "No Delphes tree in " << s.label << endl; continue;}
        ofstream ofs("/home/user/WORKSPACE/KS_detector_analysis/BRANCH_INFO.txt", ios::trunc);
        ofs.close();
        gSystem->RedirectOutput("/home/user/WORKSPACE/KS_detector_analysis/BRANCH_INFO.txt");
        t->Print("*");
        gSystem->RedirectOutput(0);

        // Set branch addresses
        //--------------------------------------------------------------------------------------------------
        
        // Particle branches
        Int_t particle_size;
        Int_t PID[10000], Status[10000], M1[10000], M2[10000];
        Float_t X[10000], Y[10000], Z[10000], T_time[10000];
        t->SetBranchAddress("Particle_size", &particle_size);
        t->SetBranchAddress("Particle.PID", PID);
        t->SetBranchAddress("Particle.Status", Status);
        t->SetBranchAddress("Particle.M1", M1);
        t->SetBranchAddress("Particle.M2", M2);
        t->SetBranchAddress("Particle.X", X);
        t->SetBranchAddress("Particle.Y", Y);
        t->SetBranchAddress("Particle.Z", Z);
        t->SetBranchAddress("Particle.T", T_time);

        // Track branches (for decay vertex)
        Int_t track_size, TrackPID[10000]; // TrackM1[10000], TrackM2[10000];
        Float_t TrackXd[10000], TrackYd[10000], TrackZd[10000];
        t->SetBranchAddress("Track_size", &track_size);
        t->SetBranchAddress("Track.PID", TrackPID);
        // t->SetBranchAddress("Track.M1", TrackM1);
        // t->SetBranchAddress("Track.M2", TrackM2);
        t->SetBranchAddress("Track.Xd", TrackXd);
        t->SetBranchAddress("Track.Yd", TrackYd);
        t->SetBranchAddress("Track.Zd", TrackZd);


        // Crate Histograms per sample
        //--------------------------------------------------------------------------------------------------
        TH1F *h_Lz = new TH1F("h_Lz","HNL decay length;L_{z} [nm];Events",100,0,10);
        TH1F *h_Lxy = new TH1F("h_Lxy","HNL transverse decay length;L_{xy} [nm];Events",100,0,10);
        TH2F *h = new TH2F(
            "h_" + s.label,
            s.label + ";X [mm];Y [mm]",
            200, -150, 150,   //250
            200, -150, 150
        );
        h->SetDirectory(0);
        h_Lxy->SetDirectory(0);
        h_Lz->SetDirectory(0);
        Long64_t nentries = t->GetEntries();
        int nDV = 0;

        // Loop over events
        //--------------------------------------------------------------------------------------------------
        for (Long64_t i = 0; i < nentries; i++) {
            t->GetEntry(i);
            for (int j = 0; j < particle_size; j++) {

                // if HNL determine its daughters
                if (abs(PID[j]) != 9900012) continue;
                bool hasHNLdaughter = false;
                for (int k=0; k<particle_size; k++) {
                    if (M1[k]==j && abs(PID[k])==9900012) {
                        hasHNLdaughter = true;
                        break;
                    }
                }
                if (hasHNLdaughter) continue;  

                // Find final-state daughters of HNL with status 1
                for (int k = 0; k < particle_size; k++) {
                    if (M1[k] != j && M2[k] != j) continue;
                    // float dx = TrackXd[k];
                    // float dy = TrackYd[k];
                    // float dz = TrackZd[k];
                    float dx = X[k]*1000;
                    float dy = Y[k]*1000;
                    float dz = Z[k]*1000;
                    float Lxy = sqrt(dx*dx + dy*dy);
                    h_Lxy->Fill(Lxy);
                    h_Lz->Fill(dz);
                    h->Fill(dx, dy);
                    nDV++;
                    break; // one DV per HNL
                }
            }
        }

        hists.push_back(h);
        histsLxy.push_back(h_Lxy);
        histsLz.push_back(h_Lz);
        f->Close();
    }
    // cout << "MOOOOOOOOOOOOOOOOOOOOO" << endl;
    // Visualization of Plots
    //--------------------------------------------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);

    int n = hists.size();
    int cols = ceil(sqrt(n));
    int rows = ceil((float)n / cols);

    // Detector outer radius (mm)
    double R_detector = 125.0;

    // Zoom region (mm)
    double zoom_range = 10;

    for (size_t i = 0; i < hists.size(); i++) {
        
        // Decode file name 
        TString label = samples[i].label;
        TString physicsTitle;
        if (label.Contains("MN50") && label.Contains("3000"))
            physicsTitle = "DV(N): m_{N} = 50 GeV, m_{W_{R}} = 3 TeV";
        if (label.Contains("MN100") && label.Contains("50000"))
            physicsTitle = "DV(N): m_{N} = 100 GeV, m_{W_{R}} = 50 TeV";
        if (label.Contains("MN100") && label.Contains("10000"))
            physicsTitle = "DV(N): m_{N} = 100 GeV, m_{W_{R}} = 10 TeV";
        if (label.Contains("MN50") && label.Contains("50000"))
            physicsTitle = "DV(N): m_{N} = 50 GeV, m_{W_{R}} = 50 TeV";
        else if (label.Contains("MN500") && label.Contains("5000"))
            physicsTitle = "DV(N): m_{N} = 500 GeV, m_{W_{R}} = 5 TeV";
        else if (label.Contains("MN500") && label.Contains("3000"))
            physicsTitle = "DV(N): m_{N} = 500 GeV, m_{W_{R}} = 3 TeV";
        else if (label.Contains("MN500") && label.Contains("50000"))
            physicsTitle = "DV(N): m_{N} = 500 GeV, m_{W_{R}} = 50 TeV";
        else if (label.Contains("MN500") && label.Contains("10000"))
            physicsTitle = "DV(N): m_{N} = 500 GeV, m_{W_{R}} = 10 TeV";
        else if (label.Contains("MN50") && label.Contains("5000"))
            physicsTitle = "DV(N): m_{N} = 50 GeV, m_{W_{R}} = 5 TeV";
        else
            physicsTitle = label;

        TCanvas *cXY = new TCanvas("cXY_"+label,"L XY",800,600);
        cXY->cd();
        histsLxy[i]->Draw();
        cXY->SaveAs(baseFigPath + outDir + "HNL L_XY" + label + ".png");

        TCanvas *cZ = new TCanvas("cZ_"+label,"L >",800,600);
        cZ->cd();
        histsLz[i]->Draw();
        cZ->SaveAs(baseFigPath + outDir + "HNL L_Z" + label + ".png");

        TCanvas *c = new TCanvas("c_"+label,"DV XY "+label,1000, 900);
        c->cd();
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.08);
        gPad->SetLogz();

        hists[i]->SetTitle(physicsTitle);
        hists[i]->GetXaxis()->SetTitle("x_{DV} / mm");
        hists[i]->GetYaxis()->SetTitle("y_{DV} /mm");
        hists[i]->GetXaxis()->SetTitleSize(0.045);
        hists[i]->GetYaxis()->SetTitleSize(0.045);
        hists[i]->GetXaxis()->SetLabelSize(0.04);
        hists[i]->GetYaxis()->SetLabelSize(0.04);
        // hists[i]->GetXaxis()->SetLabelFont(32);
        // hists[i]->GetYaxis()->SetLabelFont(32);
        // hists[i]->GetZaxis()->SetLabelFont(32);
        hists[i]->GetXaxis()->SetTitleFont(32);
        hists[i]->GetYaxis()->SetTitleFont(32);
        hists[i]->GetZaxis()->SetTitleFont(32);
        hists[i]->SetTitleFont(32); 

        hists[i]->Draw("COLZ");
        gPad->SetFixedAspectRatio();

        // Draw detector boundary
        TEllipse *detectorCircle = new TEllipse(0, 0, R_detector, R_detector);
        detectorCircle->SetFillStyle(0);
        detectorCircle->SetLineColor(kRed+1);
        detectorCircle->SetLineWidth(4);
        detectorCircle->Draw("same");

        TLegend *circleLeg = new TLegend(0.4, 0.2, 0.55, 0.35);
        circleLeg->SetBorderSize(0);
        circleLeg->SetFillStyle(0);
        circleLeg->SetTextFont(32);
        circleLeg->SetTextSize(0.02);
        circleLeg->AddEntry(detectorCircle,"Outer Pixel Tracker (R=125 mm)","l");
        // circleLeg->Draw();

        // ZOOM-IN INSET PAD
        TPad *zoomPad = new TPad("zoomPad","zoom",0.55, 0.585, 0.85, 0.915);
        zoomPad->SetFillColor(0);
        zoomPad->SetFrameLineWidth(2);
        zoomPad->Draw();
        zoomPad->cd();
        zoomPad->SetLeftMargin(0.04);
        zoomPad->SetRightMargin(0.02);
        zoomPad->SetBottomMargin(0.04);
        zoomPad->SetTopMargin(0.02);

        // gPad->SetRightMargin(0.15);
        // gPad->SetLeftMargin(0.18);
        // gPad->SetBottomMargin(0.18);
        // gPad->SetTopMargin(0.08);

        TH2F *hZoom = (TH2F*)hists[i]->Clone("hZoom_"+label);
        hZoom->GetXaxis()->SetRangeUser(-zoom_range, zoom_range);
        hZoom->GetYaxis()->SetRangeUser(-zoom_range, zoom_range);

        hZoom->SetTitle("");
        hZoom->GetXaxis()->SetTitle("");
        hZoom->GetYaxis()->SetTitle("");
        hZoom->GetZaxis()->SetTitle("");
        // hZoom->GetXaxis()->SetTitleSize(0.06);
        // hZoom->GetYaxis()->SetTitleSize(0.06);
        // hZoom->GetXaxis()->SetLabelSize(0.05);
        // hZoom->GetYaxis()->SetLabelSize(0.05);

        hZoom->Draw("COL");   // COOLZ if with colorbar
        gPad->SetFixedAspectRatio();

        TEllipse *zoomCircle = new TEllipse(0, 0, zoom_range, zoom_range);
        zoomCircle->SetFillStyle(0);
        zoomCircle->SetLineColor(kPink+1);
        zoomCircle->SetLineWidth(3);
        zoomCircle->Draw("same");
        circleLeg->AddEntry(zoomCircle, TString::Format("DV cut (R=%.1f mm)", zoom_range), "l");
        c->cd();
        circleLeg->Draw();
        // Save for each sample file
        // c->cd();
        gPad->Update();
        TPaveText *pt = (TPaveText*)gPad->GetPrimitive("title");
        if (pt) {
            pt->SetTextSize(0.045);   // try 0.04–0.06
            pt->SetTextFont(32);
        }

        c->SaveAs(baseFigPath + outDir + "DV_XY_" + label + ".png");
    }

    cout<<"Plots Generated Sucessfully! :)=<=<"<<endl;
}

    