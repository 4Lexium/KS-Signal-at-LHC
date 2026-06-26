#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPolyLine.h"
#include "TColor.h"
#include "TMath.h"
#include "TH1F.h"
#include <iostream>

double HNLWidth(double mN_GeV, double mWR_GeV) {
    // HNL Simplified Decay Width (scale factor) = 1e-5 * mN^5 / mWR^4 (in GeV)
    return 1e-5 * TMath::Power(mN_GeV, 5) / TMath::Power(mWR_GeV, 4);
}

double HNLLifetime_mm(double mN_GeV, double mWR_GeV) {
    // DECAY_LT = 1.973e-13 / Gamma   (c*tau in mm, using hbar*c ≈ 1.973e-13 GeV·mm)
    double gamma = HNLWidth(mN_GeV, mWR_GeV);
    return 1.973e-13 / gamma;
}

void plot_nominal_decaylength() {
    std::cout << "User Test 008:" << std::endl;
    gStyle->SetOptStat(0);

    // Mass range and sampling for HNL
    const int N = 500;
    const double mN_min = 5;    // HNL mass calculation range (GeV)
    const double mN_max = 1000;  // 1 TeV

    // Two M_WR values in GeV for a range bar
    const double mWR1 = 3000.0;    // 3 TeV
    const double mWR2 = 50000.0;   // 50 TeV

    // Create canvas first
    auto c = new TCanvas("c","HNL lifetime vs mass",1200,800);
    c->SetLogx();
    c->SetLogy();
    c->SetGridx();
    c->SetBottomMargin(0.25);
    // c->SetGridy();

    // Choose a nicer x range for the *frame* so log spacing by decade looks regular
    // and PDG points at a few GeV are visible:
    double frame_xmin = 2.0;    // start at 1 GeV
    double frame_xmax = 1000.0; // up to 1 TeV

    // Y range in mm
    double y_min = 1e-17;
    double y_max = 5;

    TH1F *frame = c->DrawFrame(frame_xmin, y_min, frame_xmax, y_max);
    frame->SetTitleFont(42, "t");
    frame->SetTitle("HNL Nominal Displacement");
    frame->SetTitleOffset(0.8, "t");  // "t" = main titl
    frame->GetXaxis()->SetTitle("mass [GeV]");
    frame->GetXaxis()->SetTitleOffset(1.2); 
    frame->GetYaxis()->SetTitleOffset(1.3); 
    frame->GetYaxis()->SetTitle("nominal displacement c#tau [mm]");
    frame->GetXaxis()->SetTitleFont(32);
    frame->GetYaxis()->SetTitleFont(32);

    // Build HNL curves
    auto g1 = new TGraph(N);
    auto g2 = new TGraph(N);

    for (int i = 0; i < N; ++i) {
        // try logarithmic spacing if you want a smooth band on log x:
        double x = mN_min * TMath::Power(mN_max/mN_min, double(i)/(N-1)); // log spacing
        // or keep linear spacing (this does NOT affect axis scaling, only point density)
        // double x = mN_min + (mN_max - mN_min) * double(i)/(N-1);

        double y1 = HNLLifetime_mm(x, mWR1);
        double y2 = HNLLifetime_mm(x, mWR2);

        g1->SetPoint(i, x, y1);
        g2->SetPoint(i, x, y2);
    }

    g1->SetLineColor(kRed+1);
    g1->SetLineWidth(2);
    g1->SetLineStyle(2);

    g2->SetLineColor(kBlue+1);
    g2->SetLineWidth(2);
    g2->SetLineStyle(2);

    // Build a polygon to shade area between g1 and g2
    auto band = new TGraph(2*N);
    for (int i = 0; i < N; ++i) {
        double x, y;
        g1->GetPoint(i, x, y);
        band->SetPoint(i, x, y);
    }
    for (int i = 0; i < N; ++i) {
        double x, y;
        g2->GetPoint(N-1-i, x, y);
        band->SetPoint(N + i, x, y);
    }

    band->SetFillColorAlpha(kAzure+1, 0.3);
    band->SetLineColorAlpha(kAzure+1, 0.0);

    // Draw band and curves
    band->Draw("f same");
    g1->Draw("l same");
    g2->Draw("l same");

    /*
        PDG Data (https://pdg.lbl.gov/2024/listings/particle_properties.html)
        Higgs:
            h: 125.20 GeV, 4.8*10^-14 mm
            HH: 250 GeV, 4.8*10^-14 mm
            H++: 100 GeV, 10^-4 mm
        Mesons:
            J/psi(cc): 3.097 GeV, 2.2*10^-12 mm
            Upsilon(bb): 9.46 GeV, 3.6*10^-12 mm
            B(bd): 5.28 GeV, 0.45 mm
            T(tt): 346 GeV, 1.5*10^-14 mm
        Vector gauges:
            Z: 91.19 GeV, 9*10^-17 mm 
        Diverse:
            Axions: 

    */

    const int NP = 8;
    double xP[NP] = {
        125.2,   // h
        // 250.0,   // HH
        // 100.0,   // H++
        3.097,   // J/psi
        9.46,    // Upsilon
        5.28,    // B
        346.0,   // t tbar
        91.19    // Z
    };

    double yP[NP] = {
        4.8e-14,   // h
        // 4.8e-14,   // HH
        // 1.0e-4,    // H++
        2.2e-12,   // J/psi
        3.6e-12,   // Upsilon
        4.5e-1,    // B
        1.5e-16,   // t tbar
        9.0e-17    // Z
    };

    auto gPDG = new TGraph(NP, xP, yP);
    gPDG->SetMarkerStyle(20);
    gPDG->SetMarkerSize(1.3);
    gPDG->SetMarkerColor(kBlack);
    gPDG->Draw("P same");

    // Labels for PDG points
    auto latex = new TLatex();
    latex->SetTextSize(0.036);
    latex->SetTextFont(32);
    latex->SetTextColor(kBlack);

    latex->DrawLatex(xP[0]*1.05, yP[0]*2.0,  "h");
    // latex->DrawLatex(xP[1]*1.05, yP[1]*2.0,  "HH");
    // latex->DrawLatex(xP[2]*0.72, yP[2]*0.45,  "H^{++}");
    latex->DrawLatex(xP[1]*1.1,  yP[1]*2.0,  "J/#psi(c#bar{c})");
    latex->DrawLatex(xP[2]*1.1,  yP[2]*2.0,  "#Upsilon(b#bar{b})");
    latex->DrawLatex(xP[3]*1.1, yP[3]*0.5,  "B^{0}(#bar{b}d)");
    latex->DrawLatex(xP[4]*1.0, yP[4]*2.0,  "T(t#bar{t})");
    latex->DrawLatex(xP[5]*1.05, yP[5]*2.0,  "Z");


    // --- Rotated label inside the shaded band (draw AFTER everything) ---
    // Use log-log slope of the mWR1 curve
    double x1 = 100;
    double x2 = 500;
    double y1 = HNLLifetime_mm(x1, mWR1);
    double y2 = HNLLifetime_mm(x2, mWR1);

    double angle = -21.0;

    double x_mid = 200.0;
    double y_mid = TMath::Sqrt(
        HNLLifetime_mm(x_mid, mWR1) *
        HNLLifetime_mm(x_mid, mWR2)
    );

    const int NB2 = 3;

    double xB2[NB2] = {
        750.0,   // H++
        500.0,   // Heavy H (HH)
        // 10.0,    // ALP, intermediate lifetime
        100.0    // Neutralino
    };

    double yB2[NB2] = {
        1e-10,    // H++
        1e-14,   // HH (suppressed / LLP-like scenario)
        // 1e0,     // ALP (intermediate lifetime)
        1e-3      // Neutralino (LLP-scale example)
    };

    // Asymmetric errors in x (mass) and y (ctau)

    // Lower errors in x (x - exlB)
    double exlB2[NB2] = {
        250.0,    // H++: 500 - 50 = 450 GeV
        250.0,    // HH : 500 - 50 = 450 GeV
        // 2.0,     // ALP mid:        10 - 2 = 8 GeV
        50.0     // Neutralino:     100 - 20 = 80 GeV
    };

    // Upper errors in x (x + exhB)
    double exhB2[NB2] = {
        250.0,   // H++: 500 + 100 = 600 GeV
        500.0,   // HH : 500 + 150 = 650 GeV
        // 5.0,     // ALP mid:        10 + 5 = 15 GeV
        900.0     // Neutralino:     100 + 50 = 150 GeV
    };

    // Lower errors in y (y - eylB)
    double eylB2[NB2] = {
        9.9e-11,    // H++:    1e-4   - 5e-5  = 5e-5 m
        9.9e-15,   // HH:     1e-10  - 5e-11 = 5e-11 m
        // 5e-1,    // ALP mid:1e0    - 5e-1  = 5e-1 m
        9.9e-4     // Neutralino: 1 - 0.2 = 0.8 m
    };

    // Upper errors in y (y + eyhB)
    double eyhB2[NB2] = {
        9.9e-9,    // H++:    1e-4   + 5e-4  = 6e-4 m
        9.9e-13,   // HH:     1e-10  + 5e-10 = 6e-10 m
        // 5e0,     // ALP mid:1e0    + 5e0   = 6 m
        9.9e-2      // Neutralino: 1 + 10 = 11 m
    };

    auto gBSM = new TGraphAsymmErrors(
        NB2,
        xB2,  yB2,
        exlB2, exhB2,   // x low / x high errors
        eylB2, eyhB2    // y low / y high errors
    );

    gBSM->SetMarkerStyle(21);
    gBSM->SetMarkerSize(1.4);
    gBSM->SetMarkerColor(kMagenta+2);
    gBSM->SetLineColor(kMagenta+2);
    gBSM->SetLineWidth(2);
    gBSM->Draw("P SAME");

    auto latex2 = new TLatex();

    latex2->SetTextSize(0.036);
    latex2->SetTextFont(32);
    latex2->SetTextColor(kMagenta+2);
    latex2->DrawLatex(xB2[0]*0.75, yB2[0]*1.6, "H^{++}");
    latex2->DrawLatex(xB2[1]*1.1, yB2[1]*2.0, "HH");
    // latex2->DrawLatex(xB2[2]*1.1, yB2[2]*1.5, "ALP");
    latex2->DrawLatex(xB2[2]*0.82, yB2[2]*4, "#tilde{#chi}^{0}_{1}");



//-------------------------------------------------------------------------------------

    const int NB3 = 6;

    double xB3[NB3] = {
        50.0,  
        50.0,   
        100.0,   
        100.0,
        500.0,
        500.0,    
    };

    double yB3[NB3] = {
        0.127,    
        0.0456,   
        5.88e-4,
        3.17e-5,
        8.93e-8,
        1.27e-8,
    };

    // Asymmetric errors in x (mass) and y (ctau)

    // Lower errors in x (x - exlB)
    double exlB3[NB3] = {
        0,    
        0,    
        0,
        0, 
        0,
        0, 
    };

    // Upper errors in x (x + exhB)
    double exhB3[NB3] = {
        0, 
        0,   
        0,
        0,
        0,
        0,
    };

    // Lower errors in y (y - eylB)
    double eylB3[NB3] = {
        0.067,   
        0.024, 
        3.16e-4,  
        1.71e-5,
        5.29e-8,
        0.76e-8,
    };

    // Upper errors in y (y + eyhB)
    double eyhB3[NB3] = {
        0.187,   
        0.0672, 
        8.6e-4,  
        4.63e-5,
        12.57e-8,
        1.78e-8,
    };

    auto gSIGerr = new TGraphAsymmErrors(
        NB3,
        xB3,  yB3,
        exlB3, exhB3,   // x low / x high errors
        eylB3, eyhB3    // y low / y high errors
    );

    gSIGerr->SetMarkerStyle(21);
    gSIGerr->SetMarkerSize(1.4);
    gSIGerr->SetMarkerColor(kGreen+3);
    gSIGerr->SetLineColor(kGreen+3);
    gSIGerr->SetLineWidth(2);
    gSIGerr->Draw("P SAME");

    auto latex3 = new TLatex();

    latex3->SetTextSize(0.036);
    latex3->SetTextFont(32);
    latex3->SetTextColor(kGreen+3);
    latex3->DrawLatex(xB3[0]*0.75, yB3[0]*1.8, "S1");
    latex3->DrawLatex(xB3[1]*0.75, yB3[1]*0.7, "S1'");
    latex3->DrawLatex(xB3[2]*0.75, yB3[2]*0.08, "S2");
    latex3->DrawLatex(xB3[3]*0.75, yB3[3]*0.15, "S2'");
    latex3->DrawLatex(xB3[4]*0.75, yB3[4]*0.15, "S3");
    latex3->DrawLatex(xB3[5]*0.75, yB3[5]*0.1, "S3'");


//-------------------------------------------------------------------------------------






    // Legend
    auto leg = new TLegend(0.08, 0.02, 0.75, 0.22);
    leg->SetTextFont(32);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->SetNColumns(2);
    leg->AddEntry(band, "HNL Displaced Region (unboosted)", "f");
    leg->AddEntry(gPDG, "SM PDG values (2024)", "p");
    leg->AddEntry((TObject*)nullptr, "c#tau_{N} #propto 10^{5} #times #frac{m_{W_{R}}^{4}}{m_{N}^{5}}", "");
    leg->AddEntry(gBSM, "BSM LLP benchmarks*", "lep");
    leg->AddEntry(gSIGerr, "SIG Samples** (boosted)", "lep");
    leg->Draw();

    auto bandLabel = new TLatex();
    bandLabel->SetTextSize(0.035);
    bandLabel->SetTextFont(32);
    bandLabel->SetTextColor(kAzure+3);
    bandLabel->SetTextAngle(angle);
    bandLabel->SetTextAlign(22); // center on (x_mid, y_mid)
    c->cd();                     // make sure we are on the right pad
    bandLabel->DrawLatex(x_mid, y_mid, "HNL parameter space");

    double x_W2 = 200.0;
    double x_W1 = 350.0;
    double y_W1 = HNLLifetime_mm(x_W1, mWR1)*TMath::Power(10.0, +0.6);  // red line
    double y_W2 = HNLLifetime_mm(x_W2, mWR2)*TMath::Power(10.0, -0.6);  // blue line

    auto band_W1 = new TLatex();
    band_W1->SetTextSize(0.03);
    band_W1->SetTextFont(32);
    band_W1->SetTextColor(kRed+1);
    band_W1->SetTextAngle(angle);
    band_W1->SetTextAlign(22); // center on (x_mid, y_mid)
    c->cd();                     // make sure we are on the right pad
    band_W1->DrawLatex(x_W1, y_W1, "m_{W_{R}} = 3 TeV");
    auto band_W2 = new TLatex();
    band_W2->SetTextSize(0.03);
    band_W2->SetTextFont(32);
    band_W2->SetTextColor(kBlue+1);
    band_W2->SetTextAngle(angle);
    band_W2->SetTextAlign(22); // center on (x_mid, y_mid)
    c->cd();                     // make sure we are on the right pad
    band_W2->DrawLatex(x_W2, y_W2, "m_{W_{R}} = 50 TeV");

    gPad->Update();
    TPaveText *pt = (TPaveText*)gPad->GetPrimitive("title");
    if (pt) {
        pt->SetTextSize(0.045);   
        pt->SetTextFont(32);
    }
    c->SetTicks(1,1);
    c->RedrawAxis();
    c->SaveAs("/home/user/WORKSPACE/Manual_LLP_analysis/Figures/hnl_nominal_decL.png");
    std::cout << "Nominal Decay Lenght Plot Generated Successfully!" << std::endl;
}
