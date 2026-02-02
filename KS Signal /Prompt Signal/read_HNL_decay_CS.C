#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TStyle.h>

#include <zlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>

/*
Read total cross section from <init> to </init> block of an LHE.gz file.
Example extract:
<init>
2212 2212 7.000000e+03 7.000000e+03 0 0 247000 247000 -4 2
4.324637e+01 1.690500e-01 5.596562e+01 1
1.271925e+01 6.497499e-02 5.596562e+01 2
<generator name='MadGraph5_aMC@NLO' version='2.9.25'>please cite 1405.0301 </generator>
</init>
*/

double read_xsec_from_lhe_gz(const std::string &filename) {
    gzFile file = gzopen(filename.c_str(), "rb");
    if (!file) {
        std::cerr << "ERROR: cannot open " << filename << std::endl;
        return -1.0;
    }

    char buffer[4096];
    bool in_init = false;
    double sigma_tot = 0.0;

    while (gzgets(file, buffer, sizeof(buffer))) {
        std::string line(buffer);

        if (line.find("<init>") != std::string::npos) {
            in_init = true;
            continue;
        }

        if (!in_init)
            continue;

        if (line.find("</init>") != std::string::npos ||
            line.find("<generator") != std::string::npos) {
            break;
        }

        if (line.find("2212") == 0)
            continue;

        std::istringstream iss(line);
        double sigma, err, sigmax;
        int id;

        if (iss >> sigma >> err >> sigmax >> id) {
            sigma_tot += sigma;
        }
    }

    gzclose(file);
    return sigma_tot; // pb
}

void read_HNL_decay_CS() {
    std::cout << "User Test 004: " <<std::endl;
    gStyle->SetOptStat(0);
    std::cout << "Processing HNL decay (N->ljj) cross sections" << std::endl;
    std::string base = "/storage/alexanum/WORKBENCH/runs/";

    std::vector<double> tb_values = {0, 0.1, 0.5, 0.9};
    std::vector<int>    mN_values = {50, 60, 70, 80, 90, 100, 120, 150, 200, 250, 500};
    std::vector<int>    mWR_values = {3000, 5000};

    std::vector<int> colors = {kRed, kBlue, kGreen+2, kBlack, kOrange+1};
    std::vector<int> lineStyles = {1, 2, 3, 4, 5};

    // Setup Canvas, Info box and legends
    //-----------------------------------------------------------------------------------------------
    TCanvas *c = new TCanvas("c", "HNL decay cross sections", 1000, 800);
    c->SetLogy();
    c->SetGrid();

    TPaveText *infoBox = new TPaveText(0.12, 0.82, 0.5, 0.90, "NDC");
    infoBox->SetBorderSize(0);
    infoBox->SetFillColor(0);
    infoBox->SetFillStyle(0);
    infoBox->SetTextFont(32);
    infoBox->SetTextSize(0.04);
    infoBox->SetTextAlign(12);
    infoBox->AddText("pp#rightarrowW_{R}^{#pm}#rightarrowl^{#pm}(N#rightarrowljj)");
    // infoBox->AddText("l #in (e,#mu,#tau)");
    // infoBox->AddText("m_{N_{1}} #ll m_{N_{2,3}}");

    TLegend *colorLeg = new TLegend(0.12, 0.72, 0.45, 0.82);
    // colorLeg->SetHeader("tan#beta");
    colorLeg->SetBorderSize(0);
    colorLeg->SetFillColor(0);
    colorLeg->SetFillStyle(0);
    colorLeg->SetTextFont(32);
    colorLeg->SetNColumns(2);
    colorLeg->SetTextSize(0.04);

    TLegend *styleLeg = new TLegend(0.12, 0.68, 0.48, 0.72);
    // styleLeg->SetHeader("m_{W_{R}}");
    styleLeg->SetBorderSize(0);
    styleLeg->SetFillColor(0);
    styleLeg->SetFillStyle(0);
    styleLeg->SetTextFont(32);
    styleLeg->SetNColumns(2);
    styleLeg->SetTextSize(0.04);

    // Plotting:
    // tB (colour) >> mWR (linestyle) >> mN (xaxis)
    //-----------------------------------------------------------------------------------------------
    double ymin = 1e30;
    double ymax = 0.0;
    bool first = true;

    for (size_t tb_idx = 0; tb_idx < tb_values.size(); tb_idx++) {

        double tb = tb_values[tb_idx];
        int color = colors[tb_idx % colors.size()];

        for (size_t mwr_idx = 0; mwr_idx < mWR_values.size(); mwr_idx++) {

            int mWR = mWR_values[mwr_idx];
            int lineStyle = lineStyles[mwr_idx % lineStyles.size()];

            std::vector<double> x, y;

            for (int mN : mN_values) {
                std::string tb_str;
                if (tb == 0.0 || tb == 1.0)
                    tb_str = Form("%.0f", tb);   // specifically for "0" or "1"
                else
                    tb_str = Form("%.1f", tb);   

                std::ostringstream path;
                path << base
                     << "HNL_decay2_tb"
                     << std::fixed << std::setprecision(1) << tb_str
                     << "_MN" << mN
                     << "_MWR" << mWR
                     << "/Events/run_01/unweighted_events.lhe.gz";

                double sigma = read_xsec_from_lhe_gz(path.str());
                std::cout<< sigma << std::endl;
                if (sigma <= 0) continue;

                x.push_back(mN);
                y.push_back(sigma);

                ymin = std::min(ymin, sigma);
                ymax = std::max(ymax, sigma);
            }

            if (x.empty()) continue;

            TGraph *gr = new TGraph(x.size(), x.data(), y.data());
            gr->SetLineColor(color);
            gr->SetMarkerColor(color);
            gr->SetLineStyle(lineStyle);
            gr->SetLineWidth(3);
            gr->SetMarkerStyle(20);

            if (first) {
                gr->SetTitle(
                    "Imact of Chiral mixing #xi_{W} on HNL Decay;"
                    "m_{N} [GeV];"
                    "#sigma(pp#rightarrowlN) #times #Gamma(N#rightarrowljj) [pb]"
                );
                gr->Draw("APL");
                gr->GetYaxis()->SetRangeUser(1e-6, 1e2);
                // gr->GetYaxis()->SetRangeUser(ymin/2, ymax*2);
                first = false;
            } else {
                gr->Draw("PL SAME");
            }
        }
    }

    // Populate Legends

    for (size_t i = 0; i < tb_values.size(); i++) {
        TGraph *g = new TGraph(1);
        g->SetLineColor(colors[i]);
        g->SetLineWidth(3);
        colorLeg->AddEntry(g, Form("tan#beta= %.1f", tb_values[i]), "l");
    }

    for (size_t i = 0; i < mWR_values.size(); i++) {
        TGraph *g = new TGraph(1);
        g->SetLineColor(kBlack);
        g->SetLineStyle(lineStyles[i]);
        g->SetLineWidth(3);
        styleLeg->AddEntry(g, Form("m_{W_{R}}= %dTeV", mWR_values[i]/1000), "l");
    }
// legS->AddEntry(sampS, Form("%.1f TeV", 0.002*S), "lp");
    colorLeg->Draw();
    styleLeg->Draw();
    infoBox->Draw();

    c->cd();
    gPad->Update();
    TPaveText *pt = (TPaveText*)gPad->GetPrimitive("title");
    if (pt) {
        pt->SetTextSize(0.045);   // try 0.04–0.06
        pt->SetTextFont(32);
    }
    c->SaveAs("HNL_decay2_CS.png");
}

