#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TPaveText.h>

#include <zlib.h>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <map>
#include <algorithm>

struct ScanConfig {
    int S;
    std::vector<int> mN_values;
    std::vector<int> mWR_values;
    int colorOffset;
};

std::vector<ScanConfig> scans = {
    {
        6800,
        {100, 500, 1000, 1500},
        {3000, 3500, 4000, 4500, 5000, 6000, 7000, 7500, 10000},
        0
    },
    {
        7000,
        {100, 500, 1000, 1500},
        {3000, 3500, 4000, 4500, 5000, 6000, 7000, 7500, 10000},
        0
    },
    {
        16500,
        {1000, 2000, 5000, 6000},
        {3000, 3500, 4000, 4500, 5000, 6000, 7000, 7500, 10000, 12500, 15000, 17500, 20000},
        4
    }
};


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
    if (!file) return -1.0;

    char buffer[4096];
    bool in_init = false;
    double sigma_tot = 0.0;

    while (gzgets(file, buffer, sizeof(buffer))) {
        std::string line(buffer);

        if (line.find("<init>") != std::string::npos) {
            in_init = true;
            continue;
        }
        if (!in_init) continue;
        if (line.find("</init>") != std::string::npos ||
            line.find("<generator") != std::string::npos)
            break;

        if (line.find("2212") == 0) continue;

        std::istringstream iss(line);
        double sigma, err, sigmax;
        int id;
        if (iss >> sigma >> err >> sigmax >> id)
            sigma_tot += sigma;
    }

    gzclose(file);
    return sigma_tot;
}

void read_KS_CS_mWRmNS() {
    std::vector<TGraph*> allGraphs;
    std::string base = "/storage/alexanum/WORKBENCH/runs/";
    std::cout<< "User Test 009:" <<std::endl;
    std::vector<int> colors = {
        kRed, kBlue, kGreen+2, kBlack,
        kOrange+1, kCyan+1, kViolet, kTeal+1
    };
    std::vector<int> lineStyles = {1, 2, 7};

    TCanvas *c = new TCanvas("c", "Cross sections", 1350, 1000);   //1200 8000
    c->SetLogy();
    c->SetGrid(1,1);

    TPaveText *infoBox = new TPaveText(0.55, 0.722, 0.9, 0.898, "NDC");
    infoBox->SetBorderSize(0);
    // infoBox->SetFillStyle(0);
    infoBox->SetFillColor(0);
    infoBox->SetTextFont(32);
    infoBox->SetTextSize(0.05);
    infoBox->SetTextAlign(12);  // Left aligned
    infoBox->AddText("pp #rightarrow W^{#pm(*)}_{R}#rightarrow l^{#pm}N");
    infoBox->AddText("l#in (e,#mu,#tau), m_{N_{1}}<< m_{N_{2,3}}");

    TPaveText *warning = new TPaveText(0.31, 0.405, 0.38, 0.4405, "NDC");   //0.31, 0.41, 0.38, 0.45
    warning->SetBorderSize(0);
    warning->SetFillStyle(0);
    // warning->SetFillColor(0);
    warning->SetTextFont(32);
    warning->SetTextSize(0.05);
    warning->AddText("m_{N}<~m_{W_{R}}");

    TLegend *legS = new TLegend(0.55, 0.6, 0.9, 0.72);
    legS->SetHeader("LHC com-energy#sqrt{s}");
    legS->SetBorderSize(0);
    legS->SetFillColor(0);
    legS->SetNColumns(2);
    legS->SetTextFont(32);
    legS->SetTextSize(0.04);

    TLegend *legMN = new TLegend(0.55, 0.445, 0.9, 0.59);
    legMN->SetHeader("HNL mass m_{N}");
    legMN->SetBorderSize(0);
    legMN->SetFillColor(0);
    legMN->SetNColumns(2);
    legMN->SetTextFont(32);
    legMN->SetTextSize(0.04);

    // TLegend *legMNext = new TLegend(0.15, 0.15, 0.45, 0.30);
    // legMNext->SetHeader("Extended m_{N} (16.5 TeV)");
    // legMNext->SetBorderSize(0);
    // legMNext->SetFillColor(0);
    // legMNext->SetTextFont(32);
    // legMNext->SetTextSize(0.035);

    double ymin = 1e30;
    bool frameDrawn = false;

    for (size_t s_idx = 0; s_idx < scans.size(); ++s_idx) {

        const auto& scan = scans[s_idx];
        int S = scan.S;
        int lineStyle = lineStyles[s_idx % lineStyles.size()];

        for (size_t n_idx = 0; n_idx < scan.mN_values.size(); ++n_idx) {

            int mN = scan.mN_values[n_idx];
            int color = colors[(n_idx + scan.colorOffset) % colors.size()];

            std::vector<double> x, y;

            for (int mWR : scan.mWR_values) {
                std::ostringstream path;
                path << base
                     << "KSDY_offWR_S" << S
                     << "_MN" << mN
                     << "_MWR" << mWR
                     << "/Events/run_01/unweighted_events.lhe.gz";

                double sigma = read_xsec_from_lhe_gz(path.str());
                if (sigma > 0) {
                    x.push_back(mWR / 1000.0);
                    y.push_back(sigma);
                    ymin = std::min(ymin, sigma);
                }
            }

            if (x.empty()) continue;

            TGraph *gr = new TGraph(x.size(), x.data(), y.data());
            gr->SetLineColor(color);
            gr->SetMarkerColor(color);
            gr->SetLineStyle(lineStyle);
            gr->SetLineWidth(3);
            gr->SetMarkerStyle(20 + s_idx);

            if (!frameDrawn) {
                gr->SetTitle(Form("Cross Section for Drell-Yan Production of HNL at LHC;m_{W_{R}} [TeV]; #sigma(pp#rightarrowlN) [pb]"));
                gr->Draw("APL");
                gr->GetXaxis()->SetLimits(2.5, 18);
                gr->GetYaxis()->SetRangeUser(0.5*ymin, 1e-1);  //0.5*ymin
                frameDrawn = true;
            } else {
                gr->Draw("PL SAME");
            }
            allGraphs.push_back(gr);


            if (scan.colorOffset == 0 && s_idx == 0) {
                TGraph *samp = new TGraph(1);
                samp->SetLineColor(color);
                samp->SetLineWidth(3);
                legMN->AddEntry(samp, Form("%d GeV", mN), "l");
            }

            if (scan.colorOffset > 0) {
                TGraph *samp = new TGraph(1);
                samp->SetLineColor(color);
                samp->SetLineWidth(3);
                legMN->AddEntry(samp, Form("%d GeV", mN), "l");
            }
        }

        TGraph *sampS = new TGraph(1);
        sampS->SetLineColor(kBlack);
        sampS->SetLineStyle(lineStyle);
        sampS->SetLineWidth(3);
        sampS->SetMarkerStyle(20 + s_idx);
        sampS->SetMarkerColor(kBlack);
        legS->AddEntry(sampS, Form("%.1f TeV", 0.002*S), "lp");
    }

    legS->Draw();
    legMN->Draw();
    infoBox->Draw();
    warning->Draw();
    // legMNext->Draw();
    TPad *inset = new TPad("inset", "Zoom", 0.5, 0.1,0.899, 0.45);
    inset->SetFillColor(0);
    inset->SetBorderMode(0);
    inset->SetLeftMargin(0.09);   // default is ~0.15, smaller = tighter
    inset->SetRightMargin(0.01);  // default ~0.15
    inset->SetTopMargin(0.04);    // default ~0.1
    inset->SetBottomMargin(0.08);
    inset->SetLogy();
    inset->SetGrid(1,0);
    inset->Draw();
    inset->cd();

    bool insetFrame = false;

    for (TGraph *gr : allGraphs) {

        // Clone so axis limits don't affect main plot
        TGraph *g = (TGraph*)gr->Clone();

        if (!insetFrame) {
            g->SetTitle(""); 
            g->Draw("APL");
            g->GetXaxis()->SetTitle("m_{W_{R}} [TeV]");
            g->GetYaxis()->SetTitle("#sigma [pb]");
            g->GetXaxis()->SetTitleOffset(-0.5);
            g->GetYaxis()->SetTitleOffset(-0.5);
            g->GetXaxis()->CenterTitle();
            g->GetYaxis()->CenterTitle();
            g->GetXaxis()->SetLimits(10.2, 20.2);
            g->GetYaxis()->SetRangeUser(1e-7, 1e-4);
            g->GetXaxis()->SetTitleSize(0.10);
            g->GetYaxis()->SetTitleSize(0.10);
            g->GetXaxis()->SetLabelSize(0.08);
            g->GetYaxis()->SetLabelSize(0.08);
            insetFrame = true;
        } else {
            g->Draw("PL SAME");
        }
    }
    c->cd();
    gPad->Update();
    TPaveText *pt = (TPaveText*)gPad->GetPrimitive("title");
    if (pt) {
        pt->SetTextSize(0.045);   // try 0.04–0.06
        pt->SetTextFont(32);
    }

    c->SaveAs("CS_KS_mWR_mN_S_ext.png");
}
