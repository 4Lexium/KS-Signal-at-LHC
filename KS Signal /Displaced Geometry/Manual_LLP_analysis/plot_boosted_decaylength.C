// Run in terminal with: root -l HNL_LHE_boost2.C

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>

// ROOT
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TSystem.h"

//------------------------------------------------------------------------------
// Declare File lists to scan
//------------------------------------------------------------------------------

static const std::vector<double> TB_LIST  = {0.9};   // gammabeta ~ independent of tb   0.1, 0.3, 0.5, 0.9
static const std::vector<double> MN_LIST  = {50, 70, 80, 100, 120, 150, 250};  // 50, 100, 250, 375, 500, 750, 1000
static const std::vector<double> MWR_LIST = {3000, 4000, 5000, 6000, 10000, 15000, 25000, 50000}; //3000, 4000, 5000, 6000, 10000, 15000, 25000, 50000

static const std::string BASE_PATH = "/storage/user/WORKBENCH/runs";
static const std::string PROC_NAME = "KS_delphes_decay_PythGEOGUSSE88";

static const int HNL_PDG = 9900012;

//------------------------------------------------------------------------------
// Function: build work directory for (tb, MN, MWR) combo
// ex: /storage/path.../KS_delphes_decay_PythGEOGUSSE77_tb0.1_MN50_MWR10000
// then append /Cards fpr param_card or /Events for lhe.gz
//------------------------------------------------------------------------------

std::string build_workdir(double tb, double MN, double MWR) {
    std::ostringstream dir;
    dir << BASE_PATH << "/"
        << PROC_NAME
        << "_tb"  << tb
        << "_MN"  << MN
        << "_MWR" << MWR;
    return dir.str();
}

std::string build_param_card_path(double tb, double MN, double MWR) {
    return build_workdir(tb, MN, MWR) + "/Cards/param_card.dat";
}

std::string build_lhe_path(double tb, double MN, double MWR) {
    return build_workdir(tb, MN, MWR) + "/Events/run_01/unweighted_events.lhe.gz";
}

//------------------------------------------------------------------------------
// read total width for HNL (PDG: 9900012) from param_card.dat
//------------------------------------------------------------------------------

double read_HNL_width_9900012(const std::string& paramCardPath) {
    std::ifstream fin(paramCardPath);
    if (!fin) {
        throw std::runtime_error("Cannot open param_card.dat at: " + paramCardPath);
    }

    std::string line;
    const std::string target = "DECAY  9900012";

    while (std::getline(fin, line)) {
        std::string trimmed = line;
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        if (trimmed.empty()) continue;

        if (trimmed.rfind(target, 0) == 0) {
            std::istringstream iss(trimmed);
            std::string decayWord;
            int pdg;
            double width;
            iss >> decayWord >> pdg >> width;
            if (!iss.fail() && pdg == 9900012) {
                return width; // GeV
            }
        }
    }

    throw std::runtime_error("DECAY  9900012 line not found in param_card: " + paramCardPath);
}

//------------------------------------------------------------------------------
// Error statistic for samples
//------------------------------------------------------------------------------

struct RunningStats {
    long long n = 0;
    double mean = 0.0;
    double M2 = 0.0;

    void add(double x) {
        ++n;
        double delta = x - mean;
        mean += delta / n;
        double delta2 = x - mean;
        M2 += delta * delta2;
    }

    double getMean() const {
        return (n > 0) ? mean : 0.0;
    }

    double getStdDev() const {
        return (n > 1) ? std::sqrt(M2 / (n - 1)) : 0.0;
    }
};

struct PointResult {
    double mean_bg = 0.0;
    double sigma_bg = 0.0;
    double mean_L_mm = 0.0;
    double sigma_L_mm = 0.0;
};

//------------------------------------------------------------------------------
// For one point (tb, MN, MWR): 
//      read param card: read width, compute ctau
//      read lhe: compute mean and sigma of beta*gamma
//      combine to: L = beta*gamma*c*tau (with err)
//------------------------------------------------------------------------------

PointResult analyze_point(double tb, double MN, double MWR) {
    std::string paramCardPath = build_param_card_path(tb, MN, MWR);
    std::string lheFile = build_lhe_path(tb, MN, MWR);

    std::cout << "\n=== Analyzing tb=" << tb
              << ", mN=" << MN
              << ", MWR=" << MWR << " ===\n";
    std::cout << "param_card: " << paramCardPath << "\n";
    std::cout << "LHE file:   " << lheFile       << "\n";

    PointResult result;

    // --- width and c*tau ---
    double gammaN = 0.0;
    try {
        gammaN = read_HNL_width_9900012(paramCardPath);
    } catch (const std::exception& ex) {
        std::cerr << "  ERROR reading width: " << ex.what() << "\n";
        return result;
    }

    if (gammaN <= 0.0) {
        std::cerr << "  WARNING: non-positive width " << gammaN << " at this point.\n";
        return result;
    }

    const double hbarc_mm = 1.97327e-13; // GeV * mm
    double ctau_mm = hbarc_mm / gammaN;

    std::cout << "  Gamma_N = " << gammaN << " GeV, c*tau_N = " << ctau_mm << " mm\n";

    // --- open LHE.gz via ROOT's gSystem and gunzip ---
    TString cmd = TString::Format("gunzip -c %s", lheFile.c_str());
    FILE* pipe = gSystem->OpenPipe(cmd, "r");
    if (!pipe) {
        std::cerr << "  ERROR: cannot open LHE file via: " << cmd << "\n";
        return result;
    }

    char line[4096];
    bool inEvent = false;
    bool firstEventLine = false;
    long long nEvents = 0;
    long long nHNL = 0;

    RunningStats stats_bg;
    RunningStats stats_L;

    while (fgets(line, sizeof(line), pipe)) {
        std::string s(line);
        s.erase(0, s.find_first_not_of(" \t\r\n"));
        if (!s.empty()) {
            size_t pos = s.find_last_not_of(" \t\r\n");
            if (pos != std::string::npos) s.erase(pos + 1);
        }
        if (s.empty() || s[0] == '#') continue;

        if (s.find("<event>") != std::string::npos) {
            inEvent = true;
            firstEventLine = true;
            ++nEvents;
            continue;
        }
        if (s.find("</event>") != std::string::npos) {
            inEvent = false;
            continue;
        }
        if (!inEvent) continue;

        // First line inside <event> is header
        if (firstEventLine) {
            firstEventLine = false;
            continue;
        }
        if (s[0] == '<') continue;

        // LHE particle line: id status mother1 mother2 color1 color2 px py pz E m lifetime spin
        int pdgId, status;
        int m1, m2, c1, c2;
        double px, py, pz, E, m, l, sp;

        std::istringstream iss(s);
        if (!(iss >> pdgId >> status >> m1 >> m2 >> c1 >> c2
                  >> px >> py >> pz >> E >> m >> l >> sp)) {
            continue;
        }

        if (pdgId != HNL_PDG) continue;

        double p3 = std::sqrt(px*px + py*py + pz*pz);
        double beta_gamma = p3 / MN;   // use generator mass
        double L_mm       = beta_gamma * ctau_mm;

        stats_bg.add(beta_gamma);
        stats_L.add(L_mm);
        ++nHNL;
    }

    gSystem->ClosePipe(pipe);

    if (nHNL == 0) {
        std::cerr << "  WARNING: no HNL found in LHE events at this point.\n";
        return result;
    }

    result.mean_bg    = stats_bg.getMean();
    result.sigma_bg   = stats_bg.getStdDev();
    result.mean_L_mm  = stats_L.getMean();
    result.sigma_L_mm = stats_L.getStdDev();

    std::cout << "  nEvents = " << nEvents
              << ", nHNL = " << nHNL
              << ", <beta*gamma> = " << result.mean_bg
              << " +/- " << result.sigma_bg
              << ", <L_mm> = " << result.mean_L_mm
              << " +/- " << result.sigma_L_mm << " mm\n";

    return result;
}

//------------------------------------------------------------------------------
// ROOT main function
//------------------------------------------------------------------------------

void plot_boosted_decaylength() {
    std::cout<< "User Test 008:"<<std::endl;
    size_t nTB  = TB_LIST.size();
    size_t nMN  = MN_LIST.size();
    size_t nMWR = MWR_LIST.size();

    std::vector<std::vector<std::vector<PointResult>>> results(
        nTB,
        std::vector<std::vector<PointResult>>(
            nMN,
            std::vector<PointResult>(nMWR)
        )
    );

    // Scan all points
    for (size_t itb = 0; itb < nTB; ++itb) {
        for (size_t imn = 0; imn < nMN; ++imn) {
            for (size_t imwr = 0; imwr < nMWR; ++imwr) {
                double tb  = TB_LIST[itb];
                double MN  = MN_LIST[imn];
                double MWR = MWR_LIST[imwr];
                results[itb][imn][imwr] = analyze_point(tb, MN, MWR);
            }
        }
    }

    gStyle->SetOptStat(0);
    // Style declaration: colours for WR, linestyles for tb
    std::vector<int> colors = {kRed, kOrange+7, kOrange, kOrange+3, kMagenta, kBlue, kCyan, kGreen+2};
    auto color_for_mwr = [&](size_t imwr) -> int {
        if (colors.empty()) return kBlack;
        return colors[imwr % colors.size()];
    };
    std::vector<int> styles = {1, 2}; 
    auto style_for_tb = [&](size_t itb) -> int {
        if (styles.empty()) return 1;
        return styles[itb % styles.size()];
    };

    // ------------------------------------------------------------
    // Plot 1: <beta*gamma> vs mN with error bars
    // ------------------------------------------------------------
    TCanvas* c_bg = new TCanvas("c_bg", "Boost-factor for HNL", 1000, 800);
    // c_bg->SetGrid(1);
    c_bg->SetLogx();
    TLegend* leg_bg = new TLegend(0.5, 0.50, 0.9, 0.88);
    leg_bg->SetBorderSize(0);
    leg_bg->SetTextFont(32);
    // leg_bg->SetNColumns(2);
    leg_bg->SetFillStyle(0);
    leg_bg->SetTextSize(0.035);

    // Define offset pattern for different MWR
    std::vector<double> offsets = {-0.08, -0.06, -0.04, -0.02, 0.02, 0.04, 0.06, 0.08};

    double max_bg = 0.0;
    double min_bg = 1e9;
    for (size_t itb = 0; itb < nTB; ++itb)
        for (size_t imn = 0; imn < nMN; ++imn)
            for (size_t imwr = 0; imwr < nMWR; ++imwr) {
                double val = results[itb][imn][imwr].mean_bg;
                if (val > 0.0) {
                    if (val > max_bg) max_bg = val;
                    if (val < min_bg) min_bg = val;
                }
            }
    if (max_bg <= 0.0) max_bg = 1.0;
    max_bg = max_bg +50;
    if (min_bg <= 0.0 || min_bg >= max_bg) min_bg = max_bg * 1e-2;

    bool first_bg = true;

    // Store graphs and their labels
    std::vector<TGraphErrors*> graphs;
    std::vector<std::string> labels;
    std::vector<int> indices; // Store original indices for ordering

    // First pass: create and draw graphs
    for (size_t itb = 0; itb < nTB; ++itb) {
        for (size_t imwr = 0; imwr < nMWR; ++imwr) {
            
            std::vector<double> x(nMN), y(nMN), ex(nMN), ey(nMN);
            for (size_t imn = 0; imn < nMN; ++imn) {
                double offset_factor = offsets[imwr % offsets.size()];
                double x_offset = MN_LIST[imn] * offset_factor;
                x[imn]  = MN_LIST[imn] + x_offset;
                ex[imn] = 0.0;
                y[imn]  = results[itb][imn][imwr].mean_bg;
                ey[imn] = results[itb][imn][imwr].sigma_bg;
            }
            
            TGraphErrors* gr = new TGraphErrors((int)nMN, x.data(), y.data(), ex.data(), ey.data());
            
            int col = color_for_mwr(imwr);
            int sty = style_for_tb(itb);
            
            gr->SetLineColor(col);
            gr->SetLineWidth(2);
            gr->SetLineStyle(sty);
            gr->SetMarkerColor(col);
            gr->SetMarkerStyle(20 + imwr % 5);
            gr->SetMarkerSize(0.9);
            
            if (first_bg) {
                gr->SetTitle("Sample Spread for HNL Boost-factor #beta_{N}#gamma_{N}");
                gr->GetXaxis()->SetTitle("m_{N} / GeV");
                gr->GetXaxis()->SetTitleOffset(1.2);
                gr->GetXaxis()->SetTitleFont(32);
                gr->GetYaxis()->SetTitle("#beta_{N}#gamma_{N} sample mean #pm stdev.");
                gr->GetYaxis()->SetTitleOffset(1.3);
                gr->GetYaxis()->SetTitleFont(32);
                gr->GetYaxis()->SetRangeUser(min_bg * 0.8, 80);
                gr->GetYaxis()->SetLimits(0, 100);
                gr->GetXaxis()->SetLimits(0, 1040);
                gr->Draw("ALP");
                first_bg = false;
            } else {
                gr->Draw("LP SAME");
            }
            
            graphs.push_back(gr);
            std::ostringstream lab;
            lab.setf(std::ios::fixed);
            lab.precision(0);
            lab << "m_{W_{R}} = " << MWR_LIST[imwr] / 1000.0 << " TeV";
            labels.push_back(lab.str());
            indices.push_back(imwr);
            leg_bg->AddEntry(gr, lab.str().c_str(), "lp");
        }
    }
   
    leg_bg->AddEntry((TObject*)nullptr, "mixing: tb=0.1", "");
    leg_bg->Draw();
    gPad->Update();

    TPaveText *pt = (TPaveText*)gPad->GetPrimitive("title");
    if (pt) {
        pt->SetTextSize(0.045);   
        pt->SetTextFont(32);
    }
    leg_bg->Draw();
    c_bg->SaveAs("/home/user/WORKSPACE/Manual_LLP_analysis/Figures/2HNL_LHE_mean_beta_gamma_vs_mN.png");
    std::cout << "Saved HNL_LHE_mean_beta_gamma_vs_mN.png\n";

    // ------------------------------------------------------------
    // Plot 2: <L> vs mN with error bars (log y)
    // ------------------------------------------------------------
    TCanvas* c_L = new TCanvas("c_L", "<L> vs m_{N}", 900, 700);
    c_L->SetGrid();
    c_L->SetLogy();    

    TLegend* leg_L = new TLegend(0.55, 0.60, 0.88, 0.88);
    leg_L->SetBorderSize(0);
    leg_L->SetFillStyle(0);
    leg_L->SetTextSize(0.035);

    double max_L = 0.0;
    double min_L = 1e30;
    for (size_t itb = 0; itb < nTB; ++itb)
        for (size_t imn = 0; imn < nMN; ++imn)
            for (size_t imwr = 0; imwr < nMWR; ++imwr) {
                double val = results[itb][imn][imwr].mean_L_mm;
                if (val > 0.0) {
                    if (val > max_L) max_L = val;
                    if (val < min_L) min_L = val;
                }
            }
    if (max_L <= 0.0) max_L = 1.0;
    if (min_L <= 0.0 || min_L >= max_L) min_L = max_L*1e-2;

    bool first_L = true;

    for (size_t itb = 0; itb < nTB; ++itb) {
        for (size_t imwr = 0; imwr < nMWR; ++imwr) {

            std::vector<double> x(nMN), y(nMN), ex(nMN), ey(nMN);
            for (size_t imn = 0; imn < nMN; ++imn) {
                x[imn]  = MN_LIST[imn];
                ex[imn] = 0.0;
                y[imn]  = results[itb][imn][imwr].mean_L_mm;
                ey[imn] = results[itb][imn][imwr].sigma_L_mm;
            }

            TGraphErrors* grL = new TGraphErrors((int)nMN, x.data(), y.data(), ex.data(), ey.data());

            int col = color_for_mwr(imwr);
            int sty = style_for_tb(itb);

            grL->SetLineColor(col);
            grL->SetLineWidth(2);
            grL->SetLineStyle(sty);
            grL->SetMarkerColor(col);
            grL->SetMarkerStyle(20);
            grL->SetMarkerSize(0.9);

            if (first_L) {
                grL->SetTitle("<L> = <#beta#gamma> c#tau_{N} vs m_{N};m_{N} [GeV];<L> [mm]");
                grL->GetYaxis()->SetRangeUser(min_L*0.8, max_L*2);

                grL->GetXaxis()->SetTitleSize(0.05);
                grL->GetYaxis()->SetTitleSize(0.05);
                grL->GetXaxis()->SetLabelSize(0.04);
                grL->GetYaxis()->SetLabelSize(0.04);

                grL->Draw("ALP");
                first_L = false;
            } else {
                grL->Draw("LP SAME");
            }

            std::ostringstream lab;
            lab.setf(std::ios::fixed);
            lab.precision(0);
            lab << "M_{WR}=" << MWR_LIST[imwr]/1000.0 << " TeV";
            leg_L->AddEntry(grL, lab.str().c_str(), "lp");
        }
    }
    leg_L->Draw();
    c_L->SaveAs("/home/user/WORKSPACE/Manual_LLP_analysis/Figures/2HNL_LHE_mean_L_vs_mN.png");
    std::cout << "Saved HNL_LHE_mean_L_vs_mN.png\n";

    // ------------------------------------------------------------
    // Plot 3: 2D map of <betagama> as function of mN and MWR
    // ------------------------------------------------------------
    {
        const int nPoints = nMN * nMWR;
        TGraph2D* g2_L = new TGraph2D(nPoints);
        g2_L->SetName("g2_bg2D");
        int ip = 0;
        for (size_t imn = 0; imn < nMN; ++imn) {
            for (size_t imwr = 0; imwr < nMWR; ++imwr) {
                double x = MN_LIST[imn];                  // mN
                double y = MWR_LIST[imwr] / 1000.0;       // MWR in TeV 
                double z = results[0][imn][imwr].mean_bg; // <betagama>
                g2_L->SetPoint(ip, x, y, z);
                ++ip;
            }
        }

        gStyle->SetPalette(kGreyScale);        
        gStyle->SetNumberContours(50);        // more contours for smoother gradient
        TCanvas* c_L2D = new TCanvas("c_L2D","<L> map in (m_{N}, M_{WR})",1000, 800);
        c_L2D->SetRightMargin(0.15);
        c_L2D->SetLeftMargin(0.12);
        c_L2D->SetBottomMargin(0.12);
        // c_L2D->SetLogz();
        g2_L->SetTitle("Boost-factor for HNL: #beta_{N}#gamma_{N}");  //#left(m_{N}, m_{W_{R}}#right)
        g2_L->SetNpx(100);  
        g2_L->SetNpy(100);
        g2_L->Draw("COLZ");
        c_L2D->Update(); 

        TH2D* h2 = g2_L->GetHistogram();
        if (h2) {
            h2->GetXaxis()->SetTitle("m_{N} / GeV");
            h2->GetXaxis()->SetTitleOffset(1.2); 
            h2->GetXaxis()->SetTitleFont(32);
            h2->GetYaxis()->SetTitle("m_{W_{R}} / TeV");
            h2->GetYaxis()->SetTitleOffset(1.3); 
            h2->GetYaxis()->SetTitleFont(32);
            h2->GetXaxis()->SetNdivisions(505);    
            h2->GetYaxis()->SetNdivisions(505);
        }

        // Clone the histogram for contour lines "layer"
        if (h2) {
            TH2D* h2_cont = (TH2D*)h2->Clone("h2_cont");
            const int nLines = 7;
            double levels[nLines] = {2.1, 2.7, 3.4, 4.05, 5.0, 10.0, 20.0}; 
            h2_cont->SetContour(nLines, levels);
            h2_cont->SetLineColor(kRed);
            h2_cont->SetLineWidth(2);
            h2_cont->SetLineStyle(4);
            h2_cont->Draw("CONT3 SAME");
        }

        gPad->Update();
        TPaveText *pt = (TPaveText*)gPad->GetPrimitive("title");
        if (pt) {
            pt->SetTextSize(0.045);   
            pt->SetTextFont(32);
        }
        c_L2D->SaveAs("/home/user/WORKSPACE/Manual_LLP_analysis/Figures/2HNL_BG2D_MN_MWR.png");
        std::cout << "Saved HNL_LHE_mean_L_vs_mN_vs_MWR_2D_TGraph2D.png\n";
    }

    // ------------------------------------------------------------
    // Plot 4:  2D map of <L> as function of mN and MWR 
    // ------------------------------------------------------------
    {
        const int nPoints = nMN * nMWR;
        TGraph2D* g3_L = new TGraph2D(nPoints);
        g3_L->SetName("g3_L2D");
        int ip = 0;
        for (size_t imn = 0; imn < nMN; ++imn) {
            for (size_t imwr = 0; imwr < nMWR; ++imwr) {
                double x = MN_LIST[imn];                    
                double y = MWR_LIST[imwr] / 1000.0;       
                double z = results[0][imn][imwr].mean_L_mm;  
                g3_L->SetPoint(ip, x, y, z);
                ++ip;
            }
        }

        gStyle->SetPalette(kGreyScale);         
        gStyle->SetNumberContours(50);       
        TCanvas* c2_L2D = new TCanvas("c2_L2D","<L> map in (m_{N}, M_{WR})",1000, 800);
        c2_L2D->SetRightMargin(0.15);
        c2_L2D->SetLeftMargin(0.12);
        c2_L2D->SetBottomMargin(0.12);
        c2_L2D->SetLogz();
        g3_L->SetTitle("HNL Decay Lenght: L_{N} / mm = #beta_{N}#gamma_{N}#times c#tau_{N}#; tb = 0.9");  //#left(m_{N}, m_{W_{R}}#right)
        g3_L->SetNpx(100);   
        g3_L->SetNpy(100);
        g3_L->Draw("COLZ");
        c2_L2D->Update(); 

        TH2D* h3 = g3_L->GetHistogram();
        if (h3) {
            h3->GetXaxis()->SetTitle("m_{N} / GeV");
            h3->GetXaxis()->SetTitleOffset(1.2); 
            h3->GetXaxis()->SetTitleFont(32);
            h3->GetYaxis()->SetTitle("m_{W_{R}} / TeV");
            // h3->GetZaxis()->SetTitle("L_{N} / mm");
            // h3->GetZaxis()->SetTitleFont(32);
            h3->GetYaxis()->SetTitleOffset(1.3); 
            h3->GetYaxis()->SetTitleFont(32);
            h3->GetXaxis()->SetNdivisions(505);    
            h3->GetYaxis()->SetNdivisions(505);
        }

        if (h3) {
            TH2D* h3_cont = (TH2D*)h3->Clone("h3_cont");
            const int nLines = 5;
            // double levels[nLines] = {1e-6, 1e-5, 1e-4, 1e-3, 1, 1e3};  
            double levels[nLines] = {1e-5, 1e-3, 1e-2, 1, 1e3};  
            h3_cont->SetContour(nLines, levels);
            h3_cont->SetLineColor(kRed);
            h3_cont->SetLineWidth(2);
            h3_cont->SetLineStyle(4);
            h3_cont->Draw("CONT3 SAME");
        }

        if (h3) {
            TH2D* h3_bandA = (TH2D*)h3->Clone("h3_bandA");
            double zmin = 1e-9;  
            double zmax = 1e-2;  
            // Zero out everything outside the band
            for (int ix = 1; ix <= h3_bandA->GetNbinsX(); ++ix) {
                for (int iy = 1; iy <= h3_bandA->GetNbinsY(); ++iy) {
                    double val = h3_bandA->GetBinContent(ix, iy);
                    if (val < zmin || val > zmax) {
                        h3_bandA->SetBinContent(ix, iy, 0.0);
                    } else {
                        h3_bandA->SetBinContent(ix, iy, 1.0); // constant fill
                    }
                }
            }
            // Style: hatched transparent band
            h3_bandA->SetFillStyle(3354);   // //// pattern
            h3_bandA->SetFillColor(kAzure); // hatch color
            h3_bandA->SetLineColor(kAzure);
            h3_bandA->SetLineWidth(0);
            h3_bandA->Draw("BOX SAME");
        }
        if (h3) {
            TH2D* h3_bandB = (TH2D*)h3->Clone("h3_bandB");
            double zmin = 1;   
            double zmax = 1e3;  
            for (int ix = 1; ix <= h3_bandB->GetNbinsX(); ++ix) {
                for (int iy = 1; iy <= h3_bandB->GetNbinsY(); ++iy) {
                    double val = h3_bandB->GetBinContent(ix, iy);
                    if (val < zmin || val > zmax) {
                        h3_bandB->SetBinContent(ix, iy, 0.0);
                    } else {
                        h3_bandB->SetBinContent(ix, iy, 1.0); 
                    }
                }
            }
            h3_bandB->SetFillStyle(3354);   
            h3_bandB->SetFillColor(kGreen); 
            h3_bandB->SetLineColor(kGreen);
            h3_bandB->SetLineWidth(0);
            h3_bandB->Draw("BOX SAME");
        }
        if (h3) {
            TH2D* h3_bandC = (TH2D*)h3->Clone("h3_bandC");
            double zmin = 1e3;   
            double zmax = 1e4;   
            for (int ix = 1; ix <= h3_bandC->GetNbinsX(); ++ix) {
                for (int iy = 1; iy <= h3_bandC->GetNbinsY(); ++iy) {
                    double val = h3_bandC->GetBinContent(ix, iy);
                    if (val < zmin || val > zmax) {
                        h3_bandC->SetBinContent(ix, iy, 0.0);
                    } else {
                        h3_bandC->SetBinContent(ix, iy, 1.0); 
                    }
                }
            }
            h3_bandC->SetFillStyle(3354);   
            h3_bandC->SetFillColor(kOrange+7); 
            h3_bandC->SetLineColor(kOrange+7);
            h3_bandC->SetLineWidth(0);
            h3_bandC->Draw("BOX SAME");
        }
        if (h3) {
            TH2D* h3_bandD1 = (TH2D*)h3->Clone("h3_bandD1");
            TH2D* h3_bandD2 = (TH2D*)h3->Clone("h3_bandD2");
            double zmin = 1e-2;   
            double zmax = 1;   
            for (int ix = 1; ix <= h3->GetNbinsX(); ++ix) {
                for (int iy = 1; iy <= h3->GetNbinsY(); ++iy) {
                    double val = h3->GetBinContent(ix, iy);
                    if (val < zmin || val > zmax) {
                        h3_bandD1->SetBinContent(ix, iy, 0.0);
                        h3_bandD2->SetBinContent(ix, iy, 0.0);
                    } else {
                        h3_bandD1->SetBinContent(ix, iy, 1.0); 
                        h3_bandD2->SetBinContent(ix, iy, 1.0); 
                    }
                }
            }
            h3_bandD1->SetFillStyle(3354);   
            h3_bandD1->SetFillColor(kAzure); 
            h3_bandD1->SetLineColor(kAzure);
            h3_bandD1->SetLineWidth(0);
            h3_bandD1->Draw("BOX SAME");
            h3_bandD2->SetFillStyle(3345);   
            h3_bandD2->SetFillColor(kGreen); 
            h3_bandD2->SetLineColor(kGreen);
            h3_bandD2->SetLineWidth(0);
            h3_bandD2->Draw("BOX SAME");
        }

        // TLatex latexInfo;
        // latexInfo.SetTextSize(0.035);
        // latexInfo.SetTextColor(kRed);
        // latexInfo.SetTextFont(42);
        // latexInfo.DrawLatex(350, 7, "nm");
        // latexInfo.DrawLatex(233, 8, "1 #mum");
        // latexInfo.DrawLatex(160, 13, "10#mum");
        // latexInfo.DrawLatex(108, 18, "1mm");
        // latexInfo.DrawLatex(70, 30, "1m");
        TLatex latexInfo_mm;
        latexInfo_mm.SetTextSize(0.033);
        latexInfo_mm.SetTextFont(42);
        latexInfo_mm.SetTextColor(kBlack);
        latexInfo_mm.DrawLatex(270, 36.5, "mm");  //270
        TLegend* leg_band = new TLegend(0.45, 0.7, 0.85, 0.9);
        leg_band->SetBorderSize(1);
        leg_band->SetLineColor(kBlack); 
        leg_band->SetLineWidth(1); 
        leg_band->SetFillStyle(1001);
        leg_band->SetFillColor(kWhite);
        leg_band->SetTextSize(0.03);
        leg_band->SetTextFont(32);
        TH1F* h_dummyA = new TH1F("h_dummyA","",1,0,1);
        h_dummyA->SetFillStyle(3354);     
        h_dummyA->SetFillColor(kAzure);   
        h_dummyA->SetLineColor(kBlack);
        TH1F* h_dummyB = new TH1F("h_dummyB","",1,0,1);
        h_dummyB->SetFillStyle(3354);    
        h_dummyB->SetFillColor(kGreen);  
        h_dummyB->SetLineColor(kBlack);
        TH1F* h_dummyC = new TH1F("h_dummyC","",1,0,1);
        h_dummyC->SetFillStyle(3354);     
        h_dummyC->SetFillColor(kOrange+7);   
        h_dummyC->SetLineColor(kBlack);
        // TH1F* h_dummyD = new TH1F("h_dummyC","",1,0,1);
        // h_dummyC->SetFillStyle(3354);     
        // h_dummyC->SetFillColor(kAzure);
        // h_dummyC->SetFillStyle(3345);     
        // h_dummyC->SetFillColor(kGreen);     
        // h_dummyC->SetLineColor(kBlack);
    
        leg_band->AddEntry(h_dummyA, "Prompt region#left[ < 10#mum#right]", "f");
        leg_band->AddEntry(h_dummyB, "Displaced region#left[ 0.1mm - 1m#right]", "f");
        leg_band->AddEntry(h_dummyC, "Ghost region#left[ >1m#right]", "f");
        leg_band->Draw();

        gPad->Update();
        TPaveText *pt = (TPaveText*)gPad->GetPrimitive("title");
        if (pt) {
            pt->SetTextSize(0.045);   
            pt->SetTextFont(32);
        }
        c2_L2D->SaveAs("/home/user/WORKSPACE/Manual_LLP_analysis/Figures/2HNL_L2D_MN_MWR_tb0.9.png");
        std::cout << "HNL_L2D_MN_MWR_tb0.9.png\n";
    }
    std::cout<<"Plots Generated Successfully !" <<std::endl;
}
