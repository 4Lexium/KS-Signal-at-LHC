// Run with: root -l HNL_LHE_boost_scan.C

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TMath.h"

//------------------------------------------------------------------------------
// Scan list
//------------------------------------------------------------------------------

static const std::vector<double> TB_LIST  = {0.1, 0.3};
static const std::vector<double> MN_LIST  = {100};
static const std::vector<double> MWR_LIST = {3000, 5000, 10000};

static const std::string BASE_PATH = "/storage/user/WORKBENCH/runs";
static const std::string PROC_NAME = "KS_delphes_decay_PythGEOGUSSE77";

static const int HNL_PDG = 9900012;

//------------------------------------------------------------------------------
// Helper: build work directory from (tb, MN, MWR)
// ex: /storage/.../KS_delphes_decay_PythGEOGUSSE77_tb0.1_MN50_MWR10000/filepath.XY
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
// read total width of HNL (PDG 9900012) from param_card.dat
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
// Declare Error statistic 
//------------------------------------------------------------------------------

struct RunningStats {
    long long n = 0;
    double mean = 0.0;
    double M2   = 0.0;

    void add(double x) {
        ++n;
        double delta  = x - mean;
        mean         += delta / n;
        double delta2 = x - mean;
        M2           += delta * delta2;
    }

    double getMean() const {
        return (n > 0) ? mean : 0.0;
    }

    double getStdDev() const {
        return (n > 1) ? std::sqrt(M2 / (n - 1)) : 0.0;
    }
};

struct PointResult {
    double mean_bg    = 0.0;
    double sigma_bg   = 0.0;
    double mean_L_mm  = 0.0;
    double sigma_L_mm = 0.0;
};

//------------------------------------------------------------------------------
// For one point (tb, MN, MWR): read width from paramCard -> compute ctau, read beta, gama (lab frame) from LHE -> compute L
// beta and gama contain error assosiated with kinematic distribution, ctau are theoretical values as constants.
//------------------------------------------------------------------------------

PointResult analyze_point(double tb, double MN, double MWR, std::vector<double>* L_samples) {
    std::string paramCardPath = build_param_card_path(tb, MN, MWR);
    std::string lheFile       = build_lhe_path(tb, MN, MWR);

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

    // Read lhe.gz via ROOT's gSystem and gunzip
    TString cmd = TString::Format("gunzip -c %s", lheFile.c_str());
    FILE* pipe = gSystem->OpenPipe(cmd, "r");
    if (!pipe) {
        std::cerr << "  ERROR: cannot open LHE file via: " << cmd << "\n";
        return result;
    }

    char line[4096];
    bool inEvent        = false;
    bool firstEventLine = false;
    long long nEvents   = 0;
    long long nHNL      = 0;

    RunningStats stats_bg;
    RunningStats stats_L;

    while (fgets(line, sizeof(line), pipe)) {
        std::string s(line);
        // trim
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

        // LHE particle line:
        // id status mother1 mother2 color1 color2 px py pz E m lifetime spin
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

        if (L_samples) {
            L_samples->push_back(L_mm);
        }
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


void HNL_boost() {
    std::cout << "User Test 007:" << std::endl;

    size_t nTB  = TB_LIST.size();
    size_t nMN  = MN_LIST.size();
    size_t nMWR = MWR_LIST.size();

    // -------------------------------
    // Storage
    // -------------------------------
    std::vector<std::vector<std::vector<PointResult>>> results(
        nTB, std::vector<std::vector<PointResult>>(
                 nMN, std::vector<PointResult>(nMWR)));

    std::vector<std::vector<std::vector<std::vector<double>>>> L_samples(
        nTB, std::vector<std::vector<std::vector<double>>>(
                 nMN, std::vector<std::vector<double>>(nMWR)));

    // -------------------------------
    // Scan
    // -------------------------------
    for (size_t itb = 0; itb < nTB; ++itb) {
        for (size_t imn = 0; imn < nMN; ++imn) {
            for (size_t imwr = 0; imwr < nMWR; ++imwr) {
                results[itb][imn][imwr] =
                    analyze_point(TB_LIST[itb], MN_LIST[imn], MWR_LIST[imwr],
                                  &L_samples[itb][imn][imwr]);
            }
        }
    }

    gStyle->SetOptStat(0);

    // -------------------------------
    // Robust styling
    // -------------------------------

    // color palette (MWR)
    std::vector<int> color_palette = {
        kRed+1, kBlue+1, kGreen+2, kBlack, kViolet+1
    };

    // linestyle palette (tb)
    std::vector<int> linestyle_palette = {
        1, 2, 3, 4  // solid, dashed, dotted, dash-dot
    };

    // -------------------------------
    // Plotting
    // -------------------------------
    for (size_t imn = 0; imn < nMN; ++imn) {

        double MN = MN_LIST[imn];

        // --- Find global max L
        double maxL = 0.0;
        for (size_t itb = 0; itb < nTB; ++itb) {
            for (size_t imwr = 0; imwr < nMWR; ++imwr) {
                for (double L : L_samples[itb][imn][imwr]) {
                    if (L > maxL) maxL = L;
                }
            }
        }
        if (maxL <= 0.0) maxL = 1.0;
         
        // --- Canvas
        TCanvas* c = new TCanvas(
            Form("c_mN_%.0f", MN),
            Form("L distributions mN=%.0f GeV", MN),
            900, 700
        );

        c->SetGrid();
        c->SetLogy();
        c->SetLogx();

        TLegend* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        int nbins = 50;
        double xmin = maxL * 1e-4;   // log-safe
        double xmax = maxL * 1.2;

        bool first = true;

        // --- Loop over TB and MWR
        for (size_t itb = 0; itb < nTB; ++itb) {
            double tb = TB_LIST[itb];

            int linestyle = linestyle_palette[itb % linestyle_palette.size()];

            for (size_t imwr = 0; imwr < nMWR; ++imwr) {
                double MWR = MWR_LIST[imwr];

                int color = color_palette[imwr % color_palette.size()];

                TH1D* h = new TH1D(
                    Form("h_mN%.0f_tb%.2f_MWR%.0f", MN, tb, MWR),
                    "",
                    nbins, xmin, xmax
                );

                // Fill
                for (double L : L_samples[itb][imn][imwr]) {
                    if (L > 0) h->Fill(L);
                }

                // --- Style ---
                h->SetLineColor(color);
                h->SetLineWidth(2);
                h->SetLineStyle(linestyle);

                h->SetFillColorAlpha(color, 0.35); // always filled
                h->SetFillStyle(1001);

                // --- Draw ---
                if (first) {
                    h->SetTitle(Form("L = #beta#gamma c#tau_{N} (m_{N}=%.0f GeV);L [mm];Events", MN));
                    h->SetMaximum(5e3);
                    h->Draw("HIST");
                    first = false;
                } else {
                    h->Draw("HIST SAME");
                }

                // --- Legend ---
                std::ostringstream lab;
                lab << "MWR=" << MWR/1000.0 << " TeV, tan#beta=" << tb;

                leg->AddEntry(h, lab.str().c_str(), "lf");
            }
        }

        leg->Draw();
        c->SaveAs(Form("HNL_LHE_Ldist_mN%.0f.png", MN));
        std::cout << "Saved HNL_LHE_Ldist_mN" << MN << ".png\n";
    }


    // Plot 2 normalize by the integral for each hist

    std::vector<double> TB_SELECT  = {};   // e.g. {0.1, 0.3}
    std::vector<double> MWR_SELECT = {};   // e.g. {3000, 10000}
    std::vector<double> MN_SELECT  = {}; 
    auto is_selected = [](double val, const std::vector<double>& sel) {
        if (sel.empty()) return true;
        for (double x : sel)
            if (std::fabs(x - val) < 1e-6) return true;
        return false;
    };

    for (size_t imn = 0; imn < nMN; ++imn) {
        double MN = MN_LIST[imn];
        if (!is_selected(MN, MN_SELECT)) continue;

        TCanvas* c2 = new TCanvas(
            Form("c_norm_mN%.0f", MN),
            Form("Normalized overlay mN=%.0f", MN),
            900, 700
        );

        c2->SetLogx();
        c2->SetGrid();

        TLegend* leg2 = new TLegend(0.55, 0.55, 0.88, 0.88);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);

        int nbins = 50;
        double xmax = 0.0;

        // --- global xmax across selected data
        for (size_t itb = 0; itb < nTB; ++itb) {
            double tb = TB_LIST[itb];
            if (!is_selected(tb, TB_SELECT)) continue;

            for (size_t imwr = 0; imwr < nMWR; ++imwr) {
                double MWR = MWR_LIST[imwr];
                if (!is_selected(MWR, MWR_SELECT)) continue;

                for (double L : L_samples[itb][imn][imwr]) {
                    if (L > xmax) xmax = L;
                }
            }
        }

        if (xmax <= 0) xmax = 1.0;
        double xmin = xmax * 1e-4;
        xmax *= 1.2;

        bool first = true;
        double global_max = 0.0;

        // --- MAIN LOOP (ALL SELECTED COMBINATIONS)
        for (size_t itb = 0; itb < nTB; ++itb) {

            double tb = TB_LIST[itb];
            if (!is_selected(tb, TB_SELECT)) continue;

            int linestyle = linestyle_palette[itb % linestyle_palette.size()];

            for (size_t imwr = 0; imwr < nMWR; ++imwr) {

                double MWR = MWR_LIST[imwr];
                if (!is_selected(MWR, MWR_SELECT)) continue;

                int color = color_palette[imwr % color_palette.size()];

                TH1D* h = new TH1D(
                    Form("h_norm_mN%.0f_tb%.2f_MWR%.0f", MN, tb, MWR),
                    "", nbins, xmin, xmax
                );

                // fill
                for (double L : L_samples[itb][imn][imwr]) {
                    if (L > 0) h->Fill(L);
                }

                // normalize
                double integral = h->Integral();
                if (integral > 0)
                    h->Scale(1.0 / integral);

                // style
                h->SetLineColor(color);
                h->SetLineWidth(2);
                h->SetLineStyle(linestyle);

                h->SetFillColorAlpha(color, 0.25);
                h->SetFillStyle(1001);

                global_max = std::max(global_max, h->GetMaximum());

                // draw
                if (first) {
                    h->SetTitle(
                        Form("Normalized L distributions (m_{N}=%.0f GeV);L [mm];1/N dN/dL", MN)
                    );
                    h->SetMinimum(0);
                    h->Draw("HIST");
                    first = false;
                } else {
                    h->Draw("HIST SAME");
                }

                // legend
                leg2->AddEntry(
                    h,
                    Form("MWR=%.1f TeV, tb=%.2f", MWR/1000.0, tb),
                    "lf"
                );
            }
        }

        // fix max AFTER loop
        if (!first) {
            gPad->Update();
            auto htmp = (TH1*)gPad->GetPrimitive(
                Form("h_norm_mN%.0f_tb%.2f_MWR%.0f", MN, TB_LIST[0], MWR_LIST[0])
            );
            if (htmp) htmp->SetMaximum(global_max * 1.4);
        }

        leg2->Draw();

        c2->SaveAs(Form("HNL_norm_overlay_mN%.0f.png", MN));

        std::cout << "[Plot2] Saved normalized plot for mN=" << MN << "\n";
    }
    std::cout << "Scan Completed Successfully!" << std::endl;
}

