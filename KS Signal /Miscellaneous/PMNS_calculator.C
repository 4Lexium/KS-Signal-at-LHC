#include <iostream>
#include <cmath>

#include "TROOT.h"
#include "TComplex.h"
#include "TMatrixD.h"
#include "TMath.h"

void PMNS_calculator() {  // Added parentheses for function declaration
    // Clear the ROOT canvas if any
    gROOT->Clear();
    std::cout << "User Test: 004" <<std::endl;
    // Define mixing angles (sine values)
    double s12 = 0.554;  //0.554
    double s13 = 0.148;  //0.148, 0.15
    double s23 = 0.756;  //0.756, 0.1
    double delta_cp = 3.437;

    // Calculate angles in radians
    double theta12 = asin(s12);
    double theta13 = asin(s13);
    double theta23 = asin(s23);

    // Calculate cosines
    double c12 = cos(theta12);
    double c13 = cos(theta13);
    double c23 = cos(theta23);

    // CP-violating phase factor
    TComplex phase_factor(0, delta_cp);
    TComplex exp_idelta = TComplex::Exp(phase_factor);
    
    // Compute PMNS matrix elements
    TComplex U_e1, U_e2, U_e3;
    TComplex U_mu1, U_mu2, U_mu3;
    TComplex U_tau1, U_tau2, U_tau3;
    
    // Electron row
    U_e1 = c12 * c13;
    U_e2 = s12 * c13;
    U_e3 = s13 * TComplex::Conjugate(exp_idelta);
    
    // Muon row
    U_mu1 = -s12*c23 - c12*s23*s13*exp_idelta;
    U_mu2 = c12*c23 - s12*s23*s13*exp_idelta;
    U_mu3 = s23*c13;
    
    // Tau row
    U_tau1 = s12*s23 - c12*c23*s13*exp_idelta;
    U_tau2 = -c12*s23 - s12*c23*s13*exp_idelta;
    U_tau3 = c23*c13;
    
    // Print the PMNS matrix
    std::cout << "\n==========================================" << std::endl;
    std::cout << "        PMNS-R MATRIX" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "\nInput parameters:" << std::endl;
    std::cout << "  s12 = " << s12 << std::endl;
    std::cout << "  s23 = " << s23 << std::endl;
    std::cout << "  s13 = " << s13 << std::endl;
    std::cout << "  δCP = " << delta_cp << std::endl;
    
    std::cout << "\nPMNS-R Matrix Elements:" << std::endl;
    std::cout << "U_e1  = " << Form("%8.6f", U_e1.Re()) << std::endl;
    std::cout << "U_e2  = " << Form("%8.6f", U_e2.Re()) << std::endl;
    std::cout << "U_e3  = " << Form("%8.6f + i%8.6f", U_e3.Re(), U_e3.Im()) << std::endl;
    std::cout << "U_μ1  = " << Form("%8.6f + i%8.6f", U_mu1.Re(), U_mu1.Im()) << std::endl;
    std::cout << "U_μ2  = " << Form("%8.6f + i%8.6f", U_mu2.Re(), U_mu2.Im()) << std::endl;
    std::cout << "U_μ3  = " << Form("%8.6f", U_mu3.Re()) << std::endl;
    std::cout << "U_τ1  = " << Form("%8.6f + i%8.6f", U_tau1.Re(), U_tau1.Im()) << std::endl;
    std::cout << "U_τ2  = " << Form("%8.6f + i%8.6f", U_tau2.Re(), U_tau2.Im()) << std::endl;
    std::cout << "U_τ3  = " << Form("%8.6f", U_tau3.Re()) << std::endl;

    std::cout << "Unitarity Condition: N1 Column = ";
    std::cout << (TComplex::Abs(U_e1)*TComplex::Abs(U_e1) + TComplex::Abs(U_mu1)*TComplex::Abs(U_mu1) + TComplex::Abs(U_tau1)*TComplex::Abs(U_tau1)) << std::endl;
    
    // Print magnitudes
    std::cout << "\nPMNS-R Matrix Elements:" << std::endl;
    std::cout << "        N₁          N₂          N₃" << std::endl;
    std::cout << "e " << Form("%10.6f  %10.6f  %10.6f", 
                          TComplex::Abs(U_e1), TComplex::Abs(U_e2), TComplex::Abs(U_e3)) << std::endl;
    std::cout << "μ " << Form("%10.6f  %10.6f  %10.6f", 
                          TComplex::Abs(U_mu1), TComplex::Abs(U_mu2), TComplex::Abs(U_mu3)) << std::endl;
    std::cout << "τ " << Form("%10.6f  %10.6f  %10.6f", 
                          TComplex::Abs(U_tau1), TComplex::Abs(U_tau2), TComplex::Abs(U_tau3)) << std::endl;
}