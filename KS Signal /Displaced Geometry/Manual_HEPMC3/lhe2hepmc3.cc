#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

#include <iostream>
#include <string>
#include <memory>
#include <iomanip>

using namespace Pythia8;

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " pythia8_card.dat input.lhe output.hepmc3\n";
        return 1;
    }

    std::string cardFile   = argv[1];
    std::string lheFile    = argv[2];
    std::string hepmcFile  = argv[3];

    Pythia pythia;
    pythia.readFile(cardFile);

    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = " + lheFile);

    // pythia.readString("9900012:new = nuN nuNbar 2 0 0 50 5.0e-12 0.0 0.0 $DECAY_LT");
    // pythia.readString("9900012:mayDecay = on");
    // pythia.readString("9900012:isResonance = off");
    // pythia.readString("9900012:onMode = off");
    // pythia.readString("9900012:addChannel = 1 0.125 0 11 -1 2");
    // pythia.readString("9900012:addChannel = 1 0.125 0 11 -3 4");
    // pythia.readString("9900012:addChannel = 1 0.125 0 13 -1 2");
    // pythia.readString("9900012:addChannel = 1 0.125 0 13 -3 4");
    // pythia.readString("9900012:addChannel = 1 0.125 0 -11 1 -2");
    // pythia.readString("9900012:addChannel = 1 0.125 0 -11 3 -4");
    // pythia.readString("9900012:addChannel = 1 0.125 0 -13 1 -2");
    // pythia.readString("9900012:addChannel = 1 0.125 0 -13 3 -4");
    // pythia.readString("Init:showChangedParticleData = on");

    HepMC3::Pythia8ToHepMC3 toHepMC;
    HepMC3::WriterAscii writer(hepmcFile);

    if (!pythia.init()) {
        std::cerr << "Pythia initialization failed.\n";
        return 1;
    }
    std::cout<< "User Test 007:"<<std::endl;
    const int HNL_PDG = 9900012;
    int iEvent = 0;

    std::cout << std::scientific << std::setprecision(6);

    while (pythia.next()) {
        ++iEvent;
        if (iEvent > 10) break;
        std::cout << "================ Event " << iEvent << " ================\n";

        // Find the *last* HNL in the shower chain (largest index with daughters)
        int bestIndex = -1;
        for (int i = 0; i < pythia.event.size(); ++i) {
            auto &p = pythia.event[i];
            if (p.id() != HNL_PDG) continue;
            // if (p.status() < 0) continue;
            Vec4 vProd = p.vProd();
            Vec4 vDec  = p.vDec();
            std::cout << "Evt " << iEvent
                    << " HNL i=" << i
                    << " status=" << p.status()
                    << " vProd=(" << vProd[0] << ", " << vProd[1] << ", "
                                    << vProd[2] << ", " << vProd[3] << ")"
                    << " vDec=("  << vDec[0]  << ", " << vDec[1]  << ", "
                                    << vDec[2]  << ", " << vDec[3]  << ")\n";
    
            if (p.daughter1() <= 0 || p.daughter2() < p.daughter1()) continue;
            bestIndex = i;
        }

        if (bestIndex >= 0) {
            auto &hnl = pythia.event[bestIndex];
            Vec4 vProd = hnl.vProd();
            Vec4 vDec  = hnl.vDec();
            Vec4 dV    = vDec - vProd;

            std::cout << "HNL at index " << bestIndex
                      << " status=" << hnl.status()
                      << " daughters " << hnl.daughter1()
                      << "..." << hnl.daughter2() << "\n";
            std::cout << "  vProd = (" << vProd[0] << ", "
                                      << vProd[1] << ", "
                                      << vProd[2] << ", "
                                      << vProd[3] << ") mm\n";
            std::cout << "  vDec  = (" << vDec[0] << ", "
                                      << vDec[1] << ", "
                                      << vDec[2] << ", "
                                      << vDec[3] << ") mm\n";
            std::cout << "  dV    = (" << dV[0]  << ", "
                                      << dV[1]  << ", "
                                      << dV[2]  << ", "
                                      << dV[3]  << ") mm\n";

            // Daughters
            std::cout << "  Daughters:\n";
            for (int j = hnl.daughter1(); j <= hnl.daughter2(); ++j) {
                auto &c = pythia.event[j];
                Vec4 vChildProd = c.vProd();
                std::cout << "    idx=" << j
                          << " id=" << c.id()
                          << " status=" << c.status()
                          << " vProd=(" << vChildProd[0] << ", "
                                        << vChildProd[1] << ", "
                                        << vChildProd[2] << ", "
                                        << vChildProd[3] << ") mm\n";
            }
        } else {
            std::cout << "  No decaying HNL (no daughters) in this event.\n";
        }

        auto hepmcEvent = std::make_shared<HepMC3::GenEvent>();
        toHepMC.fill_next_event(pythia, hepmcEvent.get());
        writer.write_event(*hepmcEvent);
    }
    pythia.stat();
    writer.close();
    std::cout<< "Sucess !"<<std::endl;
    return 0;
}
