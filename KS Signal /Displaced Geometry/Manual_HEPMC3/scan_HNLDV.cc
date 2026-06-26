// scan_HNLDV.cc
#include <iostream>
#include <memory>
#include <string>

// HepMC3
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/ReaderAscii.h>

// ROOT
#include "TFile.h"

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " events.hepmc3 output.root PDG\n";
        return 1;
    }

    std::string hepmcFile = argv[1];
    std::string rootFile  = argv[2];
    int targetPDG         = std::stoi(argv[3]);

    HepMC3::ReaderAscii reader(hepmcFile);
    if (reader.failed()) {
        std::cerr << "ERROR: Could not open HepMC3 file " << hepmcFile << "\n";
        return 1;
    }

    TFile *outFile = TFile::Open(rootFile.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "ERROR: Could not create ROOT file " << rootFile << "\n";
        return 1;
    }

    std::shared_ptr<HepMC3::GenEvent> evt;
    int ievt = 0;

    while (!reader.failed()) {
        evt = std::make_shared<HepMC3::GenEvent>();
        reader.read_event(*evt);
        if (reader.failed()) break;
        ++ievt;

        std::cout << "================ Event " << ievt << " ================\n";

        for (auto &p : evt->particles()) {
            if (!p) continue;
            if (p->pid() != targetPDG) continue;

            int barcode = p->id();
            std::cout << "Particle pid=" << p->pid()
                      << " barcode=" << barcode << "\n";

            // Production vertex
            std::shared_ptr<HepMC3::GenVertex> v_prod = p->production_vertex();
            if (v_prod) {
                auto x = v_prod->position();
                std::cout << "  prod vertex id=" << v_prod->id()
                          << " pos=("
                          << x.x() << ", "
                          << x.y() << ", "
                          << x.z() << ", "
                          << x.t() << ") [mm]\n";
            } else {
                std::cout << "  prod vertex: NONE\n";
            }

            // Decay vertex
            std::shared_ptr<HepMC3::GenVertex> v_end = p->end_vertex();
            if (v_end) {
                auto x = v_end->position();
                std::cout << "  decay vertex id=" << v_end->id()
                          << " pos=("
                          << x.x() << ", "
                          << x.y() << ", "
                          << x.z() << ", "
                          << x.t() << ") [mm]\n";
            } else {
                std::cout << "  decay vertex: NONE\n";
            }

            std::cout << "\n";
        }
    }

    reader.close();
    outFile->Write();
    outFile->Close();
    return 0;
}
