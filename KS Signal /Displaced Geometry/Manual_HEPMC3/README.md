Idea: we want to produce hepmc3 output, but MadGraph does a fallback to hepmc2 (if we run Pythia via run cards)

Instead use MadGraph to produce lhe (regular steps) 

Note: generate_KSDV-mg5 has updates N1 > n1 and WR > wr , *important!*

Run Pythia 8.3.17 separately using lhe2hepmc3.cc where Pythia commands can be defined in the macro. 

The produced events.hepmc3 can be directly analysef for HNL-masked vertices using scal_HNLDV.cc which works DIRECTLY on HEPMC3

Conclusion: we achieve non-zero displacement for HNL (see terminal print .png) however it does not get registered for vertex coordinates. 
