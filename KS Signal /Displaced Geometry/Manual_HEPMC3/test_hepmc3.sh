#!/usr/bin/env bash
set -e

# Paths to installations
PYTHIA_PATH=/storage/user/madgraph/MG5_aMC_v2_9_25/HEPTools/pythia8317
HEPMC3_BASE=/storage/user/madgraph/MG5_aMC_v2_9_25/HEPTools/hepmc3/HepMC3-3.3.1

# MG5 run directory to process:
WORKDIR=/storage/user/WORKBENCH/runs/KS_delphes_decay_PythGEOGUSSE4_MN50_MWR5000

# ----------------------------------------------------------------------
# 1. Compile the converter if the binary does not exist
# ----------------------------------------------------------------------

echo "Compiling lhe2hepmc3..."
g++ -O2 -std=c++11 \
  -I"$PYTHIA_PATH/include" \
  -I"$HEPMC3_BASE/include" \
  -L"$PYTHIA_PATH/lib" \
  -L"$HEPMC3_BASE/build/outputs/lib64" \
  /home/user/WORKSPACE/lhe2hepmc3.cc \
  -o /home/user/WORKSPACE/lhe2hepmc3 \
  -lpythia8 -lHepMC3
echo "Compilation done."

# ----------------------------------------------------------------------
# 2. Environment for runtime linking
# ----------------------------------------------------------------------
export LD_LIBRARY_PATH="$PYTHIA_PATH/lib:$HEPMC3_BASE/build/outputs/lib64:$LD_LIBRARY_PATH"

# ----------------------------------------------------------------------
# 3. Locate input card and LHE, set output name
# ----------------------------------------------------------------------
CARD="$WORKDIR/Cards/pythia8_card.dat"
LHE_GZ="$WORKDIR/Events/run_01/unweighted_events.lhe.gz"
LHE="$WORKDIR/Events/run_01/unweighted_events.lhe"
HEPMC="$WORKDIR/Events/run_01/events.hepmc3"

MN=50
MWR=5000
# Calculate Numerical value for HNL width: Gamma(mN, mWR)
DECAY_WIDTH=$(awk -v MN=$MN -v MWR=$MWR 'BEGIN{printf "%.12e", 1e-5*MN^5/MWR^4}')
DECAY_LT=$(awk -v GAMA=$DECAY_WIDTH 'BEGIN{printf "%.12e", 1.973e-13/GAMA}')
mkdir -p "$WORKDIR/Cards"
sed -i 's/^DECAY 9900012 .*/DECAY 9900012 0.0 # WN1/' $WORKDIR/Cards/param_card.dat
cat << PYTHIA_EOF > $WORKDIR/Cards/pythia8_card.dat
SLHA:useDecayTable = off
9900012:mayDecay = on
9900012:isResonance = false
9900012:doExternalDecay = off
9900012:tau0 = $DECAY_LT
9900012:addChannel = 1 0.125 0 11 -1 2
9900012:addChannel = 1 0.125 0 11 -3 4
9900012:addChannel = 1 0.125 0 13 -1 2
9900012:addChannel = 1 0.125 0 13 -3 4
9900012:addChannel = 1 0.125 0 -11 1 -2
9900012:addChannel = 1 0.125 0 -11 3 -4
9900012:addChannel = 1 0.125 0 -13 1 -2
9900012:addChannel = 1 0.125 0 -13 3 -4
ParticleDecays:limitTau0 = on
ParticleDecays:tau0Max = 1e6
HadronLevel:Hadronize = on
HadronLevel:Decay = on
Next:numberShowEvent = 0
PYTHIA_EOF

if [ -f "$LHE_GZ" ]; then
  echo "Found gzipped LHE: $LHE_GZ"
  gunzip -c "$LHE_GZ" > "$LHE"
elif [ ! -f "$LHE" ]; then
  echo "ERROR: No LHE file found at:"
  echo "  $LHE_GZ"
  echo "  $LHE"
  exit 1
fi

echo "Using card: $CARD"
echo "Using LHE : $LHE"
echo "Output    : $HEPMC"

# ----------------------------------------------------------------------
# 4. Run the converter
# ----------------------------------------------------------------------
/home/user/WORKSPACE/lhe2hepmc3 \
  "$CARD" \
  "$LHE" \
  "$HEPMC"
echo "Done. HepMC3 file at: $HEPMC"

# ----------------------------------------------------------------------
# 5. Analyze HNL displacement with HepMC3 + ROOT
# ----------------------------------------------------------------------
SCAN_BIN=/home/user/WORKSPACE/scan_HNLDV
HISTROOT="$WORKDIR/Events/run_01/HNL_displacement.root"

echo "Compiling scan_HNLDV..."
g++ -O2 -std=c++17 \
  -I"$HEPMC3_BASE/include" \
  -L"$HEPMC3_BASE/build/outputs/lib64" \
  /home/user/WORKSPACE/scan_HNLDV.cc \
  -o "$SCAN_BIN" \
  `root-config --cflags --libs` \
  -lHepMC3
echo "Compilation of scan_HNLDV done."

echo "Running HNL displacement analysis..."
"$SCAN_BIN" "$HEPMC" "$HISTROOT" 9900012
echo "HNL displacement histograms saved in $HISTROOT"

