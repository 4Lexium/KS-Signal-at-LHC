#!/usr/bin/env bash
set -e

#=============================
# 0. Paths and base settings
#=============================
MG5=/storage/user/madgraph/MG5_aMC_v2_9_25/bin/mg5_aMC
PYTHIA_PATH=/storage/user/madgraph/MG5_aMC_v2_9_25/HEPTools/pythia8317
HEPMC3_BASE=/storage/user/madgraph/MG5_aMC_v2_9_25/HEPTools/hepmc3/HepMC3-3.3.1

BASE=/storage/user/WORKBENCH/runs
PROC=KS_delphes_decay_PythGEOGUSSE88  

# Mass scan 
# MN_LIST="50 70 100 150 250 375 500"       #  100 500 1000 , 250 750
# MWR_LIST="3000 4000 5000 6000 10000 15000 25000 50000"  #3000 5000 10000 50000 , 4000 6000 15000 
# tb_LIST="0.1 0.3"   #0.5 0.9

MN_LIST="80 120"       #  100 500 1000 , 250 750
MWR_LIST="3000 4000 5000 6000 10000 15000 25000 50000"  #3000 5000 10000 50000 , 4000 6000 15000 
tb_LIST="0.9"   #0.5 0.9



#=============================
# 1. Generate MG process if needed
#=============================
if [ ! -d "$BASE/$PROC" ]; then
  echo "Generating process directory $PROC"
  "$MG5" /home/user/WORKSPACE/RUN_CONTROLS/generate_KSDV.mg5
else
  echo "Process directory $PROC already exists, skipping generation"
fi

#=============================
# 2. Build helper binaries once
#=============================
echo "Compiling lhe2hepmc3..."
g++ -O2 -std=c++11 \
  -I"$PYTHIA_PATH/include" \
  -I"$HEPMC3_BASE/include" \
  -L"$PYTHIA_PATH/lib" \
  -L"$HEPMC3_BASE/build/outputs/lib64" \
  /home/user/WORKSPACE/RUN_CONTROLS/Manual_HEPMC3/lhe2hepmc3.cc \
  -o /home/user/WORKSPACE/RUN_CONTROLS/Manual_HEPMC3/lhe2hepmc3 \
  -lpythia8 -lHepMC3
echo "Compilation of lhe2hepmc3 done."

echo "Compiling scan_HNLDV..."
g++ -O2 -std=c++17 \
  -I"$HEPMC3_BASE/include" \
  -L"$HEPMC3_BASE/build/outputs/lib64" \
  /home/user/WORKSPACE/RUN_CONTROLS/Manual_HEPMC3/scan_HNLDV.cc \
  -o /home/user/WORKSPACE/RUN_CONTROLS/Manual_HEPMC3/scan_HNLDV \
  $(root-config --cflags --libs) \
  -lHepMC3
echo "Compilation of scan_HNLDV done."

export LD_LIBRARY_PATH="$PYTHIA_PATH/lib:$HEPMC3_BASE/build/outputs/lib64:$LD_LIBRARY_PATH"

export MG5 PYTHIA_PATH HEPMC3_BASE BASE PROC

#=============================
# 3. Per–(MN, MWR) worker
#=============================
run_point () {
  XIbeta=$1
  MN=$2
  MWR=$3

  # Analytic width and proper lifetime (mm) for the HNL
  DECAY_WIDTH=$(awk -v MN="$MN" -v MWR="$MWR" \
      'BEGIN{printf "%.12e", 1e-5*MN^5/MWR^4}')
  DECAY_LT=$(awk -v GAMA="$DECAY_WIDTH" \
      'BEGIN{printf "%.12e", 1.973e-13/GAMA}')

  WORKDIR="$BASE/${PROC}_tb${XIbeta}_MN${MN}_MWR${MWR}"
  mkdir "$WORKDIR" 2>/dev/null || {
    echo "Directory exists, skipping $WORKDIR"
    return
  }

  echo "=== Running point: MN=$MN, MWR=$MWR → $WORKDIR ==="

  # Clone template process directory
  cp -r "$BASE/$PROC/"* "$WORKDIR/"

  #-----------------------------
  # 3a. Patch param_card and run_card
  #-----------------------------
  # Zero the HNL width so SLHA does NOT override our Pythia lifetime
#   sed -i 's/^DECAY 9900012 .*/DECAY 9900012 0.0 # WN1/' \
#     "$WORKDIR/Cards/param_card.dat"

  # Update some run_card options (optional)
  sed -i '
    s/^ptj.*/ptj = 10.0/
    s/^etaj.*/etaj = -1/
    s/^event_norm.*/event_norm = sum/
    s/^lhe_version.*/lhe_version = 3.0/
  ' "$WORKDIR/Cards/run_card.dat"

  #-----------------------------
  # 3b. Pythia8 card for external conversion
  #-----------------------------
  rm -rf "$WORKDIR/Cards/pythia*"
  mkdir -p "$WORKDIR/Cards"
  cat << PYTHIA_EOF > "$WORKDIR/Cards/mypythia8_card.dat"
SLHA:useDecayTable = off
ParticleDecays:limitTau0 = on
ParticleDecays:tau0Max = 1e6
HadronLevel:Hadronize = on
HadronLevel:Decay = on
Next:numberShowEvent = 0
9900012:m0 = $MN
9900024:m0 = $MWR
9900012:mayDecay = on
9900012:isResonance = off
9900012:doExternalDecay = off
Init:showChangedParticleData = on
PYTHIA_EOF

# 9900012:mWidth = $DECAY_WIDTH    
# 9900012:tau0 = $DECAY_LT
# 9900012:onMode = off
# 9900012:addChannel = 1 0.125 0 11 -1 2
# 9900012:addChannel = 1 0.125 0 11 -3 4
# 9900012:addChannel = 1 0.125 0 13 -1 2
# 9900012:addChannel = 1 0.125 0 13 -3 4
# 9900012:addChannel = 1 0.125 0 -11 1 -2
# 9900012:addChannel = 1 0.125 0 -11 3 -4
# 9900012:addChannel = 1 0.125 0 -13 1 -2
# 9900012:addChannel = 1 0.125 0 -13 3 -4



# # Syntax: id:new = name antiName spinType chargeType colType m0 mWidth mMin mMax tau0
# 9900012:new = nuN nuNbar 2 0 0 $MN $DECAY_WIDTH $MN $MN $DECAY_LT

  #-----------------------------
  # 3c. MG5 launch: generate LHE only (no internal Pythia/Delphes)
  #-----------------------------
  cat << LAUNCH_EOF > "$WORKDIR/launch.mg5"
launch $WORKDIR
set nevents 10000
set ebeam1 6800
set ebeam2 6800
set mass 34 $MWR
set mass 9900012 $MN
set mass 9900014 100000
set mass 9900016 100000
set zetaLR = 1
set tb = $XIbeta
set alp = 0.1
set PMNSRs12 0.554
set PMNSRs23 0.1
set PMNSRs13 0.15
compute_widths 9900012
LAUNCH_EOF
# set decay 9900012 $DECAY_WIDTH

  echo "Launching MG5 for LHE generation..."
  "$MG5" "$WORKDIR/launch.mg5"

  #-----------------------------
  # 3d. Convert LHE → HepMC3 with external Pythia8
  #-----------------------------
  CARD="$WORKDIR/Cards/mypythia8_card.dat"
  LHE_GZ="$WORKDIR/Events/run_01/unweighted_events.lhe.gz"
  LHE="$WORKDIR/Events/run_01/unweighted_events.lhe"
  HEPMC="$WORKDIR/Events/run_01/events.hepmc3"

  if [ -f "$LHE_GZ" ]; then
    echo "Found gzipped LHE: $LHE_GZ"
    gunzip -c "$LHE_GZ" > "$LHE"
  elif [ ! -f "$LHE" ]; then
    echo "ERROR: No LHE file found at:"
    echo "  $LHE_GZ"
    echo "  $LHE"
    return 1
  fi

  echo "Converting to HepMC3..."
  /home/user/WORKSPACE/RUN_CONTROLS/Manual_HEPMC3/lhe2hepmc3 \
    "$CARD" \
    "$LHE" \
    "$HEPMC" | tee "$WORKDIR/Events/run_01/pythia_full.log"
  echo "HepMC3 written to: $HEPMC"

  #-----------------------------
  # 3e. Analyze HNL displacement in HepMC3
  #-----------------------------
  # SCAN_BIN=/home/user/WORKSPACE/RUN_CONTROLS/Manual_HEPMC3/scan_HNLDV
  # HISTROOT="$WORKDIR/Events/run_01/HNL_displacement.root"

  # echo "Running HNL displacement analysis..."
  # "$SCAN_BIN" "$HEPMC" "$HISTROOT" 9900012
  # echo "HNL displacement histograms saved in $HISTROOT"
}

export -f run_point

#=============================
# 4. Run over all (MN, MWR)
#=============================
# parallel --halt soon,fail=1 -j 8 run_point ::: $tb_LIST ::: $MN_LIST ::: $MWR_LIST
parallel -j 8 run_point ::: $tb_LIST ::: $MN_LIST ::: $MWR_LIST
