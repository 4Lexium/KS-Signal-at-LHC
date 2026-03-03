set -e
# Define directories
MG5=/storage/alexanum/madgraph/MG5_aMC_v2_9_25/bin/mg5_aMC
PYTHIA_PATH=/storage/alexanum/madgraph/MG5_aMC_v2_9_25/HEPTools/pythia8
DELPHES_PATH=/storage/alexanum/madgraph/MG5_aMC_v2_9_25/Delphes
# MG5_aMC> install mg5amc_py8_interface_path
BASE=/storage/alexanum/WORKBENCH/runs
# PROC=KS_WR 
# PROC=KS_delphes_decay_2
PROC=KS_delphes_decay_PythGEO



# Recommend Independantly running: [ /storage/alexanum/madgraph/MG5_aMC_v2_9_25/bin/mg5_aMC /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/generate_KS_scan.mg5 ]
# rm -rf KS_WR_*
if [ ! -d "$BASE/$PROC" ]; then
  echo "Generating process directory $PROC"
  # $MG5 /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/remastered/generate_KS_chained.mg5
  $MG5 /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/remastered/generate_KS_DV.mg5
else
  echo "Process directory $PROC already exists, skipping generation"
fi

# MN_LIST="50 100 500 1000 1500"
# MWR_LIST="3000 4000 5000 6000 7000"
MN_LIST="100"
MWR_LIST="10000"
# tb_list="0.1 0.5"
# LHCS_LIST="16500"
export MG5 PYTHIA_PATH DELPHES_PATH BASE PROC

run_point () {
  # LHCS=$1
  MN=$1
  MWR=$2
  # Calculate Numerical value for HNL width: Gamma(mN, mWR)
  DECAY_WIDTH=$(awk -v MN=$MN -v MWR=$MWR 'BEGIN{printf "%.12e", 1e-5*MN^5/MWR^4}')
  # tanbeta=$3
  WORKDIR=$BASE/${PROC}_MN${MN}_MWR${MWR}
  rm -rf $WORKDIR
  # WORKDIR=$BASE/${PROC}_tb${tanbeta}_MN${MN}_MWR${MWR}
  cp -r $BASE/$PROC $WORKDIR
  
  # Create pythia8 card for custom settings
  mkdir -p $WORKDIR/Cards
  cat << PYTHIA_EOF > $WORKDIR/Cards/pythia8_card.dat
Tune:pp = 14
PartonLevel:MPI = on
PartonLevel:ISR = on
PartonLevel:FSR = on
HadronLevel:Hadronize = on
HadronLevel:Decay = on
SLHA:useDecayTable = on
9900012:mayDecay = on
ParticleDecays:limitTau0 = on
ParticleDecays:tau0Max = 1000000
Check:epTolErr = 1e-6
BeamRemnants:primordialKThard = 1.8
ColourReconnection:range = 1.8
PYTHIA_EOF
# 9900012:mWidth = $DECAY_WIDTH
# 9900012:mayDecay = on
# 9900012:isResonance = false
# Init:showChangedParticleData = on
# 9900012:show = on
# 9900012:tau0 = 0 manual override to PROMPT decay 
# SpaceShower:pTmaxMatch = 2
# TimeShower:pTmaxMatch = 2
# JetMatching:merge = off
# PDF:pSet = LHAPDF6:NNPDF23_lo_as_0130_qed

  # Patch run_card.dat
  sed -i '
    s/^ptj.*/ptj = 10.0/
    s/^etaj.*/etaj = -1/
    s/^event_norm.*/event_norm = sum/
    s/^lhe_version.*/lhe_version = 3.0/
  ' "$WORKDIR/Cards/run_card.dat"


  # Copy Delphes card to working directory (ensures MadGraph finds it)
  cp "$DELPHES_PATH/cards/delphes_card_ATLAS.tcl" "$WORKDIR/Cards/delphes_card.dat"

# Patch the Lepton Isolation (enter the tcl card and do the changes)
# # Electron isolation
# sed -i '/module Isolation ElectronIsolation {/,/}/ s/^\s*set DeltaRMax .*/  set DeltaRMax 0.2/' $WORKDIR/Cards/delphes_card.dat
# sed -i '/module Isolation ElectronIsolation {/,/}/ s/^\s*set PTRatioMax .*/  set PTRatioMax 0.8/' $WORKDIR/Cards/delphes_card.dat
# # Muon isolation
# sed -i '/module Isolation MuonIsolation {/,/}/ s/^\s*set DeltaRMax .*/  set DeltaRMax 0.2/' $WORKDIR/Cards/delphes_card.dat
# sed -i '/module Isolation MuonIsolation {/,/}/ s/^\s*set PTRatioMax .*/  set PTRatioMax 0.8/' $WORKDIR/Cards/delphes_card.dat

# Append (>>) HNL-friendly settings to the ATLAS card
#==================================================================
# Max proper lifetime (mm/c) and Save vertex info and extend tracker volume
# set TrackSmearing:ApplyToPileUp false    
# set ParticleDecays:tau0Max 1e6         
# set Propagator:MaxZ 10000              
# set Propagator:MaxR 5000           
cat << DELPHES_EOF >> $WORKDIR/Cards/delphes_card.dat
set ParticleDecays:tau0Max 1e6          
set Propagator:MaxZ 10000             
set Propagator:MaxR 5000            
set ParticlePropagator:UseLifetime true
set ParticlePropagator:UseProperTime true
set ParticlePropagator:SaveDisplacedVertex true   
set ParticleDecays:tau0Max 1e6
set TrackSmearing:ApplyToPileUp false   
DELPHES_EOF
    
  # Create launch script with the local card path
  # Launch script with Pythia8
  cat << LAUNCH_EOF > $WORKDIR/launch.mg5
launch $WORKDIR
set nevents 10000
set ebeam1 6800
set ebeam2 6800
set mass 9900012 $MN
set WN1 $DECAY_WIDTH
set decay 9900012 $DECAY_WIDTH
set time_of_flight 1
set mass 9900014 100000
set mass 9900016 100000
set mass 34 $MWR
set zetaLR = 1
set tb = 0.1
set alp = 0.1
set PMNSRs12 0.554
set PMNSRs23 0.1
set PMNSRs13 0.15
compute_widths 9900012
set shower Pythia8
set pythia8_path $PYTHIA_PATH
set detector Delphes
set delphes_executable DelphesHepMC
set delphes_path $DELPHES_PATH
set delphes_card delphes_card.dat
LAUNCH_EOF

  $MG5 $WORKDIR/launch.mg5
}
# nano MG5_aMC_v2_9_25/input/mg5_configuration.txt
# nano /storage/alexanum/WORKBENCH/runs/KS_delphes_decay/bin/internal/run_delphes3
# grep delphes /storage/alexanum/madgraph/MG5_aMC_v2_9_25/input/mg5_configuration.txt

# set decay 9900012 $DECAY_WIDTH
# set time_of_flight 1  (compute TOF for LHE)
export -f run_point

parallel -j 4 run_point ::: $MN_LIST ::: $MWR_LIST




# Made a symlink DelphesHepMC  -> DelphesHepMC3
# ln -sf DelphesHepMC3 DelphesHepMC
# chmod +x DelphesHepMC3
# (check) ls -l DelphesHepMC



# when erros with load .h
# nano /Delphes/Delphes.Env.sh
# added: export ROOT_INCLUDE_PATH="${DIR}/external:${ROOT_INCLUDE_PATH}"
# (reload) source /storage/alexanum/madgraph/MG5_aMC_v2_9_25/Delphes/DelphesEnv.sh 
# (check) echo $ROOT_INCLUDE_PATH
# root -l
# gSystem->Load("libDelphes.so");