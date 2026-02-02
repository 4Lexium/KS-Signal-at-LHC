set -e
# Define directories
MG5=/storage/alexanum/madgraph/MG5_aMC_v2_9_25/bin/mg5_aMC
PYTHIA_PATH=/storage/alexanum/madgraph/MG5_aMC_v2_9_25/HEPTools/pythia8
# MG5_aMC> install mg5amc_py8_interface_path
BASE=/storage/alexanum/WORKBENCH/runs
# PROC=KS_WR 
PROC=HNL_pythia_decay

# Recommend Independantly running: [ /storage/alexanum/madgraph/MG5_aMC_v2_9_25/bin/mg5_aMC /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/generate_KS_scan.mg5 ]
# rm -rf KS_WR_*
if [ ! -d "$BASE/$PROC" ]; then
  echo "Generating process directory $PROC"
  $MG5 /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/generate_HNL_decay.mg5
  # $MG5 /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/generate_KSDY.mg5
else
  echo "Process directory $PROC already exists, skipping generation"
fi

MN_LIST="100 500"
MWR_LIST="3000"
tb_list="0.1 0.5"
# LHCS_LIST="16500"
export MG5 PYTHIA_PATH BASE PROC

run_point () {
  # LHCS=$1
  MN=$1
  MWR=$2
  tanbeta=$3
  # WORKDIR=$BASE/${PROC}_S${LHCS}_MN${MN}_MWR${MWR}
  WORKDIR=$BASE/${PROC}_tb${tanbeta}_MN${MN}_MWR${MWR}
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
ParticleDecays:limitTau0 = on
ParticleDecays:tau0Max = 10.0
Check:epTolErr = 1e-6
BeamRemnants:primordialKThard = 1.8
ColourReconnection:range = 1.8
PYTHIA_EOF
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

  # Calculate Numerical value for HNL width: Gamma(mN, mWR)
  DECAY_WIDTH=$(awk -v MN=$MN -v MWR=$MWR 'BEGIN{printf "%.6e", 1e-7*MN^5/MWR^4}')

  # Launch script with Pythia8
  cat << LAUNCH_EOF > $WORKDIR/launch.mg5
launch $WORKDIR
set nevents 10000
set ebeam1 6500
set ebeam2 6500
set mass 9900012 $MN
set mass 9900014 100000
set mass 9900016 100000
set mass 34 $MWR
set decay 9900012 $DECAY_WIDTH
set zetaLR 1
set tb $tanbeta
set alp 0.1
set shower Pythia8
set pythia8_path $PYTHIA_PATH
LAUNCH_EOF

  $MG5 $WORKDIR/launch.mg5
}
# set tb = 1 (1 max mixing, 0 no gauge mixing)
# set alp = 0 (exp(ialp)=1 no damping, alp close to pi/2 max damping) 

# set PMNSRs12 1.0
# set PMNSRs23 0.0
# set PMNSRs13 0.0
# set PMNSRdel 0.0
# set PMNSRphi1 0.0
# set PMNSRphi2 0.0
# set PMNSRphi3 0.0

export -f run_point

parallel -j 4 run_point ::: $MN_LIST ::: $MWR_LIST ::: $tb_list