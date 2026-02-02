set -e
# Define directories
MG5=/storage/alexanum/madgraph/MG5_aMC_v2_9_25/bin/mg5_aMC
BASE=/storage/alexanum/WORKBENCH/runs
# PROC=KS_WR 
PROC=HNL_decay2

# Recommend Independantly running: [ /storage/alexanum/madgraph/MG5_aMC_v2_9_25/bin/mg5_aMC /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/generate_KS_scan.mg5 ]
# rm -rf KS_WR_*
if [ ! -d "$BASE/$PROC" ]; then
  echo "Generating process directory $PROC"
  $MG5 /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/generate_HNL_decay.mg5
  # $MG5 /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/generate_KSDY.mg5
else
  echo "Process directory $PROC already exists, skipping generation"
fi

MN_LIST="200 250 500"
MWR_LIST="3000 5000"
tb_list="0 0.1 0.5 0.9"
# LHCS_LIST="16500"
export MG5 BASE PROC

run_point () {
  # LHCS=$1
  MN=$1
  MWR=$2
  tanbeta=$3
  # WORKDIR=$BASE/${PROC}_S${LHCS}_MN${MN}_MWR${MWR}
  WORKDIR=$BASE/${PROC}_tb${tanbeta}_MN${MN}_MWR${MWR}
  cp -r $BASE/$PROC $WORKDIR

  cat << EOF > $WORKDIR/launch.mg5
launch $WORKDIR
set nevents 10000
set ebeam1 6500
set ebeam2 6500
set nb_core 8
set mass 9900012 $MN
set mass 9900014 100000
set mass 9900016 100000
set mass 34 $MWR
set decay 9900012 0.0000001*$MN**5/$MWR**4
set zetaLR = 1
set tb = $tanbeta
set alp = 0.1
EOF

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

# parallel -j 8 run_point ::: $LHCS_LIST ::: $MN_LIST ::: $MWR_LIST
parallel -j 14 run_point ::: $MN_LIST ::: $MWR_LIST ::: $tb_list
# parallerism 2 cores per job (saturates run-time bonus)
# madgraph 2 x 8 = 16 

# For countour analysis: use one event set but reweight:
#   reweighting  X no!
#   interpolation
#   adaptive grids