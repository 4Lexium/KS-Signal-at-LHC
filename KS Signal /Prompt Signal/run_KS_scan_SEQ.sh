set -e

MG5=/storage/alexanum/madgraph/MG5_aMC_v2_9_25/bin/mg5_aMC
BASE=/storage/alexanum/WORKBENCH/runs/KS_WR

$MG5 /home/alexanum/WORKSPACE/testing/MadGraph5/scripts/generate_KS_scan.mg5

for MN in 500
do
  for MWR in 4000
  do
    cat << EOF > launch_tmp.mg5
launch $BASE -n run_MN${MN}_MWR${MWR}
set nevents 10000
set ebeam1 7000
set ebeam2 7000

set mass 9900012 $MN
set mass 34 $MWR
EOF

    $MG5 launch_tmp.mg5
  done
done
