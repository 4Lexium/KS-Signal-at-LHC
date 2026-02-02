PROC_DIR="KS_WR_MN500_MWR3000"
RUN="run_01"

MWR_LIST=(3000 3500 4000 4500 5000 5000 6000)
MN_LIST=(100 200 300 400 500 600 700 800 900 1000)

RWGT_CARD="${PROC_DIR}/Cards/reweight_card.dat"

if [ ! -d "$PROC_DIR" ]; then
  echo "ERROR: Process directory not found: $PROC_DIR"
  exit 1
fi

if [ ! -f "${PROC_DIR}/Events/${RUN}/unweighted_events.lhe.gz" ]; then
  echo "ERROR: LHE file not found"
  exit 1
fi

echo "change rwgt_dir rwgt_MWR_MN_scan" > "$RWGT_CARD"

for MWR in "${MWR_LIST[@]}"; do
  for MN in "${MN_LIST[@]}"; do

    RWGT_NAME="MWR${MWR}_MN${MN}"
    cat >> "$RWGT_CARD" <<EOF

launch --rwgt_name=${RWGT_NAME}
  set mass 34 ${MWR}
  set mass 9900012 ${MN}

EOF

  done
done

echo "Reweight card written to $RWGT_CARD"


# Run MadGraph reweighting
cd "$PROC_DIR" || exit 1
./bin/madevent reweight ${RUN}