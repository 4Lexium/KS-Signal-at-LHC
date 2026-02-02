Navigator:

Drell-Yan Production
> generate_KSDY.mg5 (compile processess: pp>lN via on/off shell WR)

> run_KS_scan_PRL.sh
  (was used as a general launch scrip before run_HNL_decay;
  comment out PROC (line 6), .mg5 file (line 13) !make sure the outfile path insied the .mg5 mathces your new PROC
  adjust tripple for loop S, mN, mWR !note the direcotry name copies the # for loop format

Full KS process (prompt)
> generate_KS_scan.mg5 (compile process: pp > l, N > ljj via on/off shell WR (prepare 1 day for many jobs))

> generate_HNL_decay.mg5 (compile process: p p > WR, (WR > l N1, N1 > l j j) force WR in production - faster, off-shell WR+ WL in dedcay))
 
> run_HNL_decay.sh (launch code with paramters and mass-range lists, substitute COM with tb)

Analysis:
> read_KS_CS_mWR_mNS.C (reads Drell-Yan lhe.gz output plots #sigma(m_WR ; #\srt{s}, m_N)

> read_events_scan.C (reads pp>lljj output and generateø kinematic ditributions for leptons and jets)

> read_HNL_decay_CS.C (reads the DY_CS x BR(HNL > ljj) as function of mN, mWR and tb)

Other:
> gnerate_KS.txt (extra info about generating processess)

> reweight* (ingnore, neccecary to regenerate process, kinematic changes!)
