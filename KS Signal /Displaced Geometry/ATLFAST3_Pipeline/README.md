Create a directory as shown in directory.png. To avoid spam-files call all execustions from inside the temp folder (it is safe to clean)

/jobConfig contains the the script that declares the run_cards for MadGraph and Pythia, and also: manuall update of the Lhe file to add HNL lifetime (see Lifetime_formula.png for proof of principle). Original idea from: https://gitlab.cern.ch/atlas-physics/pmg/mcjoboptions/-/blob/master/537xxx/537565/MadGraphControl_HSS_Run3.py

The procedure for lhe > hepmc > EVNT > HITS > AOD > DAOD_LLP is better explained in detail in (todo)

(Make sure the files are called correctly i.e. if /jobconfig is on the same level as /temp use ../define_KS.py when executing from inside /temp. We also recomend when entering the /srv after launching an ATLAS platform to save output files (especcially HITS and AOD) in /storage or some other location with large disk space.  


