Convert EVNT -> TRUTH3:  (OR TO RUN SIG_vs_BKG)
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase  
    alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh' 
    setupATLAS -c centos7 
    asetup Athena,24.0.11 
    which Derivation_tf.py
(check with ls for EVNTFILE)
Derivation_tf.py --CA --inputEVNTFile=../ROOT_files/KS_LLP_MN70_MWR5000.EVNT.pool.root --outputDAODFile KS_LLP_MN70_MWR5000.EVNT.pool.root --formats TRUTH3  

(move it manually inside the ROOT_files folder)

Run EVENT_analysis:
(Navigate inside ATLAS_jobdir)
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
    alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
    setupATLAS -c el9
    asetup AthGeneration,23.6.63
    which athena.py
    (move to temp)
    athena.py /srv/EVENT_analysis/dump_HNL_vertices.py
