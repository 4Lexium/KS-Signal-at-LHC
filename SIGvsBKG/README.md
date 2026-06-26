multiplicity_scan.py
> reconstruction analysis of multiplicities (#jets, #leptons) and testing the general Preselection procedure (validation with dijet-system)

Comparion of Signal Models, Background MC and data from 2022 period J
>
> 



To run any of the C macro:

> navigate to directory containing Metadata JSON, DAOD ROOT Files (presumably /storage/...) you can setup the enivoment there, and run you C located at /home/XXX.C
> 
> > export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
>
> > alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
> 
> > setupATLAS -c centos7
> 
> > asetup Athena,24.0.11
> 
> > root -l -q -b /home/XXX.C
> >

Extracting the Metadata from DAOD ROOT files:

> /Metadata
