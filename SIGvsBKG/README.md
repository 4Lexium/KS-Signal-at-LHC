multiplicity_scan.py
> reconstruction analysis of multiplicities (#jets, #leptons) and testing the general Preselection procedure (validation with dijet-system)

Comparion of Signal Models, Background MC and data from 2022 period J
>
> see Samples_INFO for infomation on samples
>
> reco22.C
>
> > Signal (Model B) vs Backgroud samples
>
> data22.C
>
> > Signal vs Background vs data22 (data preprocessed using GoodEventsList, see internal helper-functions)  
>
> recoDV.py
>
> > Signal (Model A) vs Backgroud samples (using displacement branches)
>
> sgn22.C
>
> > Signal (Model A + B) vs Background samples (computing and comparing significance) 

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
> 

Extracting the Metadata from DAOD ROOT files:

> /Metadata
