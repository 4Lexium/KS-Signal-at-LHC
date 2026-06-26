> KS Signal Topology
>  - Prompt Region
>  - Displaced Region
>  - Merged Region

> Analysis Pipeline
>  - HNL production (Drell-Yan): generated with MadGraph
>  - HNL decay (prompt/displaced): MadGprah+Pythia+Delphes
>  - Combined: generation with MadGraph+Pythia -> simulation with Altfast3 (DAOD/DAOD_LLP)

> Model
>  - Complete Left-Right symmetric model
>  - acess: https://cp3.irmp.ucl.ac.be/projects/feynrules/wiki/LRSM_NLO#no1
>  - UFO: mlrsm-1.1/UFOs/mlrsm-4f

> Requirements
>  - Python 3.9.12
>  - MadGraph 2.9.25 - 3.7 (acessible: http://madgraph.phys.ucl.ac.be/cgi-bin/Register/download.cgi)
>  - Pythia8 8.3.17 
>  - Delphes 3
>  - (Pythia and Delphes can be installed inside MadGraph directory /HEPTools, then set correct path using bin> display options)
>  - ROOT 6.26/04
>  - GNU Paralell 20260222 (acess: https://zenodo.org/records/18735643)
