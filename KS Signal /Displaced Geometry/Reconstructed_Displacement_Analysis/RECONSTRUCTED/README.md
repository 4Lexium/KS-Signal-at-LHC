dump_Rec_LLP
>
> reads truth-level BSM branches and compares to reconstructed track-DCA values
>
> plots comprison to analytical Poisson distribution and correlation with transverse momenta

scanDV_testBKG
>
> Comparing the DCA method for displaced signal model and background MCs 

scanDV_testSV
>
> testing several DCA modifications: using PV coordinate correction, and SV instead of perigree. 

Execute each python script: [ python XXX.py ] in enviroment that has: Numnpy, awkward, math, ROOT, uproot.
F.ex. using VENV: 

> python3 -m venv $HOME/venvs/uproot_env
> 
> source $HOME/venvs/uproot_env/bin/activate
> 
> pip install --upgrade pip
> 
> pip install uproot awkward numpy matplotlib
> 
> pip install --force-reinstall --ignore-installed numpy==2.0.2
> 
> python -c "import uproot; print(uproot.__version__)"
> 
> python -c "import uproot, awkward as ak; print('uproot:', uproot.__version__, '| awkward:', ak.__version__)"
> 
> deactivate (to exit)
