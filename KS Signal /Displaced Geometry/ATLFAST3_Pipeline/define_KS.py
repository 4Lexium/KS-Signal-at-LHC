import os
from AthenaCommon.Include import include
from MadGraphControl.MadGraphUtils import *
keepOutput = False # Debug
def StringToFloat(s):
  if "p" in s:
    return float(s.replace("p", "."))
  return float(s)

# Block for extracting paramters from jobname (parallell job processing)
from MCJobOptionUtils.JOsupport import get_physics_short
# joName = get_physics_short()
# jobConfigParts = joName.split("_")
# masses['9900012'] = StringToFloat(jobConfigParts[4])
# masses['34'] = StringToFloat(jobConfigParts[5])

#-----------------------------------------------------------------------------------------------------
# Define Directory for LRSM paramesters including: ('mass', 'lrsminputs', 'pmnsblock')
#-----------------------------------------------------------------------------------------------------
params = {}

param_blocks = {}
decays = {}
masses = {}

# Define MWR and MN1
MWR = 10e3
MN = 50
width_N = 1e-5*MN**5/MWR**4
ctau_N = 1.973e-13/width_N #mm

masses['9900012'] = MN
masses['34'] = MWR
masses['9900014'] = 100e3
masses['9900016'] = 100e3

# Define General settings
gridpack_mode=False
# Merging settings
maxjetflavor=5
ickkw=0
nJetMax=2
ktdurham=30
dparameter=0.4
tofoption= -1

#if masses['9900012']/4 > 30:
#    ktdurham = masses['9900012']/4
#processline = "pp>lljj"
# evgenLog.info('ktdurham set to %i' %ktdurham)    
#process=processline
#evgenConfig.process = processline
#evgenLog.info('process for jet merging set to %s' %processline)

#-----------------------------------------------------------------------------------------------------
# Define MG proc for the KS signal
#-----------------------------------------------------------------------------------------------------
gen_process = '''
import model /srv/jobConfiguration/mlrsm-4f
define l+ = e+ mu+
define l- = e- mu-
define p = g u c d s u~ c~ d~ s~
define j = g u c d s u~ c~ d~ s~ 
define l = e+ e- mu+ mu- ta+ ta-
generate p p > WR+ > l+ N1, N1 > l j j 
add process p p > WR- > l- N1, N1 > l j j
output -f
'''

evt_multiplier = 1
'''
if evt_multiplier>0:
    if runArgs.nEventsPerJob>0:
        nevents=runArgs.nEventsPerJob*evt_multiplier
    else:    
        nevents=500*evt_multiplier
'''
nevents = 20000

evgenLog.info( "ControlFile: Final process card " + gen_process )

run_settings = {'lhe_version':'3.0',
                'pdlabel'    : "'lhapdf'",
                'lhaid'      : 260000,
                'ickkw'      : '0',
                'ktdurham'   : ktdurham,
                'maxjetflavor':4, # 5 flavor scheme
                'asrwgtflavor':4, # 5 flavor scheme
                'xptb':0,
                'etal':10,
                'ptl':0,
                'ptj':20.,
                'drjj':0.0,
                'xqcut':0.,
                'nhel' : 0,   #helicity
                'bwcutoff': 15.0,
                'time_of_flight': tofoption
                }

print(nevents)
run_settings['use_syst']='F'
run_settings['nevents'] = nevents

# Set up the process
process_dir = new_process(gen_process)
# Set up the run card
modify_run_card(process_dir=process_dir,runArgs=runArgs,settings=run_settings)

params['mass'] = masses

params['lrsminputs'] = {
    '1': 0.1,        # tb
    '2': 0.1,        # alp
    '7': 1.0         # zetalr
}

params['pmnsblock'] = {
    '11': 0.554,    #pmnsrs12 etc.
    '12': 0.1,
    '13': 0.15
}

modify_param_card(process_dir=process_dir,params=params)

# Generate the events
generate(process_dir=process_dir,runArgs=runArgs)


#-----------------------------------------------------------------------------------------------------
# MANUAL LIFETIME 
#-----------------------------------------------------------------------------------------------------
import random
random.seed(runArgs.randomSeed)

# Calculate lifetime function
def lifetime(avgtau = ctau_N):       
    import math
    t = random.random()
    return -1.0 * avgtau * math.log(t)
print (process_dir)

# Replace HNL's (9900012) lifetime-column inside the lhe file (manually)
unzip1 = subprocess.Popen(['gunzip',process_dir+'/Events/run_01/unweighted_events.lhe.gz'])
unzip1.wait()
oldlhe = open(process_dir+'/Events/run_01/unweighted_events.lhe','r')
newlhe = open(process_dir+'/Events/run_01/unweighted_events2.lhe','w')
init = True
for line in oldlhe:
    if init==True:
        newlhe.write(line)
        if '</init>' in line:
            init = False
    else:  
        if 'vent' in line or line.startswith("<"):
            newlhe.write(line)
            continue
        newline = line.rstrip('\n')
        columns = (' '.join(newline.split())).split()
        pdgid = int(columns[0])
        if pdgid == 9900012:
            part1 = line[:-22]
            part2 = "%.11E" % (lifetime(avgtau = ctau_N))
            part3 = line[-12:]
            newlhe.write(part1+part2+part3)
        else:
            newlhe.write(line)

oldlhe.close()
newlhe.close()
zip1 = subprocess.Popen(['gzip',process_dir+'/Events/run_01/unweighted_events2.lhe'])
zip1.wait()
shutil.move(process_dir+'/Events/run_01/unweighted_events2.lhe.gz',
            process_dir+'/Events/run_01/unweighted_events.lhe.gz')
os.remove(process_dir+'/Events/run_01/unweighted_events.lhe')



# Remember to set saveProcDir to FALSE before sending for production!!
arrange_output(process_dir=process_dir,runArgs=runArgs,lhe_version=3,saveProcDir=True)

# information about this generation
evgenLog.info('Registered generation of WR->HNL + l -> W->jj ll production, grid point '+' decoded into mass point ' + str(masses['9900012']) + ' ' + str(masses['34']))   #str(jobConfigParts)+

include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
include("Pythia8_i/Pythia8_MadGraph.py")

# Reset the number of processes for Pythia8
check_reset_proc_number(opts)

# Pythia options

pythia = genSeq.Pythia8

# genSeq.Pythia8.Commands += ["9900012:tau0 = " + str(1)]
genSeq.Pythia8.Commands += ["9900012:tau0 = " + str(lifetime(avgtau = ctau_N))]
testSeq.TestHepMC.MaxTransVtxDisp = 100000000 #in mm
testSeq.TestHepMC.MaxVtxDisp = 100000000 #in mm
testSeq.TestHepMC.MaxNonG4Energy = 100000000 #in MeV

# pythia.Commands += ["24:mMin = 0.2"]
# pythia.Commands += ["24:onMode = off"]                # switch off all W decays
# pythia.Commands += ["24:onIfAny = 11 12 13 14 15 16"] # switch on W->lnu
# #genSeq.Pythia8.Commands += [
# #    "Init:showAllParticleData = on",
# #    "Next:numberShowLHA = 10",
# #    "Next:numberShowEvent = 10",
# #    "Merging:process = guess",
#     #"9900012:spinType = 1",
# #]

#pythia.Commands += ["showChargedParticleData = on"]
'''
if nJetMax>0:
    pythia.Commands += ["Merging:mayRemoveDecayProducts = on",
                                "Merging:doKTMerging = on",
                                "Merging:ktType = 1",
                                "Merging:nJetMax = %i"%nJetMax,
                                "Merging:Process = guess",#%process,
                                "Merging:TMS = %f"%run_settings['ktdurham'],
                                "Merging:nQuarksMerge = %s" %run_settings['maxjetflavor'],
                                "Merging:Dparameter = %f"%dparameter
                                ]

# configuration for EvgenJobTransforms
evgenLog.info("Registered generation of WR -> HNL")
evgenConfig.contact = ["eirik.gramstad@cern.ch"]
evgenConfig.keywords += ["BSM", "2lepton","exotic"]
evgenConfig.description = "LRSM"
evgenConfig.generators += ["EvtGen"]

testSeq.TestHepMC.MaxTransVtxDisp = 200000 #in mm
testSeq.TestHepMC.MaxVtxDisp = 1000000 #in mm
'''

evgenLog.info("Registered generation of WR -> HNL")
evgenConfig.contact = ["eirik.gramstad@cern.ch"]
evgenConfig.keywords += ["BSM", "2lepton","exotic"]
evgenConfig.description = "LRSM"
evgenConfig.generators += ["EvtGen"]

#--------------------------------------------------------------
genSeq.Pythia8.computeEfficiency = True
