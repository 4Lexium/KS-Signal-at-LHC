evgenConfig.nEventsPerJob = 20000
evgenLog.info('evgenConfig.nEventsPerJob: {0}'.format(evgenConfig.nEventsPerJob))
include('define_KS.py')
genSeq.Pythia8.useRndmGenSvc = False