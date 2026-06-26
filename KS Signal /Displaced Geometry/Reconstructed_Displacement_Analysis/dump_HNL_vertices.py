# dump_HNL_vertices.py
# Reads an EVNT file, loops over HepMC truth, and prints
# production & decay vertices for HNLs (pdgId = 9900012).
'''
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
setupATLAS -c el9
asetup AthGeneration,23.6.63
which athena.py
athena.py /srv/EVENT_analysis/dump_HNL_vertices.py
'''

import AthenaPython.PyAthena as PyAthena
from AthenaPython.PyAthena import StatusCode

from AthenaCommon.AppMgr import ServiceMgr as svcMgr
from AthenaCommon.AlgSequence import AlgSequence

# ---------------------------------------------------------------------------
# CONFIGURE INPUT EVNT FILE
# ---------------------------------------------------------------------------
# evnt_file = "/srv/ROOT_files/mc23_13p6TeV.NEWDVTEST.MGPy8EG_LRSM_WR_HNL_5000_100_2L2J_MadSpin.EVNT.pool.root"  
# evnt_file = "/srv/ROOT_files/mc23_13p6TeV.2NEWDVTEST.MGPy8EG_LRSM_WR_HNL_10000_50_2L2J_MadSpin.EVNT.pool.root"
evnt_file = "/srv/ROOT_files/KS_LLP_MN50_MWR5000.EVNT.pool.root"
# evnt_file = "/srv/mc23_13p6TeV.EXTRAPARAMS3.MGPy8EG_LRSM_WR_HNL_5000_100_2L2J_MadSpin.EVNT.pool.root"
# just importing this sets up POOL reading in this old-style config
import AthenaPoolCnvSvc.ReadAthenaPool

svcMgr.EventSelector.InputCollections = [evnt_file]

# ---------------------------------------------------------------------------
# DEFINE PYTHON ALGORITHM
# ---------------------------------------------------------------------------

class HNLVertexDump(PyAthena.Alg):
    def __init__(self, name="HNLVertexDump", **kw):
        super(HNLVertexDump, self).__init__(name, **kw)
        self.mc_container = "GEN_EVENT"
        self.target_pdg = 9900012  # HNL PDG ID

    def initialize(self):
        self.msg.info("Initializing HNLVertexDump")
        self.msg.info("  MC container: %s", self.mc_container)
        self.msg.info("  target PDG: %d", self.target_pdg)
        return StatusCode.Success

    def execute(self):
        # Retrieve HepMC McEventCollection
        mcColl = self.evtStore[self.mc_container]

        # Simple event counter
        if not hasattr(self, "_evt"):
            self._evt = 0
        self.msg.info("----- Event %d -----", self._evt)

        for ge in mcColl:  # HepMC::GenEvent
            for p in ge.particles():  # HepMC::GenParticle
                if abs(p.pdg_id()) != self.target_pdg:
                    continue

                vProd = p.production_vertex()
                vEnd  = p.end_vertex()
                if vProd is None or vEnd is None:
                    continue

                prodPos = vProd.position()  # HepMC::FourVector
                endPos  = vEnd.position()

                xp, yp, zp = prodPos.x(), prodPos.y(), prodPos.z()
                xd, yd, zd = endPos.x(),  endPos.y(),  endPos.z()

                dx = xd - xp
                dy = yd - yp
                dz = zd - zp

                Lxy = (dx*dx + dy*dy) ** 0.5
                L   = (dx*dx + dy*dy + dz*dz) ** 0.5

                self.msg.info(
                    "HNL pdgId=%d "
                    "prod=(%.3f, %.3f, %.3f) mm  "
                    "decay=(%.3f, %.3f, %.3f) mm  "
                    "Lxy=%.3f mm  L=%.3f mm",
                    p.pdg_id(), xp, yp, zp, xd, yd, zd, Lxy, L
                )

        self._evt += 1
        return StatusCode.Success

# ---------------------------------------------------------------------------
# ADD ALGORITHM TO SEQUENCE
# ---------------------------------------------------------------------------

algSeq = AlgSequence()
algSeq += HNLVertexDump("HNLVertexDump")

# Limit number of events if you want:
from AthenaCommon.AppMgr import theApp
theApp.EvtMax = 100  # change to -1 for all events
