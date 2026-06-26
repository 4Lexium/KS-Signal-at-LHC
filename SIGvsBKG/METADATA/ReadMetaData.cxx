#include <iostream>
#include <memory>

#include <TFile.h>

// xAOD includes
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "xAODRootAccess/Init.h"
#include "xAODCore/tools/Message.h"

#include "xAODMetaData/FileMetaData.h"
#include "xAODMetaData/CutBookkeeper.h"
#include "xAODMetaData/CutBookkeeperContainer.h"

// For convenience
using std::cout;
using std::endl;

int main( int argc, char* argv[] ) {

  if ( argc < 2 ) {
    std::cerr << "Usage: " << argv[0] << " <DAOD.root>" << std::endl;
    return 1;
  }

  const char* fileName = argv[1];

  // Initialize the xAOD infrastructure
  xAOD::Init().ignore();

  // Open ROOT file
  std::unique_ptr<TFile> ifile( TFile::Open( fileName, "READ" ) );
  if ( !ifile || ifile->IsZombie() ) {
    std::cerr << "ERROR: Could not open file " << fileName << std::endl;
    return 1;
  }

  // TEvent to access metadata and events
  xAOD::TEvent event( xAOD::TEvent::kBranchAccess );

  // Connect the file to TEvent
  if ( !event.readFrom( ifile.get() ).isSuccess() ) {
    std::cerr << "ERROR: Failed to read from file " << fileName << std::endl;
    return 1;
  }

  cout << "File: " << fileName << endl;

  // ----------------------------------------------------------------------
  // 1) Access xAOD::FileMetaData
  // ----------------------------------------------------------------------
  const xAOD::FileMetaData* fmd = nullptr;
  if ( event.retrieveMetaInput( fmd, "FileMetaData" ).isSuccess() ) {

    cout << "=== FileMetaData ===" << endl;

    // Example properties (check which keys exist in your release)
    std::string val_str;
    int val_int = 0;

    if ( fmd->value( xAOD::FileMetaData::projectName, val_str ) ) {
      cout << "  projectName       : " << val_str << endl;
    }
    if ( fmd->value( xAOD::FileMetaData::productionRelease, val_str ) ) {
      cout << "  productionRelease : " << val_str << endl;
    }
    if ( fmd->value( xAOD::FileMetaData::amgTag, val_str ) ) {
      cout << "  AMI/AMG tag       : " << val_str << endl;
    }
    if ( fmd->value( xAOD::FileMetaData::dataType, val_int ) ) {
      cout << "  dataType          : " << val_int << " (0: data, 1: MC)" << endl;
    }
    if ( fmd->value( xAOD::FileMetaData::runNumber, val_int ) ) {
      cout << "  runNumber         : " << val_int << endl;
    }

  } else {
    std::cerr << "WARNING: Could not retrieve FileMetaData" << std::endl;
  }

  // ----------------------------------------------------------------------
  // 2) Access CutBookkeepers (for cross sections / event counts)
  // ----------------------------------------------------------------------
  const xAOD::CutBookkeeperContainer* allCBK = nullptr;
  if ( event.retrieveMetaInput( allCBK, "CutBookkeepers" ).isSuccess() ) {

    cout << "=== CutBookkeepers ===" << endl;
    cout << "  size = " << allCBK->size() << endl;

    for ( const xAOD::CutBookkeeper* cbk : *allCBK ) {

      const std::string& name = cbk->name();
      const std::string& desc = cbk->description();
      const int cycle        = cbk->cycle();
      const double nAccepted = cbk->nAcceptedEvents();
      const double kFactor   = cbk->kFactor();
      const double xSec      = cbk->crossSection();
      const double genEff    = cbk->efficiency();

      cout << "  CBK: " << name
           << " (cycle " << cycle << ")" << endl;
      cout << "    description   : " << desc << endl;
      cout << "    nAccepted     : " << nAccepted << endl;
      cout << "    crossSection  : " << xSec << " pb" << endl;
      cout << "    kFactor       : " << kFactor << endl;
      cout << "    efficiency    : " << genEff << endl;
    }

    // Often you want the "MC generator" (all events) cut-bookkeeper.
    // A common pattern is to select the one with name "AllExecutedEvents"
    // and the highest cycle:
    const xAOD::CutBookkeeper* cbkAll = nullptr;
    int maxCycle = -1;
    for ( const xAOD::CutBookkeeper* cbk : *allCBK ) {
      if ( cbk->name() == "AllExecutedEvents" && cbk->cycle() > maxCycle ) {
        maxCycle = cbk->cycle();
        cbkAll   = cbk;
      }
    }
    if ( cbkAll ) {
      cout << "  Selected CBK \"AllExecutedEvents\" (cycle "
           << cbkAll->cycle() << "):" << endl;
      cout << "    nAccepted     : " << cbkAll->nAcceptedEvents() << endl;
      cout << "    sumOfEventWts : " << cbkAll->sumOfEventWeights() << endl;
      cout << "    xSection      : " << cbkAll->crossSection() << " pb" << endl;
    }

  } else {
    std::cerr << "WARNING: Could not retrieve CutBookkeepers" << std::endl;
  }

  // ----------------------------------------------------------------------
  // 3) (Optional) Access event-level metadata like EventFormat
  // ----------------------------------------------------------------------
  const xAOD::EventFormat* evtFmt = event.eventFormat();
  if ( evtFmt ) {
    cout << "=== EventFormat ===" << endl;
    cout << "  number of branches: " << evtFmt->size() << endl;
  }

  cout << "Done." << endl;
  return 0;
}
