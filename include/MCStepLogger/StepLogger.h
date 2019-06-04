// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


#ifndef STEPLOGGER_H_
#define STEPLOGGER_H_

#include "MCStepLogger/StepInfo.h"
#include "MCStepLogger/MetaInfo.h"
#include "MCStepLogger/StepLoggerUtilities.h"
#include <TBranch.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <TVirtualMCApplication.h>
#include <iostream>



namespace o2
{

class StepLogger
{
  int stepcounter = 0;

  std::set<int> trackset;
  std::set<int> pdgset;
  std::map<int, int> volumetosteps;
  std::map<int, int> volumetohits;
  std::map<int, std::string> idtovolname;
  std::map<int, int> volumetoNSecondaries;            // number of secondaries created in this volume
  std::map<std::pair<int, int>, int> volumetoProcess; // mapping of volumeid x processID to secondaries produced

  std::vector<StepInfo> container;
  bool mTTreeIO = false;
  bool mVolMapInitialized = false;
  // To add a hit to the volume a step was made, set in addStep and used if mTTreeIO == false
  int mCurrentVolId = -1;

 public:

  StepLogger()
  {
    // check if streaming or interactive
    // configuration done via env variable
    if (std::getenv("MCSTEPLOG_TTREE")) {
      mTTreeIO = true;
    }
  }

  static StepLogger& Instance()
  {
    static StepLogger logger;
    return logger;
  }

  void addStep(TVirtualMC* mc)
  {
    if (!mVolMapInitialized) {
      // try to load the volumename -> modulename mapping
      mcsteploggerutilities::initVolumeMap();
      StepInfo::lookupstructures.initSensitiveVolLookup(mcsteploggerutilities::getSensitiveVolFile());
      mVolMapInitialized = true;
    }

    if (mTTreeIO) {
      container.emplace_back(mc);
    } else {
      assert(mc);
      stepcounter++;

      auto stack = mc->GetStack();
      assert(stack);
      trackset.insert(stack->GetCurrentTrackNumber());
      pdgset.insert(mc->TrackPid());
      int copyNo;
      mCurrentVolId = mc->CurrentVolID(copyNo);

      TArrayI procs;
      mc->StepProcesses(procs);

      if (volumetosteps.find(mCurrentVolId) == volumetosteps.end()) {
        volumetosteps.insert(std::pair<int, int>(mCurrentVolId, 0));
        // At the same time extend hit map
        volumetohits.insert(std::pair<int, int>(mCurrentVolId, 0));
      } else {
        volumetosteps[mCurrentVolId]++;
      }
      if (idtovolname.find(mCurrentVolId) == idtovolname.end()) {
        idtovolname.insert(std::pair<int, std::string>(mCurrentVolId, std::string(mc->CurrentVolName())));
      }

      // for the secondaries
      if (volumetoNSecondaries.find(mCurrentVolId) == volumetoNSecondaries.end()) {
        volumetoNSecondaries.insert(std::pair<int, int>(mCurrentVolId, mc->NSecondaries()));
      } else {
        volumetoNSecondaries[mCurrentVolId] += mc->NSecondaries();
      }

      // for the processes
      for (int i = 0; i < mc->NSecondaries(); ++i) {
        auto process = mc->ProdProcess(i);
        auto p = std::pair<int, int>(mCurrentVolId, process);
        if (volumetoProcess.find(p) == volumetoProcess.end()) {
          volumetoProcess.insert(std::pair<std::pair<int, int>, int>(p, 1));
        } else {
          volumetoProcess[p]++;
        }
      }
    }
  }

  void addHit(int detId, int stepId = -1)
  {
    if(container.size() == 0) {
      std::cerr << "Empty step container, exit...\n";
      exit(1);
    }
    if (mTTreeIO) {
      if(stepId < 0) {
        container.back().detectorHitId = detId;
      } else if(stepId < static_cast<int>(container.size())) {
        container[stepId].detectorHitId = detId;
      } else {
        std::cerr << "Step ID " << stepId << " could not be found. Skip hit\n";
      }
    } else if(mCurrentVolId >= 0) {
      volumetohits[mCurrentVolId]++;
    }
  }

  void clear()
  {
    stepcounter = 0;
    trackset.clear();
    pdgset.clear();
    volumetosteps.clear();
    volumetohits.clear();
    idtovolname.clear();
    volumetoNSecondaries.clear();
    volumetoProcess.clear();
    mCurrentVolId = -1;
    if (mTTreeIO) {
      container.clear();
    }
    StepInfo::resetCounter();
  }

  // prints list of processes for volumeID
  void printProcesses(int volid)
  {
    for (auto& p : volumetoProcess) {
      if (p.first.first == volid) {
        std::cerr << "P[" << TMCProcessName[p.first.second] << "]:" << p.second << "\t";
      }
    }
  }

  void flush()
  {
    if (!mTTreeIO) {
      std::cerr << "[STEPLOGGER]: did " << stepcounter << " steps \n";
      std::cerr << "[STEPLOGGER]: transported " << trackset.size() << " different tracks \n";
      std::cerr << "[STEPLOGGER]: transported " << pdgset.size() << " different types \n";
      // summarize steps per volume
      for (auto& p : volumetosteps) {
        std::cerr << "[STEPLOGGER]: VolName " << idtovolname[p.first] << " COUNT steps " << p.second << " SECONDARIES "
                  << volumetoNSecondaries[p.first] << " COUNT hist " << volumetohits[p.first] << " ";
        // loop over processes
        printProcesses(p.first);
        std::cerr << "\n";
      }
      std::cerr << "[STEPLOGGER]: ----- END OF EVENT ------\n";
    } else {
      mcsteploggerutilities::flushToTTree("Steps", &container);
      mcsteploggerutilities::flushToTTree("Lookups", &StepInfo::lookupstructures);
      // we need to reset some parts of the lookupstructures for the next event
      StepInfo::lookupstructures.clearTrackLookups();
    }
    clear();
  }
};

} // namespace o2

#endif /* STEPLOGGER_H_ */
