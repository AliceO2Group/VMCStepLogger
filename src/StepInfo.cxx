// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//  @file   StepInfo.cxx
//  @author Sandro Wenzel
//  @since  2017-06-29
//  @brief  structures encapsulating information about MC stepping

#include "MCStepLogger/StepInfo.h"
#include <TArrayI.h>
#include <TParticle.h>
#include <TVirtualMC.h>
#include <chrono>

#include <TDatabasePDG.h>
#include <TGeoManager.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TMCProcess.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

ClassImp(o2::StepInfo);
ClassImp(o2::MagCallInfo);

namespace o2
{
// construct directly using virtual mc
StepInfo::StepInfo(TVirtualMC* mc)
{
  assert(mc);

  // init base time point
  if (stepcounter == -1) {
    starttime = std::chrono::high_resolution_clock::now();
  }
  stepcounter++;
  stepid = stepcounter;

  auto stack = mc->GetStack();

  trackID = stack->GetCurrentTrackNumber();
  lookupstructures.insertPDG(trackID, mc->TrackPid());

  auto id = mc->CurrentVolID(copyNo);
  volId = id;

  auto parentID = trackID < stack->GetNprimary() ? -1 : stack->GetCurrentParentTrackNumber();
  lookupstructures.insertParent(trackID, parentID);

  // try to resolve the module via external map
  // keep information in faster vector once looked up
  auto volname = mc->CurrentVolName();
  lookupstructures.insertVolName(volId, volname);

  if (volnametomodulemap && volnametomodulemap->size() > 0 && volId >= 0) {
    if (lookupstructures.getModuleAt(volId) == nullptr) {
      // lookup in map
      auto iter = volnametomodulemap->find(volname);
      if (iter != volnametomodulemap->end()) {
        lookupstructures.insertModuleName(volId, iter->second);
      }
      else {
        // std::cout << "VOL NOT FOUND .. GO UP UNTIL WE FIND A KNOWN VOLUME NAME " << volname << "\n";
        // trying to look upward
        int up = 1;
        int upcopy;
        int upVolID;
        const char* upVolName;
        int limit = 20;
        do {
          upVolID = mc->CurrentVolOffID(up, upcopy);
          upVolName = mc->CurrentVolOffName(up);
          up++;
          auto iter2 = volnametomodulemap->find(upVolName);
          if (iter2 != volnametomodulemap->end()) {
            lookupstructures.insertModuleName(volId, iter2->second);
            // std::cout << "FIXING TO " << iter2->second;
            break;
          }
        } while (up < limit);
        //
      }
    }
  }

  // sensitive volume info
  if (lookupstructures.volidtoissensitive.size() > 0) {
    insensitiveRegion = lookupstructures.volidtoissensitive[volId];
  }

  double xd, yd, zd;
  mc->TrackPosition(xd, yd, zd);
  x = xd;
  y = yd;
  z = zd;
  step = mc->TrackStep();
  maxstep = mc->MaxStep();
  E = mc->Etot();
  auto now = std::chrono::high_resolution_clock::now();
  // cputimestamp = std::chrono::duration_cast<std::chrono::nanoseconds>(now - starttime).count();
  nsecondaries = mc->NSecondaries();

  if (nsecondaries > 0) {
    lookupstructures.setProducedSecondary(trackID, true);
  }
  
  if (mc->IsTrackExiting()) {
    lookupstructures.setCrossedBoundary(trackID, true);
  }
  
  prodprocess = stack->GetCurrentTrack()->GetUniqueID();

  TArrayI procs;
  mc->StepProcesses(procs);
  nprocessesactive = procs.GetSize();

  // was track stopped due to energy limit ??
  stopped = mc->IsTrackStop();
  exited = mc->IsTrackExiting();
  entered = mc->IsTrackEntering();
  newtrack = mc->IsNewTrack();
}

const char* StepInfo::getProdProcessAsString() const {
  return TMCProcessName[prodprocess];
}
  
std::chrono::time_point<std::chrono::high_resolution_clock> StepInfo::starttime;
int StepInfo::stepcounter = -1;
std::map<std::string, std::string>* StepInfo::volnametomodulemap = nullptr;
std::vector<std::string*> StepInfo::volidtomodulevector;
StepLookups StepInfo::lookupstructures;

MagCallInfo::MagCallInfo(TVirtualMC* mc, float ax, float ay, float az, float aBx, float aBy, float aBz)
  : x{ ax }, y{ ay }, z{ az }, B{ std::sqrt(aBx * aBx + aBy * aBy + aBz * aBz) }
{
  stepcounter++;
  id = stepcounter;
  stepid = StepInfo::stepcounter;
}

int MagCallInfo::stepcounter = -1;

// try to init a map of volID to sensitive/ornot
// by using a list of sensitive volume names (given in a file)
// the current ROOT geometry loaded
bool StepLookups::initSensitiveVolLookup(const std::string& filename)
{
  if (!gGeoManager) {
    std::cerr << "[MCSTEPLOG] : Cannot setup sensitive lookup since GeoManager not found \n";
    return false;
  }
  // get list of all TGeoVolumes
  auto vlist = gGeoManager->GetListOfVolumes();
  volidtoissensitive.resize(vlist->GetEntries(), false);

  // open for reading or fail
  std::ifstream ifs;
  ifs.open(filename);
  if (ifs.is_open()) {
    std::string line;
    std::vector<int> ids;
    while (std::getline(ifs, line)) {
      std::istringstream ss(line);
      std::string token;
      // split the line into key + value
      int counter = 0;
      std::string keyvalue[2] = { "NULL", "NULL" };
      while (counter < 2 && std::getline(ss, token, ':')) {
        if (!token.empty()) {
          keyvalue[counter] = token;
          counter++;
        }
      }
      // mark this as sensitive
      int index = atoi(keyvalue[0].c_str());
      assert(index >= 0);
      assert(index < volidtoissensitive.size());
      volidtoissensitive[index] = true;
    }
    return true;
  }
  return false;
}

} // namespace o2
