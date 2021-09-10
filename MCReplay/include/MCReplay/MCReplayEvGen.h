// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef MC_REPLAY_EVGEN_H
#define MC_REPLAY_EVGEN_H

#include <string>

class TVirtualMCStack;
class TFile;
class TBranch;

namespace mcreplay
{

class MCReplayEvGen
{

 public:
  MCReplayEvGen(const std::string& filename, const std::string& treename);
  MCReplayEvGen() = default;
  virtual ~MCReplayEvGen();

  // set filename and treename where to find the steps
  void setStepFilename(const std::string& filename)
  {
    mStepLoggerFilename = filename;
  }
  void setStepTreename(const std::string& treename)
  {
    mStepLoggerTreename = treename;
  }

  bool init();

  bool next(TVirtualMCStack* stack);

 private:
  // File and tree name to process
  std::string mStepLoggerFilename;
  std::string mStepLoggerTreename;

  // pointer to opened step file
  TFile* mStepFile = nullptr;

  // branches in the MCStepLogger file
  TBranch* mStepBranch = nullptr;
  TBranch* mLookupBranch = nullptr;

  // count events
  int mEventCounter = 0;
  int mEventsAvailable = 0;

  // flag if initialised
  bool mIsInitialised = false;

  ClassDefNV(MCReplayEvGen, 1);
};
} // end namespace mcreplay

#endif /* MC_REPLAY_EVGEN_H */
