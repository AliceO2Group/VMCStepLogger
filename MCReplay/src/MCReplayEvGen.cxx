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

#include <vector>
#include <iostream>

#include <TMCProcess.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TVirtualMCStack.h>

#include "MCStepLogger/StepInfo.h"
#include "MCReplay/MCReplayEvGen.h"

ClassImp(mcreplay::MCReplayEvGen);

using namespace mcreplay;

MCReplayEvGen::MCReplayEvGen(const std::string& filename, const std::string& treename)
  : mStepLoggerFilename{filename}, mStepLoggerTreename{treename}
{
}

MCReplayEvGen::~MCReplayEvGen()
{
  if (mStepFile) {
    mStepFile->Close();
  }
}

bool MCReplayEvGen::init()
{
  mStepFile = TFile::Open(mStepLoggerFilename.c_str(), "READ");
  auto tree = (TTree*)mStepFile->Get(mStepLoggerTreename.c_str());
  if (!tree) {
    std::cerr << "Cannot retrieve tree " << mStepLoggerTreename << " in file " << mStepLoggerFilename << "\n";
    return false;
  }
  mStepBranch = tree->GetBranch("Steps");
  mLookupBranch = tree->GetBranch("Lookups");
  if (!mStepBranch || !mLookupBranch) {
    std::cerr << "Cannot get branches \"Steps\" and \"Lookups\"\n";
    return false;
  }
  mEventsAvailable = mStepBranch->GetEntries();
  mIsInitialised = true;
  return true;
}

bool MCReplayEvGen::next(TVirtualMCStack* stack)
{
  if (!mIsInitialised) {
    std::cerr << "Not yet initialised\n";
    return false;
  }
  if (mEventCounter >= mEventsAvailable) {
    std::cerr << "Ran out of events, only " << mEventsAvailable << " available\n";
    return false;
  }

  // flag if we found at least one primary
  bool foundPrimary{false};

  // just provide for pushing to stack, nothing else
  int stackTrackID;
  std::vector<o2::StepInfo>* steps = nullptr;
  o2::StepLookups* lookups = nullptr;

  mStepBranch->SetAddress(&steps);
  mLookupBranch->SetAddress(&lookups);

  mStepBranch->GetEntry(mEventCounter);
  mLookupBranch->GetEntry(mEventCounter);

  for (const auto& step : *steps) {
    if (!step.newtrack || lookups->tracktoparent[step.trackID] >= 0) {
      continue;
    }
    std::cout << "Push primary " << step.trackID << " with PDG " << lookups->tracktopdg[step.trackID] << " to stack" << std::endl;
    foundPrimary = true;
    stack->PushTrack(1, -1, lookups->tracktopdg[step.trackID], step.px, step.py, step.pz, step.E, step.x, step.y, step.z, step.t, 1., 1., 1., TMCProcess(step.prodprocess), stackTrackID, 1., 0);
  }
  // No longer needed
  delete steps;
  delete lookups;
  mEventCounter++;
  return foundPrimary;
}
