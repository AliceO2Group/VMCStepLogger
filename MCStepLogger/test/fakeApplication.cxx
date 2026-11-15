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

// A stand-in for FairMCApplication that counts how often its original methods are reached

#include <TVirtualMCApplication.h>

class FairMCApplication : public TVirtualMCApplication
{
 public:
  FairMCApplication() : TVirtualMCApplication("fake", "fake application") {}
  void ConstructGeometry() override;
  void InitGeometry() override {}
  void GeneratePrimaries() override {}
  void BeginEvent() override {}
  void BeginPrimary() override {}
  void PreTrack() override {}
  void Stepping() override;
  void PostTrack() override {}
  void FinishPrimary() override {}
  void FinishEvent() override;
};

namespace
{
int nConstructGeometry = 0;
int nStepping = 0;
int nFinishEvent = 0;
} // namespace

void FairMCApplication::ConstructGeometry() { ++nConstructGeometry; }
void FairMCApplication::Stepping() { ++nStepping; }
void FairMCApplication::FinishEvent() { ++nFinishEvent; }

extern "C" TVirtualMCApplication* makeFakeApplication() { return new FairMCApplication(); }

extern "C" void fakeApplicationCounts(int* construct, int* stepping, int* finish)
{
  *construct = nConstructGeometry;
  *stepping = nStepping;
  *finish = nFinishEvent;
}
