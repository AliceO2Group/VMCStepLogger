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

#include "TGeoManager.h"

#include "MCReplay/MCReplayGenericStack.h"
#include "MCReplay/MCReplayGenericApplication.h"

ClassImp(mcreplay::MCReplayGenericApplication);

using namespace mcreplay;

MCReplayGenericApplication::MCReplayGenericApplication(const std::string& geoFilename, const std::string& geoKeyname, const std::string& stepFilename, const std::string& stepTreename)
  : TVirtualMCApplication{"MCReplayGenericApplication", "MCReplayGenericApplication"},
    mGeoFilename{geoFilename},
    mGeoKeyname{geoKeyname},
    mGeoManager{nullptr},
    mStack{nullptr},
    mPrimGen{stepFilename, stepTreename}
{
}

void MCReplayGenericApplication::ConstructGeometry()
{
  if (!gGeoManager) {
    // construct one if not present
    new TGeoManager();
  }

  mGeoManager = gGeoManager;

  if (!mGeoFilename.empty() && !mGeoKeyname.empty()) {
    // Otherwise assume that the geometry was constructed already by the user
    mGeoManager->Import(mGeoFilename.c_str(), mGeoKeyname.c_str());
  }
}

void MCReplayGenericApplication::BeginEvent()
{
  mStack->newEvent();
}

void MCReplayGenericApplication::GeneratePrimaries()
{
  if (!mPrimGen.next(mStack)) {
    ::Warning("MCReplayGenericApplication::GeneratePrimaries", "Could not retrieve primaries");
  }
}
