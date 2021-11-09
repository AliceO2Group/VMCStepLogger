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
#include "MCReplay/MCReplayEvGen.h"
#include "MCReplay/MCReplayGenericApplication.h"

ClassImp(mcreplay::MCReplayGenericApplication);

using namespace mcreplay;

MCReplayGenericApplication::MCReplayGenericApplication(const std::string& geoFilename, const std::string& geoKeyname)
  : TVirtualMCApplication{"MCReplayGenericApplication", "MCReplayGenericApplication"},
    mGeoFilename{geoFilename},
    mGeoKeyname{geoKeyname},
    mGeoManager{nullptr},
    mStack{nullptr},
    mPrimGen{nullptr}
{
  if (mGeoFilename.empty()) {
    // we only require the file name to be non-empty. In that case the geo manager would take the first key
    // under which a geometry can be found
    ::Fatal("MCReplayGenericApplication::ctor", "Need file and key name where to find TGeo geometry");
  }

  mGeoManager = TGeoManager::Import(mGeoFilename.c_str(), mGeoKeyname.c_str());

  if (!mGeoManager) {
    ::Fatal("MCReplayGenericApplication::ctor", "Could not load geometry from path %s under key %s", mGeoFilename.c_str(), mGeoKeyname.c_str());
  }
}

void MCReplayGenericApplication::ConstructGeometry()
{
  ::Info("MCReplayGenericApplication::ConstructGeometry", "geometry already constructed");
}

void MCReplayGenericApplication::BeginEvent()
{
  mStack->newEvent();
}

void MCReplayGenericApplication::GeneratePrimaries()
{
  if (!mPrimGen->next(mStack)) {
    ::Warning("MCReplayGenericApplication::GeneratePrimaries", "Could not retrieve primaries");
  }
}
