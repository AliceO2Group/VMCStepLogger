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

#include "TMCReplay/GenericStack.h"
#include "TMCReplay/GenericApplication.h"

ClassImp(tmcreplay::GenericApplication);

using namespace tmcreplay;

GenericApplication::GenericApplication(const std::string& geoFilename, const std::string& geoKeyname)
  : TVirtualMCApplication{"GenericApplication", "GenericApplication"},
    fGeoFilename{geoFilename},
    fGeoKeyname{geoKeyname},
    fGeoManager{nullptr}
{
}

void GenericApplication::ConstructGeometry()
{
  if (!gGeoManager) {
    // construct one if not present
    new TGeoManager();
  }

  fGeoManager = gGeoManager;

  if (!fGeoFilename.empty() && !fGeoKeyname.empty()) {
    // Otherwise assume that the geometry was constructed already by the user
    fGeoManager->Import(fGeoFilename.c_str(), fGeoKeyname.c_str());
  }
}

void GenericApplication::BeginEvent()
{
  fStack->newEvent();
}
