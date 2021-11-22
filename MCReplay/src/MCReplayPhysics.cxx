
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

#include <cmath>
#include "MCReplay/MCReplayPhysics.h"

bool mcreplay::physics::isPhoton(int pdg)
{
  return pdg == 22;
}

bool mcreplay::physics::isElectronPositron(int pdg)
{
  return std::abs(pdg) == 11;
}

bool mcreplay::physics::isMuonAntiMuon(int pdg)
{
  return std::abs(pdg) == 13;
}

bool mcreplay::physics::isHadron(int pdg)
{
  // Taken from Pythia8 ParticleDataEntry::isHadron

  if (pdg <= 100 || (pdg >= 1000000 && pdg <= 9000000) || pdg >= 9900000) {
    return false;
  }
  if (pdg == 130 || pdg == 310) {
    return true;
  }
  if (pdg % 10 == 0 || (pdg / 10) % 10 == 0 || (pdg / 100) % 10 == 0) {
    return false;
  }
  return true;
}
