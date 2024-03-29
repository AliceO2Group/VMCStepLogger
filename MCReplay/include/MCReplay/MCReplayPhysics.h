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

#ifndef MC_REPLAY_PHYSICS_H
#define MC_REPLAY_PHYSICS_H

#include <array>
#include <algorithm>
#include <cstring>

// Fix process and cut names
namespace mcreplay
{
namespace physics
{
// standard GEANT process parameter names
constexpr std::array<const char*, 15> namesProcesses = {"PAIR", "COMP", "PHOT", "PFIS", "DRAY", "ANNI", "BREM", "HADR", "MUNU", "DCAY", "LOSS", "MULS", "CKOV", "RAYL", "LABS"};
// standard GEANT cut parameter names (1-11)
// additional Replay cut parameter names (12-n):
// CUTALL: energy cut applying to all tracks independent of PDG
constexpr std::array<const char*, 12> namesCuts = {"CUTGAM", "CUTELE", "CUTNEU", "CUTHAD", "CUTMUO", "BCUTE", "BCUTM", "DCUTE", "DCUTM", "PPCUTM", "TOFMAX", "CUTALL"};
constexpr const char* unknownParam = "UNKNOWN";

template <typename T, std::size_t N>
int paramToIndex(const std::array<T, N>& allParams, const char* paramName)
{
  const auto& it = std::find_if(allParams.begin(), allParams.end(), [&paramName](const char* comp) { return std::strcmp(comp, paramName) == 0; });
  if (it == allParams.end()) {
    return -1;
  }
  return it - allParams.begin();
}

template <typename T, std::size_t N>
const char* indexToParam(const std::array<T, N>& allParams, std::size_t index)
{
  if (index >= allParams.size()) {
    return physics::unknownParam;
  }
  return allParams[index];
}

bool isPhoton(int pdg);

bool isElectronPositron(int pdg);

bool isMuonAntiMuon(int pdg);

bool isHadron(int pdg);

} // namespace physics
} // end namespace mcreplay

#endif /* MC_REPLAY_PHYSICS_H */
