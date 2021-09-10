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

#include <TParticle.h>

#include "MCReplay/MCReplayGenericStack.h"

ClassImp(mcreplay::MCReplayGenericStack);

using namespace mcreplay;

MCReplayGenericStack::~MCReplayGenericStack()
{
  clear();
}

void MCReplayGenericStack::PushTrack(Int_t toBeDone, Int_t parent, Int_t pdg,
                                     Double_t px, Double_t py, Double_t pz, Double_t e,
                                     Double_t vx, Double_t vy, Double_t vz, Double_t tof,
                                     Double_t polx, Double_t poly, Double_t polz,
                                     TMCProcess mech, Int_t& ntr, Double_t weight,
                                     Int_t is)
{
  // We are not setting the track number but use the same indexing used before.
  // Hence, the passed value of ntr is assumed to be already the correct ID
  // This is important to ensure the mapping of track IDs to volumes and produced number of secondaries etc. is the same in the original simulation and the replay.
  // Otherwise, for instance a steplogging of the replay would have slight differences in the lookup.
  auto particle = new TParticle(pdg, is, parent, ntr, -1, -1, px, py, pz, e, vx, vy, vz, tof);
  particle->SetPolarisation(polx, poly, polz);
  particle->SetWeight(weight);
  particle->SetUniqueID(mech);
  if (parent < 0) {
    mNPrimaries++;
  }
  insertParticle(particle, ntr);
}

void MCReplayGenericStack::SetCurrentTrack(Int_t trackNumber)
{
  mCurrentTrackId = trackNumber;
  mCurrentParticle = mParticles[trackNumber];
  mCurrentParentTrackId = mCurrentParticle->GetFirstMother();
}

void MCReplayGenericStack::clear()
{
  for (auto& p : mParticles) {
    delete p;
  }
  mParticles.clear();
  mCurrentTrackId = -1;
  mCurrentParentTrackId = -1;
  mNPrimaries = 0;
}

void MCReplayGenericStack::newEvent()
{
  clear();
}

void MCReplayGenericStack::insertParticle(TParticle* particle, int id)
{
  if (mParticles.size() <= id) {
    mParticles.resize(id + 1, nullptr);
  }
  mParticles[id] = particle;
}
