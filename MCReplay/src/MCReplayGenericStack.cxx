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

#include <iostream>

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
  ntr = mNTracks++;
  auto particle = new TParticle(pdg, is, parent, ntr, -1, -1, px, py, pz, e, vx, vy, vz, tof);
  particle->SetPolarisation(polx, poly, polz);
  particle->SetWeight(weight);
  particle->SetUniqueID(mech);
  if (parent < 0) {
    mNPrimaries++;
  }
  insertParticle(particle, ntr);
  if (toBeDone == 1) {
    mStack.push(ntr);
  }
}

void MCReplayGenericStack::SetCurrentTrack(Int_t trackNumber)
{
  mCurrentTrackId = trackNumber;
  mCurrentParticle = mParticles[trackNumber];
  mCurrentParentTrackId = mCurrentParticle->GetFirstMother();
}

TParticle* MCReplayGenericStack::PopNextTrack(Int_t& itrack)
{
  if (!mStack.empty()) {
    itrack = mStack.top();
    mStack.pop();
    return mParticles[itrack];
  }
  itrack = -1;
  return nullptr;
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
  mNTracks = 0;

  if (!mStack.empty()) {
    std::cerr << "WARNING: There were still particles to be tracked on the stack. However, they will be removed now.\n";
  }
  while (!mStack.empty()) {
    mStack.pop();
  }
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
