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

#ifndef MC_REPLAY_GENERIC_STACK_H
#define MC_REPLAY_GENERIC_STACK_H

#include <TVirtualMCStack.h>

class TParticle;

namespace mcreplay
{

class MCReplayGenericStack : public TVirtualMCStack
{

 public:
  MCReplayGenericStack() = default;
  virtual ~MCReplayGenericStack();

  virtual void PushTrack(Int_t toBeDone, Int_t parent, Int_t pdg,
                         Double_t px, Double_t py, Double_t pz, Double_t e,
                         Double_t vx, Double_t vy, Double_t vz, Double_t tof,
                         Double_t polx, Double_t poly, Double_t polz,
                         TMCProcess mech, Int_t& ntr, Double_t weight,
                         Int_t is) override;

  virtual TParticle* PopNextTrack(Int_t& itrack) override
  {
    return nullptr;
  }
  virtual TParticle* PopPrimaryForTracking(Int_t i) override
  {
    // not implemented
    return nullptr;
  }

  virtual void SetCurrentTrack(Int_t trackNumber) override;

  /// Total number of tracks
  virtual Int_t GetNtrack() const override
  {
    // TODO to be implemented properly
    return mParticles.size();
  }

  /// Total number of primary tracks
  virtual Int_t GetNprimary() const override
  {
    // TODO to be implemented properly
    return mNPrimaries;
  }

  /// Current track particle
  virtual TParticle* GetCurrentTrack() const override
  {
    return mCurrentParticle;
  }

  /// Current track number
  virtual Int_t GetCurrentTrackNumber() const override
  {
    return mCurrentTrackId;
  }

  /// Number of the parent of the current track
  virtual Int_t GetCurrentParentTrackNumber() const override
  {
    // TODO to be implemented properly
    return mCurrentParentTrackId;
  }

  void newEvent();

 private:
  void clear();
  void insertParticle(TParticle* particle, int id);

 private:
  // current track ID
  Int_t mCurrentTrackId = -1;
  // current parent track ID
  Int_t mCurrentParentTrackId = -1;
  // number of primaries
  Int_t mNPrimaries = 0;
  // current TMCParticle pointer
  TParticle* mCurrentParticle = nullptr;
  // all particles ever pushed
  // TODO we can make that unique_ptrs?!
  std::vector<TParticle*> mParticles;

  ClassDefOverride(MCReplayGenericStack, 1);
};
} // end namespace mcreplay

#endif /* MC_REPLAY_GENERIC_STACK_H */
