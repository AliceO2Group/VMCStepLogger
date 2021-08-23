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

#include <algorithm>
#include <iostream>

// For now include all TGeo headers here
#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoTrd1.h>
#include <TGeoTrd2.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoSphere.h>
#include <TGeoPara.h>
#include <TGeoPgon.h>
#include <TGeoPcon.h>
#include <TGeoEltu.h>
#include <TGeoHype.h>
#include <TGeoArb8.h>

#include "TMCReplay/Engine.h"

ClassImp(tmcreplay::Engine);

namespace vmcsl = o2::mcstepanalysis;
using namespace tmcreplay;

Engine::Engine(const std::string& filename, const std::string& treename)
  : TVirtualMC{"ReplayEngine", "ReplayEngine", kTRUE}, fIsRootGeometrySupported{kTRUE}, fIsRunStopped{kFALSE}, fIsEventStopped{kFALSE}, fIsTrackStopped{kFALSE}, fFilename{filename}, fTreename{treename}, fCurrentStepInfo{nullptr}, fCurrentMagCallInfo{nullptr}, fCurrentLookups{nullptr}, fCurrentEvent{0}, fCurrentStep{nullptr}, fCurrentTrackLength{0.}, fProcessesGlobal(physics::namesProcesses.size(), -1), fCutsGlobal(physics::namesCuts.size(), -1.), fcurrentProcesses{nullptr}, fcurrentCuts{nullptr}, fGeoManager{nullptr}
{
}

Engine::~Engine()
{
  for (auto& m : fProcesses) {
    delete m;
  }
  for (auto& m : fCuts) {
    delete m;
  }
}

void Engine::Init()
{
  fApplication->ConstructGeometry();
  fApplication->MisalignGeometry();
  fApplication->InitGeometry();
  // gGeoManager must be valid pointer now
  fGeoManager = gGeoManager;
}

Bool_t Engine::GetTransformation(const TString& volumePath, TGeoHMatrix& matrix)
{
  // We have to preserve the modeler state
  fGeoManager->PushPath();
  if (!fGeoManager->cd(volumePath.Data())) {
    fGeoManager->PopPath();
    return kFALSE;
  }
  matrix = *fGeoManager->GetCurrentMatrix();
  fGeoManager->PopPath();
  return kTRUE;
}

Bool_t Engine::GetShape(const TString& volumePath, TString& shapeType, TArrayD& par)
{
  Int_t npar;
  fGeoManager->PushPath();
  if (!fGeoManager->cd(volumePath.Data())) {
    fGeoManager->PopPath();
    return kFALSE;
  }
  auto vol = fGeoManager->GetCurrentVolume();
  fGeoManager->PopPath();
  if (!vol) {
    return kFALSE;
  }
  auto shape = vol->GetShape();
  TClass* class_type = shape->IsA();
  if (class_type == TGeoBBox::Class()) {
    shapeType = "BOX";
    npar = 3;
    par.Set(npar);
    auto box = (TGeoBBox*)shape;
    par.AddAt(box->GetDX(), 0);
    par.AddAt(box->GetDY(), 1);
    par.AddAt(box->GetDZ(), 2);
    return kTRUE;
  }
  if (class_type == TGeoTrd1::Class()) {
    shapeType = "TRD1";
    npar = 4;
    par.Set(npar);
    auto trd1 = (TGeoTrd1*)shape;
    par.AddAt(trd1->GetDx1(), 0);
    par.AddAt(trd1->GetDx2(), 1);
    par.AddAt(trd1->GetDy(), 2);
    par.AddAt(trd1->GetDz(), 3);
    return kTRUE;
  }
  if (class_type == TGeoTrd2::Class()) {
    shapeType = "TRD2";
    npar = 5;
    par.Set(npar);
    auto trd2 = (TGeoTrd2*)shape;
    par.AddAt(trd2->GetDx1(), 0);
    par.AddAt(trd2->GetDx2(), 1);
    par.AddAt(trd2->GetDy1(), 2);
    par.AddAt(trd2->GetDy2(), 3);
    par.AddAt(trd2->GetDz(), 4);
    return kTRUE;
  }
  if (class_type == TGeoTrap::Class()) {
    shapeType = "TRAP";
    npar = 11;
    par.Set(npar);
    auto trap = (TGeoTrap*)shape;
    Double_t tth = TMath::Tan(trap->GetTheta() * TMath::DegToRad());
    par.AddAt(trap->GetDz(), 0);
    par.AddAt(tth * TMath::Cos(trap->GetPhi() * TMath::DegToRad()), 1);
    par.AddAt(tth * TMath::Sin(trap->GetPhi() * TMath::DegToRad()), 2);
    par.AddAt(trap->GetH1(), 3);
    par.AddAt(trap->GetBl1(), 4);
    par.AddAt(trap->GetTl1(), 5);
    par.AddAt(TMath::Tan(trap->GetAlpha1() * TMath::DegToRad()), 6);
    par.AddAt(trap->GetH2(), 7);
    par.AddAt(trap->GetBl2(), 8);
    par.AddAt(trap->GetTl2(), 9);
    par.AddAt(TMath::Tan(trap->GetAlpha2() * TMath::DegToRad()), 10);
    return kTRUE;
  }
  if (class_type == TGeoTube::Class()) {
    shapeType = "TUBE";
    npar = 3;
    par.Set(npar);
    auto tube = (TGeoTube*)shape;
    par.AddAt(tube->GetRmin(), 0);
    par.AddAt(tube->GetRmax(), 1);
    par.AddAt(tube->GetDz(), 2);
    return kTRUE;
  }
  if (class_type == TGeoTubeSeg::Class()) {
    shapeType = "TUBS";
    npar = 5;
    par.Set(npar);
    auto tubs = (TGeoTubeSeg*)shape;
    par.AddAt(tubs->GetRmin(), 0);
    par.AddAt(tubs->GetRmax(), 1);
    par.AddAt(tubs->GetDz(), 2);
    par.AddAt(tubs->GetPhi1(), 3);
    par.AddAt(tubs->GetPhi2(), 4);
    return kTRUE;
  }
  if (class_type == TGeoCone::Class()) {
    shapeType = "CONE";
    npar = 5;
    par.Set(npar);
    auto cone = (TGeoCone*)shape;
    par.AddAt(cone->GetDz(), 0);
    par.AddAt(cone->GetRmin1(), 1);
    par.AddAt(cone->GetRmax1(), 2);
    par.AddAt(cone->GetRmin2(), 3);
    par.AddAt(cone->GetRmax2(), 4);
    return kTRUE;
  }
  if (class_type == TGeoConeSeg::Class()) {
    shapeType = "CONS";
    npar = 7;
    par.Set(npar);
    auto cons = (TGeoConeSeg*)shape;
    par.AddAt(cons->GetDz(), 0);
    par.AddAt(cons->GetRmin1(), 1);
    par.AddAt(cons->GetRmax1(), 2);
    par.AddAt(cons->GetRmin2(), 3);
    par.AddAt(cons->GetRmax2(), 4);
    par.AddAt(cons->GetPhi1(), 5);
    par.AddAt(cons->GetPhi2(), 6);
    return kTRUE;
  }
  if (class_type == TGeoSphere::Class()) {
    shapeType = "SPHE";
    npar = 6;
    par.Set(npar);
    auto sphe = (TGeoSphere*)shape;
    par.AddAt(sphe->GetRmin(), 0);
    par.AddAt(sphe->GetRmax(), 1);
    par.AddAt(sphe->GetTheta1(), 2);
    par.AddAt(sphe->GetTheta2(), 3);
    par.AddAt(sphe->GetPhi1(), 4);
    par.AddAt(sphe->GetPhi2(), 5);
    return kTRUE;
  }
  if (class_type == TGeoPara::Class()) {
    shapeType = "PARA";
    npar = 6;
    par.Set(npar);
    auto para = (TGeoPara*)shape;
    par.AddAt(para->GetX(), 0);
    par.AddAt(para->GetY(), 1);
    par.AddAt(para->GetZ(), 2);
    par.AddAt(para->GetTxy(), 3);
    par.AddAt(para->GetTxz(), 4);
    par.AddAt(para->GetTyz(), 5);
    return kTRUE;
  }
  if (class_type == TGeoPgon::Class()) {
    shapeType = "PGON";
    auto pgon = (TGeoPgon*)shape;
    Int_t nz = pgon->GetNz();
    const Double_t* rmin = pgon->GetRmin();
    const Double_t* rmax = pgon->GetRmax();
    const Double_t* z = pgon->GetZ();
    npar = 4 + 3 * nz;
    par.Set(npar);
    par.AddAt(pgon->GetPhi1(), 0);
    par.AddAt(pgon->GetDphi(), 1);
    par.AddAt(pgon->GetNedges(), 2);
    par.AddAt(pgon->GetNz(), 3);
    for (Int_t i = 0; i < nz; i++) {
      par.AddAt(z[i], 4 + 3 * i);
      par.AddAt(rmin[i], 4 + 3 * i + 1);
      par.AddAt(rmax[i], 4 + 3 * i + 2);
    }
    return kTRUE;
  }
  if (class_type == TGeoPcon::Class()) {
    shapeType = "PCON";
    auto pcon = (TGeoPcon*)shape;
    Int_t nz = pcon->GetNz();
    const Double_t* rmin = pcon->GetRmin();
    const Double_t* rmax = pcon->GetRmax();
    const Double_t* z = pcon->GetZ();
    npar = 3 + 3 * nz;
    par.Set(npar);
    par.AddAt(pcon->GetPhi1(), 0);
    par.AddAt(pcon->GetDphi(), 1);
    par.AddAt(pcon->GetNz(), 2);
    for (Int_t i = 0; i < nz; i++) {
      par.AddAt(z[i], 3 + 3 * i);
      par.AddAt(rmin[i], 3 + 3 * i + 1);
      par.AddAt(rmax[i], 3 + 3 * i + 2);
    }
    return kTRUE;
  }
  if (class_type == TGeoEltu::Class()) {
    shapeType = "ELTU";
    npar = 3;
    par.Set(npar);
    auto eltu = (TGeoEltu*)shape;
    par.AddAt(eltu->GetA(), 0);
    par.AddAt(eltu->GetB(), 1);
    par.AddAt(eltu->GetDz(), 2);
    return kTRUE;
  }
  if (class_type == TGeoHype::Class()) {
    shapeType = "HYPE";
    npar = 5;
    par.Set(npar);
    auto hype = (TGeoHype*)shape;
    par.AddAt(TMath::Sqrt(hype->RadiusHypeSq(0., kTRUE)), 0);
    par.AddAt(TMath::Sqrt(hype->RadiusHypeSq(0., kFALSE)), 1);
    par.AddAt(hype->GetDZ(), 2);
    par.AddAt(hype->GetStIn(), 3);
    par.AddAt(hype->GetStOut(), 4);
    return kTRUE;
  }
  if (class_type == TGeoGtra::Class()) {
    shapeType = "GTRA";
    npar = 12;
    par.Set(npar);
    auto trap = (TGeoGtra*)shape;
    Double_t tth = TMath::Tan(trap->GetTheta() * TMath::DegToRad());
    par.AddAt(trap->GetDz(), 0);
    par.AddAt(tth * TMath::Cos(trap->GetPhi() * TMath::DegToRad()), 1);
    par.AddAt(tth * TMath::Sin(trap->GetPhi() * TMath::DegToRad()), 2);
    par.AddAt(trap->GetH1(), 3);
    par.AddAt(trap->GetBl1(), 4);
    par.AddAt(trap->GetTl1(), 5);
    par.AddAt(TMath::Tan(trap->GetAlpha1() * TMath::DegToRad()), 6);
    par.AddAt(trap->GetH2(), 7);
    par.AddAt(trap->GetBl2(), 8);
    par.AddAt(trap->GetTl2(), 9);
    par.AddAt(TMath::Tan(trap->GetAlpha2() * TMath::DegToRad()), 10);
    par.AddAt(trap->GetTwistAngle(), 11);
    return kTRUE;
  }
  if (class_type == TGeoCtub::Class()) {
    shapeType = "CTUB";
    npar = 11;
    par.Set(npar);
    auto ctub = (TGeoCtub*)shape;
    const Double_t* lx = ctub->GetNlow();
    const Double_t* tx = ctub->GetNhigh();
    par.AddAt(ctub->GetRmin(), 0);
    par.AddAt(ctub->GetRmax(), 1);
    par.AddAt(ctub->GetDz(), 2);
    par.AddAt(ctub->GetPhi1(), 3);
    par.AddAt(ctub->GetPhi2(), 4);
    par.AddAt(lx[0], 5);
    par.AddAt(lx[1], 6);
    par.AddAt(lx[2], 7);
    par.AddAt(tx[0], 8);
    par.AddAt(tx[1], 9);
    par.AddAt(tx[2], 10);
    return kTRUE;
  }
  ::Error("GetShape", "Getting shape parameters for shape %s not implemented", shape->ClassName());
  return kFALSE;
}

Bool_t Engine::GetMaterial(const TString& volumeName, TString& name, Int_t& imat, Double_t& a, Double_t& z, Double_t& density,
                           Double_t& radl, Double_t& inter, TArrayD& par)
{
  auto vol = fGeoManager->GetVolume(volumeName.Data());
  if (!vol) {
    return kFALSE;
  }
  auto med = vol->GetMedium();
  if (!med) {
    return kFALSE;
  }
  auto mat = med->GetMaterial();
  imat = mat->GetUniqueID();
  name = mat->GetName();
  name = name.Strip(TString::kTrailing, '$');
  a = mat->GetA();
  z = mat->GetZ();
  density = mat->GetDensity();
  radl = mat->GetRadLen();
  inter = mat->GetIntLen(); // WARNING: THIS IS NOT COMPUTED NATIVELY BY TGEO
  par.Set(0);               // NO USER PARAMETERS STORED IN TGEO
  return kTRUE;
}

Bool_t Engine::GetMedium(const TString& volumeName, TString& name, Int_t& imed, Int_t& nmat, Int_t& isvol, Int_t& ifield,
                         Double_t& fieldm, Double_t& tmaxfd, Double_t& stemax, Double_t& deemax, Double_t& epsil, Double_t& stmin,
                         TArrayD& par)
{
  auto vol = fGeoManager->GetVolume(volumeName.Data());
  if (!vol)
    return kFALSE;
  TGeoMedium* med = vol->GetMedium();
  if (!med)
    return kFALSE;
  auto mat = med->GetMaterial();
  nmat = mat->GetUniqueID();
  imed = med->GetId();
  name = med->GetName();
  name = name.Strip(TString::kTrailing, '$');
  par.Set(0); // NO USER PARAMETERS IN TGEO
  isvol = (Int_t)med->GetParam(0);
  ifield = (Int_t)med->GetParam(1);
  fieldm = med->GetParam(2);
  tmaxfd = med->GetParam(3);
  stemax = med->GetParam(4);
  deemax = med->GetParam(5);
  epsil = med->GetParam(6);
  stmin = med->GetParam(7);
  return kTRUE;
}

Int_t Engine::VolId(const char* volName) const
{
  auto uid{fGeoManager->GetUID(volName)};
  if (uid < 0) {
    ::Error("VolId: Volume %s not found\n", volName);
    return -1; // TODO This was returning 0 before. Is it possible that 0 actually refers to a valid volume?
  }
  return uid;
}

const char* Engine::VolName(Int_t id) const
{
  auto volume{fGeoManager->GetVolume(id)};
  if (!volume) {
    ::Error("VolName", "volume with id=%d does not exist", id);
    return "NULL";
  }
  return volume->GetName();
}

Int_t Engine::MediumId(const char* mediumName) const
{
  auto medium{fGeoManager->GetMedium(mediumName)};
  if (medium) {
    return medium->GetId();
  }
  return -1;
}

Int_t Engine::NofVolumes() const
{
  return fGeoManager->GetListOfVolumes()->GetEntriesFast() - 1;
}

Int_t Engine::VolId2Mate(Int_t id) const
{
  auto volume{fGeoManager->GetVolume(id)};
  if (!volume) {
    ::Error("VolId2Mate", "volume with id=%d does not exist", id);
    return 0;
  }
  auto med{volume->GetMedium()};
  if (!med) {
    return 0; // TODO Is that right or could 0 in any case be a valid material ID?
  }
  return med->GetId();
}

Int_t Engine::NofVolDaughters(const char* volName) const
{
  auto volume{fGeoManager->GetVolume(volName)};
  if (!volume) {
    ::Error("NofVolDaughters", "Volume %s not found.", volName);
    return 0;
  }
  return volume->GetNdaughters();
}

const char* Engine::VolDaughterName(const char* volName, Int_t i) const
{
  // Get volume
  auto volume{fGeoManager->GetVolume(volName)};
  if (!volume) {
    ::Error("VolDaughterName", "Volume %s not found.", volName);
    return "";
  }

  // Check index
  if (i < 0 || i >= volume->GetNdaughters()) {
    ::Error("VolDaughterName", "Volume %s Index out of limits", volName);
    return "";
  }
  // Return node's volume name
  return volume->GetNode(i)->GetVolume()->GetName();
}

Int_t Engine::VolDaughterCopyNo(const char* volName, Int_t i) const
{
  // Get volume
  auto volume{fGeoManager->GetVolume(volName)};
  if (!volume) {
    ::Error("VolDaughterName", "Volume %s not found.", volName);
    return 0;
  }

  // Check index
  if (i < 0 || i >= volume->GetNdaughters()) {
    ::Error("VolDaughterName", "Volume %s Index out of limits", volName);
    return 0;
  }

  // Return node's copyNo
  return volume->GetNode(i)->GetNumber();
}

bool Engine::isPrimary(int trackId) const
{
  // TODO These checks should't be necessary when dealing with a sane MCStepLogger file
  if (trackId > -1 && trackId < fCurrentLookups->tracktoparent.size()) {
    return fCurrentLookups->tracktoparent[trackId] < 0;
  }
  return false;
}

int Engine::getMediumId(int volId) const
{
  if (volId > -1 && volId < fCurrentLookups->volidtomedium.size()) {
    if (fCurrentLookups->volidtomedium[volId] &&
        fCurrentLookups->volidtomedium[volId]->size() != 0) {
      // extract medium name, from that the TGeoMedium and from that finally the ID
      auto mediumName{*(fCurrentLookups->volidtomedium[volId])};
      auto medium{gGeoManager->GetMedium(mediumName.c_str())};
      if (medium) {
        return medium->GetId();
      }
    }
  }
  return -1;
}

Bool_t Engine::SetProcess(const char* flagName, Int_t flagValue)
{
  return insertProcessOrCut(fProcessesGlobal, physics::namesCuts, flagName, flagValue);
}

Bool_t Engine::SetCut(const char* cutName, Double_t cutValue)
{
  return insertProcessOrCut(fCutsGlobal, physics::namesCuts, cutName, cutValue);
}

Int_t Engine::CurrentVolID(Int_t& copyNo) const
{
  copyNo = fCurrentStep->copyNo;
  return fCurrentStep->volId;
}

Int_t Engine::CurrentVolOffID(Int_t off, Int_t& copyNo) const
{
  if (off < 0 || off > fGeoManager->GetLevel()) {
    return 0;
  }
  if (off == 0) {
    return CurrentVolID(copyNo);
  }
  auto node{fGeoManager->GetMother(off)};
  if (!node) {
    return 0;
  }
  copyNo = node->GetNumber();
  return node->GetVolume()->GetNumber();
}

const char* Engine::CurrentVolName() const
{
  return fCurrentLookups->volidtovolname[fCurrentStep->volId]->c_str();
}

const char* Engine::CurrentVolOffName(Int_t off) const
{
  if (off < 0 || off > fGeoManager->GetLevel()) {
    return 0;
  }
  if (off == 0) {
    return CurrentVolName();
  }
  TGeoNode* node = fGeoManager->GetMother(off);
  if (!node) {
    return 0;
  }
  return node->GetVolume()->GetName();
}

const char* Engine::CurrentVolPath()
{
  return fGeoManager->GetPath();
}

void Engine::Gstpar(Int_t itmed, const char* param, Double_t parval)
{
  if (insertProcessOrCut(fProcesses, physics::namesProcesses, fProcessesGlobal, itmed, param, (int)parval)) {
    return;
  }
  insertProcessOrCut(fCuts, physics::namesCuts, fCutsGlobal, itmed, param, parval);
}

void Engine::loadCurrentCutsAndProcesses(int volId)
{
  fcurrentProcesses = &fProcessesGlobal;
  fcurrentCuts = &fCutsGlobal;
  auto mediumId = getMediumId(volId);
  if (mediumId > -1) {
    if (fProcesses[mediumId]) {
      fcurrentProcesses = fProcesses[mediumId];
    }
    if (fCuts[mediumId]) {
      fcurrentCuts = fCuts[mediumId];
    }
  }
}

bool Engine::keepDueToProcesses(const o2::StepInfo& step) const
{
  return true;
}

bool Engine::keepDueToCuts(const o2::StepInfo& step) const
{
  if (std::abs(fCurrentLookups->tracktopdg[step.trackID]) == 22 && (*fcurrentCuts)[0] >= 0 && (*fcurrentCuts)[0] >= step.E) {
    return false;
  }
  if (std::abs(fCurrentLookups->tracktopdg[step.trackID]) == 11 && (*fcurrentCuts)[1] >= 0 && (*fcurrentCuts)[1] >= step.E) {
    return false;
  }
  return true;
}

bool Engine::keepStep(const o2::StepInfo& step) const
{
  if (!keepDueToProcesses(step)) {
    // TODO could combine both conditions with ||, however, maybe we want to extract the exact reason why this track was stopped
    return false;
  }
  if (!keepDueToCuts(step)) {
    // TODO could combine both conditions with ||, however, maybe we want to extract the exact reason why this track was stopped
    return false;
  }
  return true; // to be implemented
}

void Engine::ProcessEvent(Int_t eventId)
{
  // Push all particles to the stack already (we are not using it but an empty stack might cause the application to immediately stop without particles on the stack)
  // FIXME

  fApplication->BeginEvent();

  std::vector<bool> skipTrack(fCurrentLookups->tracktopdg.size(), false);
  int currentTrackId = -1;
  bool currentIsPrimary = false;
  int previousVolId = -1;
  // If we replay and in particular when certain tracks are killed during replay we have to make sure to obay the indexing og the user's stack
  std::vector<int> userTrackId(fCurrentLookups->tracktopdg.size(), -1);

  // push all primaries already to be consistent with VMC behaviour
  for (auto& step : *fCurrentStepInfo) {
    // by default just assign the previous track ID, the stack might then decide to do something else
    userTrackId[step.trackID] = step.trackID;
    if (step.newtrack && fCurrentLookups->tracktoparent[step.trackID] < 0) {
      fMCStack->PushTrack(0, -1, fCurrentLookups->tracktopdg[step.trackID], -1., -1., -1.,
                          step.E, step.x, step.y, step.z, -1., -1., -1., -1., TMCProcess(step.prodprocess), userTrackId[step.trackID], 1., -1);
    }
  }

  for (const auto& step : *fCurrentStepInfo) {
    // loop over all steps of one event

    if (fIsEventStopped) {
      fIsEventStopped = kFALSE;
      return;
    }

    if (fIsTrackStopped || skipTrack[step.trackID] || (fCurrentLookups->tracktoparent[step.trackID] > -1 && skipTrack[fCurrentLookups->tracktoparent[step.trackID]]) || !keepStep(step)) {
      // even if that track is not flagged to be skipped, the parent track might be, so set it to true also in that case to recursively guarantee that we skip the whole history of a killed track.
      // In that case, this is a child track being a potential parent to other tracks which must not be transported.
      skipTrack[step.trackID] = true;
      continue;
    }

    // Set this before any method from the application is called
    fCurrentStep = const_cast<o2::StepInfo*>(&step);

    if (step.volId != previousVolId) {
      // find the correct set of cuts and processes for this volume
      loadCurrentCutsAndProcesses(step.volId);
      previousVolId = step.volId;
      // TODO Can we find a cheaper way of doing that?
      fGeoManager->FindNode(step.volId);
    }

    if (currentTrackId != step.trackID) {

      if (currentTrackId > -1) {
        // there was a track before
        if (currentIsPrimary) {
          fApplication->FinishPrimary();
        }
        fApplication->PostTrack();
      }
      fIsTrackStopped = kFALSE;
      fCurrentTrackLength = 0.;

      currentTrackId = step.trackID;
      // If primary, it should have the same ID as during original simulation since primaries are always pushed at the beginning all at once. The same we do here (see above)
      if (!isPrimary(currentTrackId)) {
        // by default just assign the previous track ID, the stack might then decide to do something else
        // here we need also be careful to set the correct parent comlying with the user stack indexing
        userTrackId[step.trackID] = step.trackID;
        fMCStack->PushTrack(0, userTrackId[fCurrentLookups->tracktoparent[step.trackID]], fCurrentLookups->tracktopdg[step.trackID], -1., -1., -1.,
                            step.E, step.x, step.y, step.z, -1., -1., -1., -1., TMCProcess(step.prodprocess), userTrackId[step.trackID], 1., -1);
      }

      // Need to set from the mapped user track ID because we don't know the internals of the indexing of that stack but we have to comply with it
      fMCStack->SetCurrentTrack(userTrackId[step.trackID]);

      fApplication->PreTrack();

      currentIsPrimary = false;
      if (isPrimary(currentTrackId)) {
        fApplication->BeginPrimary();
        currentIsPrimary = true;
      }
    }

    // TODO That seems to be the correct way of implementing.
    // However, should we do it before or after calling VMC::Stepping?
    fCurrentTrackLength += step.step;

    // TODO Maybe needed to find out which hit logic is used now by O2?!
    fApplication->Stepping();
  }

  // Finish last track
  if (currentIsPrimary) {
    fApplication->FinishPrimary();
  }
  fApplication->PostTrack();

  // Finish this event
  fApplication->FinishEvent();
}

Bool_t Engine::ProcessRun(Int_t nevent)
{
  fCurrentEvent = 0;

  o2::mcstepanalysis::ROOTIOUtilities treeReader{fFilename};

  if (!treeReader.changeToTTree(fTreename)) {
    Engine::Fatal("Cannot find TTree %s", fTreename.c_str());
  }

  if (treeReader.nEntries() == 0) {
    ::Warning("Engine::ProcessRun", "MCStepLogger tree empty, no events to simulate", fTreename.c_str());
    return kFALSE;
  }

  fMCStack = GetStack();

  if (nevent < 0) {
    // just process all events
    nevent = treeReader.nEntries();
  }

  if (!treeReader.setBranch("Steps", &fCurrentStepInfo) || !treeReader.setBranch("Calls", &fCurrentMagCallInfo) || !treeReader.setBranch("Lookups", &fCurrentLookups)) {
    treeReader.close();
    Engine::Fatal("Engine::ProcessRun", "Cannot find required branches in TTree %s", fTreename.c_str());
  }

  while (treeReader.processTTree() && fCurrentEvent < nevent) {
    ::Info("Engine::ProcessRun", "Process event %d", fCurrentEvent + 1);
    if (fIsRunStopped) {
      fIsRunStopped = kFALSE;
      return kTRUE;
    }
    ProcessEvent(fCurrentEvent);
    fCurrentEvent++;
  }

  treeReader.close();

  return kTRUE;
}
