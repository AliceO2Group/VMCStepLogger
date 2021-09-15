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
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TArrayD.h>
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

#include "TParticle.h"

#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TIterator.h>

#include "MCReplay/MCReplayEngine.h"

ClassImp(mcreplay::MCReplayEngine);

using namespace mcreplay;

MCReplayEngine::MCReplayEngine(const std::string& filename, const std::string& treename)
  : TVirtualMC{"MCReplayEngine", "MCReplayEngine", kTRUE},
    mStepLoggerFilename{filename},
    mStepLoggerTreename{treename},
    mProcessesGlobal(physics::namesProcesses.size(), -1),
    mCutsGlobal(physics::namesCuts.size(), -1.)
{
}

MCReplayEngine::MCReplayEngine() : MCReplayEngine("", "") {}

MCReplayEngine::~MCReplayEngine()
{
  // These are all pointers we own
  // TODO Actually make them unique pointers
  for (auto& m : mProcesses) {
    delete m;
  }
  for (auto& m : mCuts) {
    delete m;
  }
  mStepFile->Close();
}

void MCReplayEngine::Init()
{
  fApplication->AddParticles();
  fApplication->AddIons();
  // must cache it here cause during geometry construction VMC methods might be used
  mGeoManager = gGeoManager;
  // now just run all necessary steps to have the geometry available
  fApplication->ConstructGeometry();
  fApplication->MisalignGeometry();
  fApplication->ConstructOpGeometry();
  fApplication->ConstructSensitiveDetectors();
  fApplication->InitGeometry();
}

void MCReplayEngine::printCurrentCuts() const
{
  std::cout << "CURRENT CUTS\n";
  for (int i = 0; i < mCurrentCuts->size(); i++) {
    std::cout << "  -> " << physics::namesCuts[i] << ": " << (*mCurrentCuts)[i] << "\n";
  }
  std::cout << "---" << std::endl;
}

void MCReplayEngine::adaptToTGeoName(const char* nameIn, char* nameOut) const
{
  auto l = strlen(nameIn);
  if (l >= 79) {
    l = 79;
  }
  for (int i = 0; i < l; i++) {
    nameOut[i] = nameIn[i];
  }
  nameOut[l] = 0;
}

ULong_t MCReplayEngine::makeHash(const o2::StepInfo& step) const
{
  auto asLong = [](decltype(step.x) x) {
    return (ULong_t) * (reinterpret_cast<ULong_t*>(&x));
  };

  ULong_t hash;

  hash = asLong(step.x);
  hash ^= asLong(step.y);
  hash ^= asLong(step.z);
  hash ^= asLong(step.t);
  hash ^= asLong(step.px);
  hash ^= asLong(step.py);
  hash ^= asLong(step.pz);
  hash += (ULong_t)mCurrentLookups->tracktopdg[step.trackID];
  return hash;
}

Double_t* MCReplayEngine::makeDoubleArray(Float_t* arrIn, int np) const
{
  Double_t* arrOut = np > 0 ? new Double_t[np] : nullptr;
  for (int i = 0; i < np; i++) {
    arrOut[i] = arrIn[i];
  }
  return arrOut;
}

void MCReplayEngine::Material(Int_t& kmat, const char* name, Double_t a, Double_t z, Double_t dens, Double_t radl, Double_t absl, Float_t* buf, Int_t nwbuf)
{
  Double_t* bufDouble = nullptr;
  Material(kmat, name, a, z, dens, radl, absl, bufDouble, -1);
}

void MCReplayEngine::Material(Int_t& kmat, const char* name, Double_t a, Double_t z, Double_t dens, Double_t radl, Double_t absl, Double_t* buf, Int_t nwbuf)
{
  kmat = ++mMaterialCounter;
  mGeoManager->Material(name, a, z, dens, kmat, radl, absl);
}

void MCReplayEngine::Mixture(Int_t& kmat, const char* name, Double_t* a, Double_t* z, Double_t dens, Int_t nlmat, Double_t* wmat)
{
  if (nlmat < 0) {
    nlmat = -nlmat;
    Double_t amol = 0;
    Int_t i;
    for (i = 0; i < nlmat; i++) {
      amol += a[i] * wmat[i];
    }
    for (i = 0; i < nlmat; i++) {
      wmat[i] *= a[i] / amol;
    }
  }
  kmat = ++mMaterialCounter;
  mGeoManager->Mixture(name, a, z, dens, nlmat, wmat, kmat);
}

void MCReplayEngine::Mixture(Int_t& kmat, const char* name, Float_t* a, Float_t* z, Double_t dens, Int_t nlmat, Float_t* wmat)
{
  auto aDouble = makeDoubleArray(a, TMath::Abs(nlmat));
  auto zDouble = makeDoubleArray(z, TMath::Abs(nlmat));
  auto wmatDouble = makeDoubleArray(wmat, TMath::Abs(nlmat));
  Mixture(kmat, name, aDouble, zDouble, dens, nlmat, wmatDouble);
  for (int i = 0; i < TMath::Abs(nlmat); i++) {
    a[i] = aDouble[i];
    z[i] = zDouble[i];
    wmat[i] = wmatDouble[i];
  }
  delete[] aDouble;
  delete[] zDouble;
  delete[] wmatDouble;
}

void MCReplayEngine::Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, Double_t stemax, Double_t deemax, Double_t epsil, Double_t stmin, Float_t* ubuf, Int_t nbuf)
{
  auto ubufDouble = makeDoubleArray(ubuf, nbuf);
  Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin, ubufDouble, nbuf);
  delete[] ubufDouble;
}

void MCReplayEngine::Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, Double_t stemax, Double_t deemax, Double_t epsil, Double_t stmin, Double_t* ubuf, Int_t nbuf)
{
  kmed = ++mMediumCounter;
  mGeoManager->Medium(name, kmed, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
}

Int_t MCReplayEngine::Gsvolu(const char* name, const char* shape, Int_t nmed, Double_t* upar, Int_t np)
{
  char nameOut[80];
  adaptToTGeoName(name, nameOut);
  char shapeOut[5];
  adaptToTGeoName(shape, shapeOut);
  auto vol = mGeoManager->Volume(nameOut, shapeOut, nmed, upar, np);
  if (!vol) {
    ::Fatal("Gsvolu", "Could not construct volume %s", name);
    return -1;
  }
  return vol->GetNumber();
}

Int_t MCReplayEngine::Gsvolu(const char* name, const char* shape, Int_t nmed, Float_t* upar, Int_t np)
{
  auto uparDouble = makeDoubleArray(upar, np);
  auto id = Gsvolu(name, shape, nmed, uparDouble, np);
  delete[] uparDouble;
  return id;
}

void MCReplayEngine::Gsdvn(const char* name, const char* mother, Int_t ndiv, Int_t iaxis)
{
  char nameOut[80];
  char nameMotherOut[80];
  adaptToTGeoName(name, nameOut);
  adaptToTGeoName(mother, nameMotherOut);
  mGeoManager->Division(nameOut, nameMotherOut, iaxis, ndiv, 0, 0, 0, "n");
}

void MCReplayEngine::Gsdvn2(const char* name, const char* mother, Int_t ndiv, Int_t iaxis, Double_t c0i, Int_t numed)
{
  char nameOut[80];
  char nameMotherOut[80];
  adaptToTGeoName(name, nameOut);
  adaptToTGeoName(mother, nameMotherOut);
  mGeoManager->Division(nameOut, nameMotherOut, iaxis, ndiv, c0i, 0, numed, "nx");
}

void MCReplayEngine::Gsdvt(const char* name, const char* mother, Double_t step, Int_t iaxis, Int_t numed, Int_t ndvmx)
{
  char nameOut[80];
  char nameMotherOut[80];
  adaptToTGeoName(name, nameOut);
  adaptToTGeoName(mother, nameMotherOut);
  mGeoManager->Division(nameOut, nameMotherOut, iaxis, 0, 0, step, numed, "s");
}

void MCReplayEngine::Gsdvt2(const char* name, const char* mother, Double_t step, Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx)
{
  char nameOut[80];
  char nameMotherOut[80];
  adaptToTGeoName(name, nameOut);
  adaptToTGeoName(mother, nameMotherOut);
  mGeoManager->Division(nameOut, nameMotherOut, iaxis, 0, c0, step, numed, "sx");
}

void MCReplayEngine::Gspos(const char* name, Int_t nr, const char* mother, Double_t x, Double_t y, Double_t z, Int_t irot, const char* konly)
{
  char nameOut[80];
  char nameMotherOut[80];
  adaptToTGeoName(name, nameOut);
  adaptToTGeoName(mother, nameMotherOut);
  TString onlyString{konly};
  onlyString.ToLower();
  bool isOnly = onlyString.Contains("only") ? true : false;
  Double_t* upar = nullptr;
  mGeoManager->Node(nameOut, nr, nameMotherOut, x, y, z, irot, isOnly, upar);
}

void MCReplayEngine::Gsposp(const char* name, Int_t nr, const char* mother, Double_t x, Double_t y, Double_t z, Int_t irot, const char* konly, Double_t* upar, Int_t np)
{
  char nameOut[80];
  char nameMotherOut[80];
  adaptToTGeoName(name, nameOut);
  adaptToTGeoName(mother, nameMotherOut);
  TString onlyString{konly};
  onlyString.ToLower();
  bool isOnly = onlyString.Contains("only") ? true : false;
  mGeoManager->Node(nameOut, nr, nameMotherOut, x, y, z, irot, isOnly, upar, np);
}

void MCReplayEngine::Gsposp(const char* name, Int_t nr, const char* mother, Double_t x, Double_t y, Double_t z, Int_t irot, const char* konly, Float_t* upar, Int_t np)
{
  auto uparDouble = makeDoubleArray(upar, np);
  Gsposp(name, nr, mother, x, y, z, irot, konly, uparDouble, np);
  delete[] uparDouble;
}

void MCReplayEngine::Matrix(Int_t& krot, Double_t thetaX, Double_t phiX, Double_t thetaY, Double_t phiY, Double_t thetaZ, Double_t phiZ)
{
  krot = mGeoManager->GetListOfMatrices()->GetEntriesFast();
  mGeoManager->Matrix(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ);
}

Bool_t MCReplayEngine::GetTransformation(const TString& volumePath, TGeoHMatrix& matrix)
{
  // We have to preserve the modeler state
  mGeoManager->PushPath();
  if (!mGeoManager->cd(volumePath.Data())) {
    mGeoManager->PopPath();
    return kFALSE;
  }
  matrix = *mGeoManager->GetCurrentMatrix();
  mGeoManager->PopPath();
  return kTRUE;
}

Bool_t MCReplayEngine::GetShape(const TString& volumePath, TString& shapeType, TArrayD& par)
{
  Int_t npar;
  mGeoManager->PushPath();
  if (!mGeoManager->cd(volumePath.Data())) {
    mGeoManager->PopPath();
    return kFALSE;
  }
  auto vol = mGeoManager->GetCurrentVolume();
  mGeoManager->PopPath();
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

Bool_t MCReplayEngine::GetMaterial(const TString& volumeName, TString& name, Int_t& imat, Double_t& a, Double_t& z, Double_t& density, Double_t& radl, Double_t& inter, TArrayD& par)
{
  auto vol = mGeoManager->GetVolume(volumeName.Data());
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

Bool_t MCReplayEngine::GetMedium(const TString& volumeName, TString& name, Int_t& imed, Int_t& nmat, Int_t& isvol, Int_t& ifield, Double_t& fieldm, Double_t& tmaxfd, Double_t& stemax, Double_t& deemax, Double_t& epsil, Double_t& stmin, TArrayD& par)
{
  auto vol = mGeoManager->GetVolume(volumeName.Data());
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

Int_t MCReplayEngine::VolId(const char* volName) const
{
  auto uid{mGeoManager->GetUID(volName)};
  if (uid < 0) {
    ::Error("VolId: Volume %s not found\n", volName);
    return -1; // TODO This was returning 0 before. Is it possible that 0 actually refers to a valid volume?
  }
  return uid;
}

const char* MCReplayEngine::VolName(Int_t id) const
{
  auto volume{mGeoManager->GetVolume(id)};
  if (!volume) {
    ::Error("VolName", "volume with id=%d does not exist", id);
    return "NULL";
  }
  return volume->GetName();
}

Int_t MCReplayEngine::MediumId(const char* mediumName) const
{
  auto medium{mGeoManager->GetMedium(mediumName)};
  if (medium) {
    return medium->GetId();
  }
  return -1;
}

Int_t MCReplayEngine::NofVolumes() const
{
  return mGeoManager->GetListOfVolumes()->GetEntriesFast() - 1;
}

Int_t MCReplayEngine::VolId2Mate(Int_t id) const
{
  auto volume{mGeoManager->GetVolume(id)};
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

Int_t MCReplayEngine::NofVolDaughters(const char* volName) const
{
  auto volume{mGeoManager->GetVolume(volName)};
  if (!volume) {
    ::Error("NofVolDaughters", "Volume %s not found.", volName);
    return 0;
  }
  return volume->GetNdaughters();
}

const char* MCReplayEngine::VolDaughterName(const char* volName, Int_t i) const
{
  // Get volume
  auto volume{mGeoManager->GetVolume(volName)};
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

Int_t MCReplayEngine::VolDaughterCopyNo(const char* volName, Int_t i) const
{
  // Get volume
  auto volume{mGeoManager->GetVolume(volName)};
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

bool MCReplayEngine::isPrimary(int trackId) const
{
  // TODO These checks should't be necessary when dealing with a sane MCStepLogger file
  if (trackId > -1 && trackId < mCurrentLookups->tracktoparent.size()) {
    return mCurrentLookups->tracktoparent[trackId] < 0;
  }
  return false;
}

int MCReplayEngine::getMediumId(int volId) const
{
  if (volId > -1 && volId < mCurrentLookups->volidtomedium.size()) {
    if (mCurrentLookups->volidtomedium[volId] &&
        mCurrentLookups->volidtomedium[volId]->size() != 0) {
      // extract medium name, from that the TGeoMedium and from that finally the ID
      auto mediumName{*(mCurrentLookups->volidtomedium[volId])};
      auto medium{gGeoManager->GetMedium(mediumName.c_str())};
      if (medium) {
        return medium->GetId();
      }
    }
  }
  return -1;
}

Bool_t MCReplayEngine::SetProcess(const char* flagName, Int_t flagValue)
{
  return insertProcessOrCut(mProcessesGlobal, physics::namesCuts, flagName, flagValue);
}

Bool_t MCReplayEngine::SetCut(const char* cutName, Double_t cutValue)
{
  return insertProcessOrCut(mCutsGlobal, physics::namesCuts, cutName, cutValue);
}

Int_t MCReplayEngine::CurrentVolID(Int_t& copyNo) const
{
  if (mCurrentStep->outside) {
    return 0;
  }
  copyNo = mCurrentStep->copyNo;
  return mCurrentStep->volId;
}

Int_t MCReplayEngine::CurrentVolOffID(Int_t off, Int_t& copyNo) const
{
  if (off < 0 || off > mGeoManager->GetLevel()) {
    return 0;
  }
  if (off == 0) {
    return CurrentVolID(copyNo);
  }
  auto node{mGeoManager->GetMother(off)};
  if (!node) {
    return 0;
  }
  copyNo = node->GetNumber();
  return node->GetVolume()->GetNumber();
}

const char* MCReplayEngine::CurrentVolName() const
{
  if (mCurrentStep->outside) {
    return mGeoManager->GetTopVolume()->GetName();
  }
  return mCurrentLookups->volidtovolname[mCurrentStep->volId]->c_str();
}

const char* MCReplayEngine::CurrentVolOffName(Int_t off) const
{
  if (off < 0 || off > mGeoManager->GetLevel()) {
    return 0;
  }
  if (off == 0) {
    return CurrentVolName();
  }
  TGeoNode* node = mGeoManager->GetMother(off);
  if (!node) {
    return 0;
  }
  return node->GetVolume()->GetName();
}

const char* MCReplayEngine::CurrentVolPath()
{
  return mGeoManager->GetPath();
}

void MCReplayEngine::Gmtod(Double_t* xm, Double_t* xd, Int_t iflag)
{
  if (iflag == 1) {
    gGeoManager->MasterToLocal(xm, xd);
  } else {
    gGeoManager->MasterToLocalVect(xm, xd);
  }
}

void MCReplayEngine::Gmtod(Float_t* xm, Float_t* xd, Int_t iflag)
{
  Double_t xmDouble[3], xdDouble[3];
  for (int i = 0; i < 3; i++) {
    xmDouble[i] = xm[i];
  }
  Gmtod(xmDouble, xdDouble, iflag);
  for (int i = 0; i < 3; i++) {
    xd[i] = xdDouble[i];
  }
}

void MCReplayEngine::Gdtom(Double_t* xd, Double_t* xm, Int_t iflag)
{
  if (iflag == 1) {
    gGeoManager->LocalToMaster(xd, xm);
  } else {
    gGeoManager->LocalToMasterVect(xd, xm);
  }
}

void MCReplayEngine::Gdtom(Float_t* xd, Float_t* xm, Int_t iflag)
{
  Double_t xmDouble[3], xdDouble[3];
  for (int i = 0; i < 3; i++) {
    xdDouble[i] = xd[i];
  }
  Gdtom(xdDouble, xmDouble, iflag);
  for (int i = 0; i < 3; i++) {
    xm[i] = xmDouble[i];
  }
}

void MCReplayEngine::Gstpar(Int_t itmed, const char* param, Double_t parval)
{
  if (insertProcessOrCut(mProcesses, physics::namesProcesses, mProcessesGlobal, itmed, param, (int)parval)) {
    return;
  }
  insertProcessOrCut(mCuts, physics::namesCuts, mCutsGlobal, itmed, param, parval);
}

void MCReplayEngine::loadCurrentCutsAndProcesses(int volId)
{
  mCurrentProcesses = &mProcessesGlobal;
  mCurrentCuts = &mCutsGlobal;
  auto mediumId = getMediumId(volId);
  if (mediumId > -1) {
    if (mProcesses.size() > mediumId && mProcesses[mediumId]) {
      mCurrentProcesses = mProcesses[mediumId];
    }
    if (mCuts.size() > mediumId && mCuts[mediumId]) {
      mCurrentCuts = mCuts[mediumId];
    }
  }
}

bool MCReplayEngine::keepDueToProcesses(const o2::StepInfo& step) const
{
  return true;
}

bool MCReplayEngine::keepDueToCuts(const o2::StepInfo& step) const
{
  if ((*mCurrentCuts)[11] > 0. && step.E < (*mCurrentCuts)[11]) {
    // check global energy cut
    return false;
  }
  return true;
}

bool MCReplayEngine::keepStep(const o2::StepInfo& step) const
{
  return keepDueToCuts(step) && keepDueToProcesses(step);
}

bool MCReplayEngine::initRun()
{
  if (mIsInitialised) {
    return true;
  }

  mStepFile = TFile::Open(mStepLoggerFilename.c_str(), "READ");
  if (mStepFile == nullptr) {
    ::Error("::MCReplayEngine::initRun", "Cannot open file %s", mStepLoggerFilename.c_str());
    return false;
  }

  auto tree = (TTree*)mStepFile->Get(mStepLoggerTreename.c_str());
  mStepBranch = tree->GetBranch("Steps");
  mLookupBranch = tree->GetBranch("Lookups");
  if (!mStepBranch || !mLookupBranch) {
    ::Error("::MCReplayEngine::initRun", "Cannot get branches \"Steps\" and \"Lookups\"");
    return false;
  }

  mStack = GetStack();

  mIsInitialised = true;
  return true;
}

void MCReplayEngine::transportUserHitSecondary()
{
  // For now we just start and finish this secondary. This is only done to pretend the transport for the user stack
  fApplication->PreTrack();
  fApplication->PostTrack();
}

void MCReplayEngine::ProcessEvent(Int_t eventId)
{
  if (!initRun()) {
    return;
  }

  // prepare
  std::vector<o2::StepInfo>* steps = nullptr;
  mCurrentLookups = nullptr;
  mStepBranch->SetAddress(&steps);
  mLookupBranch->SetAddress(&mCurrentLookups);
  mStepBranch->GetEvent(mCurrentEvent);
  mLookupBranch->GetEvent(mCurrentEvent);

  // whether or not to skip certain tracks
  std::vector<bool> skipTrack(mCurrentLookups->tracktopdg.size(), false);
  // we need to make sure we follow the indexing of the user stack. During the original simulation, there might have been more tracks pushed than transported. In the replay case, we only have the tracks that have been originally transported. Hence, the indexing this time might be different.
  std::vector<int> userTrackId(mCurrentLookups->tracktopdg.size(), -1);
  // some caching to be able to run pre- and post-hooks at the right time
  int currentTrackId = -1;
  // some caching to be able to update the TGeoManager state
  int previousVolId = -1;
  int previousCopyNo = -1;

  fApplication->BeginEvent();
  fApplication->GeneratePrimaries();

  // Remember how many primaries are expected to be transported
  auto expectedNPrimaries = mStack->GetNprimary();

  /* SOME REMARKS

  1. Steps are expected to be grouped together track-by-track (that comes from the MCStepLogger and will only work properly for the serial simulation run)
  2. Therefore, if a step with a different track ID is reached, it is assumed that the previous track is finished
  3. Since in the original reference un, some secondaries might have been pushed to the user stack but not transported, we comply with the track ID assignement of the user stack in a replay run
     See also comments in the corresponding section below
  4. Whenever a new node (different volume ID or copy number) is entered, we let the TGeoManager find that node (see further comments below)
  5. The current approach definitely works for GEANT3 in which case the next primary is popped and transported only after all secondaries of the previous primary have been transported.
     In that case TVirtualMCApplication::FinishPrimary() is called after the last secondary and before the next primary
     ==> Might need some refinement for GEANT4
     ==> TODO Most definitely, the geometry construction between GEANT3 and GEANT4 differs, so we find different volume IDs, in particular for sensitive volumes ==> NEeds fixing/clarification

  In case no additional cuts are applied during the replay:
  1. At least with GEANT3 reference runs, steps are exactly reproduced comparing the reference run and a replay
  2. hits are exactly reproduced for each replay run
  */

  unsigned int nStepsKept{0};

  for (const auto& step : *steps) {

    // TODO This should not happen. (see comment in header file for these flags)
    if (mIsEventStopped) {
      mIsEventStopped = false;
      ::Warning("MCReplayEngine::ProcessEvent", "Event %d was stopped. That should usually not happen", mCurrentEvent);
      delete steps;
      delete mCurrentLookups;
      return;
    }

    // skip if flagged and increment
    if (skipTrack[step.trackID]) {
      continue;
    }
    if (mCurrentLookups->tracktoparent[step.trackID] > -1 && skipTrack[mCurrentLookups->tracktoparent[step.trackID]]) {
      // skip recursively in case parent was skipped, only affects secondaries of course
      skipTrack[step.trackID] = true;
      continue;
    }

    if (step.volId != previousVolId || step.copyNo != previousCopyNo) {
      // find the correct set of cuts and processes for this volume
      loadCurrentCutsAndProcesses(step.volId);
      previousVolId = step.volId;
      previousCopyNo = step.copyNo;
      // We need the navigator in the current TGeoNode in order to perform tasks like CurrentVolOffID.
      // TODO Can we find a cheaper way of doing that?
      // It does not seem to be possible to just get a TGeoNode even if the volume ID and copy number are known...
      // Indeed, that could be a cheaper way but for the moment we stick to the TGeoManager
      // NOTE/TODO Actually, this could MAYBE go wrong hence finding the previous volume, for instance in case the step was made close to/at the boundary
      mGeoManager->FindNode(step.x, step.y, step.z);
    }

    // NOW PERFORM EVERYTHIN NECESSARY TO START / FINISH A TRACK

    if (currentTrackId == step.trackID && mIsTrackStopped) {
      // as the warning says, it should not happen but could (e.g. due to double vs. float precision)
      // but mainly due to the fact that the RNG (TRandom) in the reference run and in this replay do not have the same seeding
      ::Warning("MCReplayEngine::ProcessEvent", "Track %d was stopped (PDG %d). That should usually not happen", userTrackId[step.trackID], mCurrentLookups->tracktopdg[step.trackID]);
      mIsTrackStopped = false;
    }

    if (currentTrackId != step.trackID) {
      // Found a new track

      auto prim = isPrimary(step.trackID);

      if (currentTrackId > -1) {
        // only invoke if there has been at least 1 track before
        fApplication->PostTrack();
        if (prim) {
          // finish previous primary only in case all its secondary tracks have been transported and a new primary is about to be transported
          fApplication->FinishPrimary();
        }
      }

      if (!prim) {
        // decide to skip this secondary immediately, primaries must be popped first since otherwise we might run into inconsitencies with the user stack
        // therefore, it looks as if this secondary never appeared
        if (!keepStep(step)) {
          skipTrack[step.trackID] = true;
          continue;
        }

        // only push track if it is a secondary
        // by default just assign the MCStepLogger track ID, but the user stack might then decide to do something else
        userTrackId[step.trackID] = step.trackID;
        mStack->PushTrack(1, userTrackId[mCurrentLookups->tracktoparent[step.trackID]], mCurrentLookups->tracktopdg[step.trackID], step.px, step.py, step.pz, step.E, step.x, step.y, step.z, step.t, 1., 1., 1., TMCProcess(step.prodprocess), userTrackId[step.trackID], 1., 0);
        mStack->PopNextTrack(userTrackId[step.trackID]);
      } else {
        // We ignore all tracks which were pushed somehow during hit creation, for instance based on an RNG. If we wouldn't do that, it would lead to an incoherent stack compared to the logged steps
        // That indeed happens in case of the O2 TRD
        int popTrackId;
        mStack->PopNextTrack(popTrackId);
        // NOTE Only in case of primaires we assume their indexing to be 0...nPrimaries-1
        while (popTrackId >= expectedNPrimaries) {
          // In such a case there are secondaries pushed during hit creation...
          transportUserHitSecondary();
          mStack->PopNextTrack(popTrackId);
        }
        userTrackId[step.trackID] = popTrackId;
      }

      currentTrackId = step.trackID;
      mStack->SetCurrentTrack(userTrackId[step.trackID]);
      mCurrentTrackLength = 0.;

      if (prim) {
        fApplication->BeginPrimary();
      }
      fApplication->PreTrack();

      if (!keepStep(step)) {
        // This is for primaries since secondaries would have been skipped already above
        skipTrack[step.trackID] = true;
        continue;
      }

      // We fix the RNG so in case the hits are somehow based on that, the hits among different Replay runs are exactly reproduced
      gRandom->SetSeed(makeHash(step));
    }

    mCurrentTrackLength += step.step;
    mCurrentStep = const_cast<o2::StepInfo*>(&step);

    nStepsKept++;
    fApplication->Stepping();
  }

  // Finish last track and the entire event
  fApplication->PostTrack();
  fApplication->FinishPrimary();

  fApplication->FinishEvent();

  auto nStepsSkipped = steps->size() - nStepsKept;
  std::cout << "Original number, skipped, kept, skipped fraction and kept fraction of steps: " << steps->size() << " " << nStepsSkipped << " " << nStepsKept << " " << static_cast<float>(nStepsSkipped) / steps->size() << " " << static_cast<float>(nStepsSkipped) / steps->size() << "\n";

  delete steps;
  delete mCurrentLookups;

  mCurrentEvent++;
}

Bool_t MCReplayEngine::ProcessRun(Int_t nevent)
{
  for (int i = 0; i < nevent; i++) {
    ProcessEvent(mCurrentEvent);
  }
  return kTRUE;
}
