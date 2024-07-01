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

#ifndef MC_REPLAY_ENGINE_H
#define MC_REPLAY_ENGINE_H

#include <vector>
#include <string>

#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMCProcess.h>

#include <TVirtualMC.h>

#include "MCStepLogger/StepInfo.h"
#include "MCReplay/MCReplayPhysics.h"

class TVirtualMCStack;
class TBranch;
class TFile;
class TGeoManager;
class TGeoVolume;
class TGeoNode;

namespace mcreplay
{

typedef bool (*user_keep_step_type)(const o2::StepInfo&, const o2::StepLookups*);

class MCReplayEngine : public TVirtualMC
{

 public:
  MCReplayEngine(const std::string& filename, const std::string& treename);
  MCReplayEngine();
  MCReplayEngine(MCReplayEngine const&) = delete;
  MCReplayEngine& operator=(MCReplayEngine const&) = delete;

  /// For now just default destructor
  virtual ~MCReplayEngine();

  //
  // All the derived stuff
  //

  //
  // ------------------------------------------------
  // methods for building/management of geometry
  // ------------------------------------------------
  //

  /// Info about supporting geometry defined via Root
  Bool_t IsRootGeometrySupported() const override
  {
    return kTRUE;
  }

  //
  // functions from GCONS
  // ------------------------------------------------
  //

  /// Define a material
  /// - kmat   number assigned to the material
  /// - name   material name
  /// - a      atomic mass in au
  /// - z      atomic number
  /// - dens   density in g/cm3
  /// - absl   absorption length in cm;
  ///               if >=0 it is ignored and the program
  ///               calculates it, if <0. -absl is taken
  /// - radl   radiation length in cm
  ///               if >=0 it is ignored and the program
  ///               calculates it, if <0. -radl is taken
  /// - buf    pointer to an array of user words
  /// - nwbuf  number of user words
  void Material(Int_t& kmat, const char* name, Double_t a, Double_t z, Double_t dens, Double_t radl, Double_t absl, Float_t* buf, Int_t nwbuf) override;

  /// The same as previous but in double precision
  void Material(Int_t& kmat, const char* name, Double_t a, Double_t z, Double_t dens, Double_t radl, Double_t absl, Double_t* buf, Int_t nwbuf) override;

  /// Define a mixture or a compound
  /// with a number kmat composed by the basic nlmat materials defined
  /// by arrays a, z and wmat
  ///
  /// If nlmat > 0 then wmat contains the proportion by
  /// weights of each basic material in the mixture.
  ///
  /// If nlmat < 0 then wmat contains the number of atoms
  /// of a given kind into the molecule of the compound.
  /// In this case, wmat in output is changed to relative
  /// weights.
  void Mixture(Int_t& kmat, const char* name, Float_t* a, Float_t* z, Double_t dens, Int_t nlmat, Float_t* wmat) override;

  /// The same as previous but in double precision
  void Mixture(Int_t& kmat, const char* name, Double_t* a, Double_t* z, Double_t dens, Int_t nlmat, Double_t* wmat) override;

  /// Define a medium.
  /// - kmed      tracking medium number assigned
  /// - name      tracking medium name
  /// - nmat      material number
  /// - isvol     sensitive volume flag
  /// - ifield    magnetic field:
  ///                  - ifield = 0 if no magnetic field;
  ///                  - ifield = -1 if user decision in guswim;
  ///                  - ifield = 1 if tracking performed with g3rkuta;
  ///                  - ifield = 2 if tracking performed with g3helix;
  ///                  - ifield = 3 if tracking performed with g3helx3.
  /// - fieldm    max. field value (kilogauss)
  /// - tmaxfd    max. angle due to field (deg/step)
  /// - stemax    max. step allowed
  /// - deemax    max. fraction of energy lost in a step
  /// - epsil     tracking precision (cm)
  /// - stmin     min. step due to continuous processes (cm)
  /// - ubuf      pointer to an array of user words
  /// - nbuf      number of user words
  void Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, Double_t stemax, Double_t deemax, Double_t epsil, Double_t stmin, Float_t* ubuf, Int_t nbuf) override;

  /// The same as previous but in double precision
  void Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, Double_t stemax, Double_t deemax, Double_t epsil, Double_t stmin, Double_t* ubuf, Int_t nbuf) override;

  /// Define a rotation matrix
  /// - krot     rotation matrix number assigned
  /// - thetaX   polar angle for axis X
  /// - phiX     azimuthal angle for axis X
  /// - thetaY   polar angle for axis Y
  /// - phiY     azimuthal angle for axis Y
  /// - thetaZ   polar angle for axis Z
  /// - phiZ     azimuthal angle for axis Z
  void Matrix(Int_t& krot, Double_t thetaX, Double_t phiX, Double_t thetaY, Double_t phiY, Double_t thetaZ, Double_t phiZ) override;

  /// Change the value of cut or mechanism param
  /// to a new value parval for tracking medium itmed.
  /// In Geant3, the  data  structure JTMED contains the standard tracking
  /// parameters (CUTS and flags to control the physics processes)  which
  /// are used  by default for all tracking media.
  /// It is possible to redefine individually with this function any of these
  /// parameters for a given tracking medium.
  /// - itmed   tracking medium number
  /// - param   is a character string (variable name)
  /// - parval  must be given as a floating point.
  void Gstpar(Int_t itmed, const char* param, Double_t parval) override;

  //
  // functions from GGEOM
  // ------------------------------------------------
  //

  /// Create a new volume
  /// - name   Volume name
  /// - shape  Volume type
  /// - nmed   Tracking medium number
  /// - np     Number of shape parameters
  /// - upar   Vector containing shape parameters
  Int_t Gsvolu(const char* name, const char* shape, Int_t nmed, Float_t* upar, Int_t np) override;

  /// The same as previous but in double precision
  Int_t Gsvolu(const char* name, const char* shape, Int_t nmed, Double_t* upar, Int_t np) override;

  /// Create a new volume by dividing an existing one.
  /// It divides a previously defined volume
  /// - name   Volume name
  /// - mother Mother volume name
  /// - ndiv   Number of divisions
  /// - iaxis  Axis value:
  ///               X,Y,Z of CAXIS will be translated to 1,2,3 for IAXIS.
  void Gsdvn(const char* name, const char* mother, Int_t ndiv, Int_t iaxis) override;

  /// Create a new volume by dividing an existing one.
  /// Divide mother into ndiv divisions called name
  /// along axis iaxis starting at coordinate value c0i.
  /// The new volume created will be medium number numed.
  void Gsdvn2(const char* name, const char* mother, Int_t ndiv, Int_t iaxis, Double_t c0i, Int_t numed) override;

  /// Create a new volume by dividing an existing one
  /// Divide mother into divisions called name along
  /// axis iaxis in steps of step. If not exactly divisible
  /// will make as many as possible and will center them
  /// with respect to the mother. Divisions will have medium
  /// number numed. If numed is 0, numed of mother is taken.
  /// ndvmx is the expected maximum number of divisions
  /// (If 0, no protection tests are performed in Geant3)
  void Gsdvt(const char* name, const char* mother, Double_t step, Int_t iaxis, Int_t numed, Int_t ndvmx) override;

  /// Create a new volume by dividing an existing one
  /// Divides mother into divisions called name along
  /// axis iaxis starting at coordinate value c0 with step
  /// size step.
  /// The new volume created will have medium number numed.
  /// If numed is 0, numed of mother is taken.
  /// ndvmx is the expected maximum number of divisions
  /// (If 0, no protection tests are performed in Geant3)
  void Gsdvt2(const char* name, const char* mother, Double_t step, Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx) override;

  /// Flag volume name whose contents will have to be ordered
  /// along axis iax, by setting the search flag to -iax
  /// (Geant3 only)
  void Gsord(const char* name, Int_t iax) override
  {
    Warning("Gsord", "Not supported");
  }

  /// Position a volume into an existing one.
  /// It positions a previously defined volume in the mother.
  /// - name   Volume name
  /// - nr     Copy number of the volume
  /// - mother Mother volume name
  /// - x      X coord. of the volume in mother ref. sys.
  /// - y      Y coord. of the volume in mother ref. sys.
  /// - z      Z coord. of the volume in mother ref. sys.
  /// - irot   Rotation matrix number w.r.t. mother ref. sys.
  /// - konly  ONLY/MANY flag
  void Gspos(const char* name, Int_t nr, const char* mother, Double_t x, Double_t y, Double_t z, Int_t irot, const char* konly = "ONLY") override;

  /// Place a copy of generic volume name with user number
  ///  nr inside mother, with its parameters upar(1..np)
  void Gsposp(const char* name, Int_t nr, const char* mother, Double_t x, Double_t y, Double_t z, Int_t irot, const char* konly, Float_t* upar, Int_t np) override;

  /// The same as previous but in double precision
  void Gsposp(const char* name, Int_t nr, const char* mother, Double_t x, Double_t y, Double_t z, Int_t irot, const char* konly, Double_t* upar, Int_t np) override;

  /// Helper function for resolving MANY.
  /// Specify the ONLY volume that overlaps with the
  /// specified MANY and has to be substracted.
  /// (Geant4 only)
  void Gsbool(const char* onlyVolName, const char* manyVolName) override
  {
    Warning("Gsbool", "Not supported");
  }

  /// Define the tables for UV photon tracking in medium itmed.
  /// Please note that it is the user's responsibility to
  /// provide all the coefficients:
  /// - itmed       Tracking medium number
  /// - npckov      Number of bins of each table
  /// - ppckov      Value of photon momentum (in GeV)
  /// - absco       Absorption coefficients
  ///                     - dielectric: absorption length in cm
  ///                     - metals    : absorption fraction (0<=x<=1)
  /// - effic       Detection efficiency for UV photons
  /// - rindex      Refraction index (if=0 metal)
  void SetCerenkov(Int_t itmed, Int_t npckov, Float_t* ppckov, Float_t* absco, Float_t* effic, Float_t* rindex, Bool_t aspline = false, Bool_t rspline = false) override
  {
    Warning("SetCerenkov", "Not supported");
  }

  /// The same as previous but in double precision
  void SetCerenkov(Int_t itmed, Int_t npckov, Double_t* ppckov, Double_t* absco, Double_t* effic, Double_t* rindex, Bool_t aspline = false, Bool_t rspline = false) override
  {
    Warning("SetCerenkov", "Not supported");
  }

  //
  // functions for definition of surfaces
  // and material properties for optical physics
  // ------------------------------------------------
  //

  /// Define the optical surface
  /// - name           surface name
  /// - model          selection of model (see #EMCOpSurfaceModel values)
  /// - surfaceType    surface type (see #EMCOpSurfaceType values)
  /// - surfaceFinish  surface quality (see #EMCOpSurfaceType values)
  /// - sigmaAlpha     an unified model surface parameter
  /// (Geant4 only)
  void DefineOpSurface(const char* name, EMCOpSurfaceModel model, EMCOpSurfaceType surfaceType, EMCOpSurfaceFinish surfaceFinish, Double_t sigmaAlpha) override
  {
    Warning("DefineOpSurface", "Not supported");
  }

  /// Define the optical surface border
  /// - name        border surface name
  /// - vol1Name    first volume name
  /// - vol1CopyNo  first volume copy number
  /// - vol2Name    second volume name
  /// - vol2CopyNo  second volume copy number
  /// - opSurfaceName  name of optical surface which this border belongs to
  /// (Geant4 only)
  void SetBorderSurface(const char* name, const char* vol1Name, int vol1CopyNo, const char* vol2Name, int vol2CopyNo, const char* opSurfaceName) override
  {
    Warning("SetBorderSurface", "Not supported");
  }

  /// Define the optical skin surface
  /// - name        skin surface name
  /// - volName     volume name
  /// - opSurfaceName  name of optical surface which this border belongs to
  /// (Geant4 only)
  void SetSkinSurface(const char* name, const char* volName, const char* opSurfaceName) override
  {
    Warning("SetSkinSurface", "Not supported");
  }

  /// Define material property via a table of values
  /// - itmed         tracking medium id
  /// - propertyName  property name
  /// - np            number of bins of the table
  /// - pp            value of photon momentum (in GeV)
  /// - values        property values
  /// (Geant4 only)
  void SetMaterialProperty(Int_t itmed, const char* propertyName, Int_t np, Double_t* pp, Double_t* values, Bool_t createNewKey = false, Bool_t spline = false) override
  {
    Warning("SetMaterialProperty", "Not supported");
  }

  /// Define material property via a value
  /// - itmed         tracking medium id
  /// - propertyName  property name
  /// - value         property value
  /// (Geant4 only)
  void SetMaterialProperty(Int_t itmed, const char* propertyName, Double_t value) override
  {
    Warning("SetMaterialProperty", "Not supported");
  }

  /// Define optical surface property via a table of values
  /// - surfaceName   optical surface name
  /// - propertyName  property name
  /// - np            number of bins of the table
  /// - pp            value of photon momentum (in GeV)
  /// - values        property values
  /// (Geant4 only)
  void SetMaterialProperty(const char* surfaceName, const char* propertyName, Int_t np, Double_t* pp, Double_t* values, Bool_t createNewKey = false, Bool_t spline = false) override
  {
    Warning("SetMaterialProperty", "Not supported");
  }

  //
  // functions for access to geometry
  // ------------------------------------------------
  //

  /// Return the transformation matrix between the volume specified by
  /// the path volumePath and the top or master volume.
  Bool_t GetTransformation(const TString& volumePath, TGeoHMatrix& matrix) override;

  /// Return the name of the shape (shapeType)  and its parameters par
  /// for the volume specified by the path volumePath .
  Bool_t GetShape(const TString& volumePath, TString& shapeType, TArrayD& par) override;

  /// Return the material parameters for the material specified by
  /// the material Id
  Bool_t GetMaterial(Int_t imat, TString& name, Double_t& a, Double_t& z, Double_t& density, Double_t& radl, Double_t& inter, TArrayD& par) override
  {
    Warning("GetMaterial", "Not yet implemented");
    return kFALSE;
  }

  /// Return the material parameters for the volume specified by
  /// the volumeName.
  Bool_t GetMaterial(const TString& volumeName, TString& name, Int_t& imat, Double_t& a, Double_t& z, Double_t& density, Double_t& radl, Double_t& inter, TArrayD& par) override;

  /// Return the medium parameters for the volume specified by the
  /// volumeName.
  Bool_t GetMedium(const TString& volumeName, TString& name, Int_t& imed, Int_t& nmat, Int_t& isvol, Int_t& ifield, Double_t& fieldm, Double_t& tmaxfd, Double_t& stemax, Double_t& deemax, Double_t& epsil, Double_t& stmin, TArrayD& par) override;

  /// Write out the geometry of the detector in EUCLID file format
  /// - filnam  file name - will be with the extension .euc                 *
  /// - topvol  volume name of the starting node
  /// - number  copy number of topvol (relevant for gsposp)
  /// - nlevel  number of  levels in the tree structure
  ///                to be written out, starting from topvol
  /// (Geant3 only)
  /// Deprecated
  void WriteEuclid(const char* filnam, const char* topvol, Int_t number, Int_t nlevel) override
  {
    Warning("WriteEuclid", "Not supported");
  }

  /// Set geometry from Root (built via TGeo)
  void SetRootGeometry() override {}

  /// Activate the parameters defined in tracking media
  /// (DEEMAX, STMIN, STEMAX), which are, be default, ignored.
  /// In Geant4 case, only STEMAX is taken into account.
  /// In FLUKA, all tracking media parameters are ignored.
  void SetUserParameters(Bool_t isUserParameters) override
  {
    Warning("SetUserParameters", "Not supported");
  }

  //
  // get methods
  // ------------------------------------------------
  //

  /// Return the unique numeric identifier for volume name volName
  Int_t VolId(const char* volName) const override;

  /// Return the volume name for a given volume identifier id
  const char* VolName(Int_t id) const override;

  /// Return the unique numeric identifier for medium name mediumName
  Int_t MediumId(const char* mediumName) const override;

  /// Return total number of volumes in the geometry
  Int_t NofVolumes() const override;

  /// Return material number for a given volume id
  Int_t VolId2Mate(Int_t id) const override;

  /// Return number of daughters of the volume specified by volName
  Int_t NofVolDaughters(const char* volName) const override;

  /// Return the name of i-th daughter of the volume specified by volName
  const char* VolDaughterName(const char* volName, Int_t i) const override;

  /// Return the copyNo of i-th daughter of the volume specified by volName
  Int_t VolDaughterCopyNo(const char* volName, Int_t i) const override;

  //
  // ------------------------------------------------
  // methods for sensitive detectors
  // ------------------------------------------------
  //

  // Set a sensitive detector to a volume
  // - volName - the volume name
  // - sd - the user sensitive detector
  void SetSensitiveDetector(const TString& volName, TVirtualMCSensitiveDetector* sd) override
  {
    Warning("SetSensitiveDetector", "Not supported");
  }

  // Get a sensitive detector of a volume
  // - volName - the volume name
  TVirtualMCSensitiveDetector* GetSensitiveDetector(const TString& volName) const override
  {
    Warning("GetSensitiveDetector", "Not supported");
    return nullptr;
  }

  // The scoring option:
  // if true, scoring is performed only via user defined sensitive detectors and
  // MCApplication::Stepping is not called
  void SetExclusiveSDScoring(Bool_t exclusiveSDScoring) override
  {
    Warning("SetExclusiveSDScoring", "Not supported");
  }

  //
  // ------------------------------------------------
  // methods for physics management
  // ------------------------------------------------
  //

  //
  // set methods
  // ------------------------------------------------
  //

  /// Set transport cuts for particles
  Bool_t SetCut(const char* cutName, Double_t cutValue) override;

  /// Set process control
  Bool_t SetProcess(const char* flagName, Int_t flagValue) override;

  /// Set a user defined particle
  /// Function is ignored if particle with specified pdg
  /// already exists and error report is printed.
  /// - pdg           PDG encoding
  /// - name          particle name
  /// - mcType        VMC Particle type
  /// - mass          mass [GeV]
  /// - charge        charge [eplus]
  /// - lifetime      time of life [s]
  /// - pType         particle type as in Geant4
  /// - width         width [GeV]
  /// - iSpin         spin
  /// - iParity       parity
  /// - iConjugation  conjugation
  /// - iIsospin      isospin
  /// - iIsospinZ     isospin - #rd component
  /// - gParity       gParity
  /// - lepton        lepton number
  /// - baryon        baryon number
  /// - stable        stability
  /// - shortlived    is shorlived?
  /// - subType       particle subType as in Geant4
  /// - antiEncoding  anti encoding
  /// - magMoment     magnetic moment
  /// - excitation    excitation energy [GeV]
  Bool_t DefineParticle(Int_t pdg, const char* name, TMCParticleType mcType, Double_t mass, Double_t charge, Double_t lifetime) override
  {
    Warning("DefineParticle", "Not yet implemented");
    return kFALSE;
  }

  /// Set a user defined particle
  /// Function is ignored if particle with specified pdg
  /// already exists and error report is printed.
  /// - pdg           PDG encoding
  /// - name          particle name
  /// - mcType        VMC Particle type
  /// - mass          mass [GeV]
  /// - charge        charge [eplus]
  /// - lifetime      time of life [s]
  /// - pType         particle type as in Geant4
  /// - width         width [GeV]
  /// - iSpin         spin
  /// - iParity       parity
  /// - iConjugation  conjugation
  /// - iIsospin      isospin
  /// - iIsospinZ     isospin - #rd component
  /// - gParity       gParity
  /// - lepton        lepton number
  /// - baryon        baryon number
  /// - stable        stability
  /// - shortlived    is shorlived?
  /// - subType       particle subType as in Geant4
  /// - antiEncoding  anti encoding
  /// - magMoment     magnetic moment
  /// - excitation    excitation energy [GeV]
  Bool_t DefineParticle(Int_t pdg, const char* name, TMCParticleType mcType, Double_t mass, Double_t charge, Double_t lifetime,
                        const TString& pType, Double_t width, Int_t iSpin, Int_t iParity, Int_t iConjugation,
                        Int_t iIsospin, Int_t iIsospinZ, Int_t gParity, Int_t lepton, Int_t baryon, Bool_t stable,
                        Bool_t shortlived = kFALSE, const TString& subType = "", Int_t antiEncoding = 0, Double_t magMoment = 0.0,
                        Double_t excitation = 0.0) override
  {
    Warning("DefineParticle", "Not yet implemented");
    return kFALSE;
  }

  /// Set a user defined ion.
  /// - name          ion name
  /// - Z             atomic number
  /// - A             atomic mass
  /// - Q             charge [eplus}
  /// - excitation    excitation energy [GeV]
  /// - mass          mass  [GeV] (if not specified by user, approximative
  ///                 mass is calculated)
  Bool_t DefineIon(const char* name, Int_t Z, Int_t A, Int_t Q, Double_t excEnergy, Double_t mass = 0.) override
  {
    Warning("DefineIon", "Not yet implemented");
    return kFALSE;
  }

  /// Set a user phase space decay for a particle
  /// -  pdg           particle PDG encoding
  /// -  bratios       the array with branching ratios (in %)
  /// -  mode[6][3]    the array with daughters particles PDG codes  for each
  ///                 decay channel
  Bool_t SetDecayMode(Int_t pdg, Float_t bratio[6], Int_t mode[6][3]) override
  {
    Warning("SetDecayMode", "Not yet implemented");
    return kFALSE;
  }

  /// Calculate X-sections
  /// (Geant3 only)
  /// Deprecated
  Double_t Xsec(char*, Double_t, Int_t, Int_t) override
  {
    Warning("Xsec", "Not yet implemented");
    return -1.;
  }

  //
  // particle table usage
  // ------------------------------------------------
  //

  /// Return MC specific code from a PDG and pseudo ENDF code (pdg)
  Int_t IdFromPDG(Int_t pdg) const override
  {
    Warning("IdFromPDG", "Not yet implemented");
    return -1;
  }

  /// Return PDG code and pseudo ENDF code from MC specific code (id)
  Int_t PDGFromId(Int_t id) const override
  {
    Warning("PDGFromId", "Not yet implemented");
    return -1;
  }

  //
  // get methods
  // ------------------------------------------------
  //

  /// Return name of the particle specified by pdg.
  TString ParticleName(Int_t pdg) const override
  {
    Warning("ParticleName", "Not yet implemented");
    return TString();
  }

  /// Return mass of the particle specified by pdg.
  Double_t ParticleMass(Int_t pdg) const override
  {
    Warning("ParticleMass", "Not yet implemented");
    return -1.;
  }

  /// Return charge (in e units) of the particle specified by pdg.
  Double_t ParticleCharge(Int_t pdg) const override
  {
    Warning("ParticleCharge", "Not yet implemented");
    return -1.;
  }

  /// Return life time of the particle specified by pdg.
  Double_t ParticleLifeTime(Int_t pdg) const override
  {
    Warning("ParticleLifeTime", "Not yet implemented");
    return -1.;
  }

  /// Return VMC type of the particle specified by pdg.
  TMCParticleType ParticleMCType(Int_t pdg) const override
  {
    Warning("ParticleMCType", "Not yet implemented");
    return TMCParticleType();
  }
  //
  // ------------------------------------------------
  // methods for step management
  // ------------------------------------------------
  //

  //
  // action methods
  // ------------------------------------------------
  //

  /// Stop the transport of the current particle and skip to the next
  void StopTrack() override
  {
    // Actually, that should never be called but can happen in a replay due to double (original) vs. float (replay), or similar
    mIsTrackStopped = true;
  }

  /// Stop simulation of the current event and skip to the next
  void StopEvent() override
  {
    StopTrack();
    mIsEventStopped = true;
  }

  /// Stop simulation of the current event and set the abort run flag to true
  void StopRun() override
  {
    StopEvent();
    mIsRunStopped = true;
  }

  //
  // set methods
  // ------------------------------------------------
  //

  /// Set the maximum step allowed till the particle is in the current medium
  void SetMaxStep(Double_t) override
  {
    // Warning("SetMaxStep", "Not yet implemented");
  }

  /// Set the maximum number of steps till the particle is in the current medium
  void SetMaxNStep(Int_t) override
  {
    Warning("SetMaxNStep", "Not yet implemented");
  }

  /// Force the decays of particles to be done with Pythia
  /// and not with the Geant routines.
  void SetUserDecay(Int_t pdg) override
  {
    Warning("SetUserDecay", "Not yet implemented");
  }

  /// Force the decay time of the current particle
  void ForceDecayTime(Float_t) override
  {
    Warning("ForceDecayTime", "Not yet implemented");
  }

  //
  // tracking volume(s)
  // ------------------------------------------------
  //

  /// Return the current volume ID and copy number
  Int_t CurrentVolID(Int_t& copyNo) const override;

  /// Return the current volume off upward in the geometrical tree
  /// ID and copy number
  Int_t CurrentVolOffID(Int_t off, Int_t& copyNo) const override;

  /// Return the current volume name
  const char* CurrentVolName() const override;

  /// Return the current volume off upward in the geometrical tree
  /// name and copy number'
  /// if name=0 no name is returned
  const char* CurrentVolOffName(Int_t off) const override;

  /// Return the path in geometry tree for the current volume
  const char* CurrentVolPath() override;

  /// If track is on a geometry boundary, fill the normal vector of the crossing
  /// volume surface and return true, return false otherwise
  Bool_t CurrentBoundaryNormal(Double_t& x, Double_t& y, Double_t& z) const override
  {
    // TODO Not yet implemented, needs some investigation how it is implemente din TGeant3 or TGeant4
    Warning("CurrentBoundaryNormal", "Not yet implemented");
    return kFALSE;
  }

  /// Return the parameters of the current material during transport
  Int_t CurrentMaterial(Float_t& a, Float_t& z, Float_t& dens, Float_t& radl, Float_t& absl) const override
  {
    Warning("CurrentMaterial", "Not yet implemented");
    return -1;
  }

  //// Return the number of the current medium
  Int_t CurrentMedium() const override
  {
    return getMediumId(mCurrentStep->volId);
  }
  // new function (to replace GetMedium() const)

  /// Return the number of the current event
  Int_t CurrentEvent() const override
  {
    return mCurrentEvent;
  }

  /// Computes coordinates xd in daughter reference system
  /// from known coordinates xm in mother reference system.
  /// - xm    coordinates in mother reference system (input)
  /// - xd    coordinates in daughter reference system (output)
  /// - iflag
  ///   - IFLAG = 1  convert coordinates
  ///   - IFLAG = 2  convert direction cosines
  void Gmtod(Float_t* xm, Float_t* xd, Int_t iflag) override;

  /// The same as previous but in double precision
  void Gmtod(Double_t* xm, Double_t* xd, Int_t iflag) override;

  /// Computes coordinates xm in mother reference system
  /// from known coordinates xd in daughter reference system.
  /// - xd    coordinates in daughter reference system (input)
  /// - xm    coordinates in mother reference system (output)
  /// - iflag
  ///   - IFLAG = 1  convert coordinates
  ///   - IFLAG = 2  convert direction cosines
  void Gdtom(Float_t* xd, Float_t* xm, Int_t iflag) override;

  /// The same as previous but in double precision
  void Gdtom(Double_t* xd, Double_t* xm, Int_t iflag) override;

  /// Return the maximum step length in the current medium
  Double_t MaxStep() const override
  {
    return mCurrentStep->maxstep;
  }

  /// Return the maximum number of steps allowed in the current medium
  Int_t GetMaxNStep() const override
  {
    Warning("MaxStep", "Not yet implemented");
    return -1.;
  }

  //
  // get methods
  // tracking particle
  // dynamic properties
  // ------------------------------------------------
  //

  /// Return the current position in the master reference frame of the
  /// track being transported
  void TrackPosition(TLorentzVector& position) const override
  {
    // Time not yet implemented in MCStepLogger
    position.SetXYZT(mCurrentStep->x, mCurrentStep->y, mCurrentStep->z, mCurrentStep->t);
  }

  /// Only return spatial coordinates (as double)
  void TrackPosition(Double_t& x, Double_t& y, Double_t& z) const override
  {
    x = mCurrentStep->x;
    y = mCurrentStep->y;
    z = mCurrentStep->z;
  }

  /// Only return spatial coordinates (as float)
  void TrackPosition(Float_t& x, Float_t& y, Float_t& z) const override
  {
    x = mCurrentStep->x;
    y = mCurrentStep->y;
    z = mCurrentStep->z;
  }

  /// Return the direction and the momentum (GeV/c) of the track
  /// currently being transported
  void TrackMomentum(TLorentzVector& momentum) const override
  {
    momentum.SetXYZT(mCurrentStep->px, mCurrentStep->py, mCurrentStep->pz, mCurrentStep->E);
  }

  /// Return the direction and the momentum (GeV/c) of the track
  /// currently being transported (as double)
  void TrackMomentum(Double_t& px, Double_t& py, Double_t& pz, Double_t& etot) const override
  {
    px = mCurrentStep->px;
    py = mCurrentStep->py;
    pz = mCurrentStep->pz;
    etot = mCurrentStep->E;
  }

  /// Return the direction and the momentum (GeV/c) of the track
  /// currently being transported (as float)
  void TrackMomentum(Float_t& px, Float_t& py, Float_t& pz, Float_t& etot) const override
  {
    px = mCurrentStep->px;
    py = mCurrentStep->py;
    pz = mCurrentStep->pz;
    etot = mCurrentStep->E;
  }

  /// Return the length in centimeters of the current step (in cm)
  Double_t TrackStep() const override
  {
    return mCurrentStep->step;
  }

  /// Return the length of the current track from its origin (in cm)
  Double_t TrackLength() const override
  {
    return mCurrentTrackLength;
  }

  /// Return the current time of flight of the track being transported
  Double_t TrackTime() const override
  {
    return mCurrentStep->t;
  }

  /// Return the energy lost in the current step
  Double_t Edep() const override
  {
    return mCurrentStep->edep;
  }

  /// Return the non-ionising energy lost (NIEL) in the current step
  // TODO To be implemented???
  Double_t NIELEdep() const override
  {
    Warning("NIELEdep", "Not yet implemented");
    return -1.;
  }

  /// Return the current step number
  // TODO To be implemented
  Int_t StepNumber() const override
  {
    Warning("StepNumber", "Not yet implemented");
    return -1;
  }

  /// Get the current weight
  // TODO To be implemented
  Double_t TrackWeight() const override
  {
    Warning("TrackWeight", "Not yet implemented");
    return -1.;
  }

  /// Get the current polarization
  // TODO To be implemented
  void TrackPolarization(Double_t& polX, Double_t& polY, Double_t& polZ) const override
  {
    Warning("TrackPolarization", "Not yet implemented");
  }

  /// Get the current polarization
  // TODO To be implemented
  void TrackPolarization(TVector3& pol) const override
  {
    Warning("TrackPolarization", "Not yet implemented");
  }

  //
  // get methods
  // tracking particle
  // static properties
  // ------------------------------------------------
  //

  /// Return the PDG of the particle transported
  Int_t TrackPid() const override
  {
    return mCurrentLookups->tracktopdg[mCurrentStep->trackID];
  }

  /// Return the charge of the track currently transported
  Double_t TrackCharge() const override
  {
    return mCurrentLookups->tracktocharge[mCurrentStep->trackID];
  }

  /// Return the mass of the track currently transported
  Double_t TrackMass() const override
  {
    return mCurrentLookups->tracktomass[mCurrentStep->trackID];
  }

  /// Return the total energy of the current track
  Double_t Etot() const override
  {
    // TODO make sure E is actually the energy and NOT the energy deposit
    return mCurrentStep->E;
  }

  //
  // get methods - track status
  // ------------------------------------------------
  //

  /// Return true when the track performs the first step
  Bool_t IsNewTrack() const override
  {
    return mCurrentStep->newtrack;
  }

  /// Return true if the track is not at the boundary of the current volume
  Bool_t IsTrackInside() const override
  {
    return mCurrentStep->inside;
  }

  /// Return true if this is the first step of the track in the current volume
  Bool_t IsTrackEntering() const override
  {
    return mCurrentStep->entered;
  }

  /// Return true if this is the last step of the track in the current volume
  Bool_t IsTrackExiting() const override
  {
    return mCurrentStep->exited;
  }

  /// Return true if the track is out of the setup
  Bool_t IsTrackOut() const override
  {
    return mCurrentStep->outside;
  }

  /// Return true if the current particle has disappeared
  /// either because it decayed or because it underwent
  /// an inelastic collision
  Bool_t IsTrackDisappeared() const override
  {
    return mCurrentStep->disappeared;
  }

  /// Return true if the track energy has fallen below the threshold
  Bool_t IsTrackStop() const override
  {
    return mCurrentStep->stopped;
  }

  /// Return true if the current particle is alive and will continue to be
  /// transported
  Bool_t IsTrackAlive() const override
  {
    return !mCurrentStep->stopped;
  }

  //
  // get methods - secondaries
  // ------------------------------------------------
  //

  /// Return the number of secondary particles generated in the current step
  Int_t NSecondaries() const override
  {
    return mCurrentStep->nsecondaries;
  }

  /// Return the parameters of the secondary track number isec produced
  /// in the current step
  void GetSecondary(Int_t isec, Int_t& particleId, TLorentzVector& position, TLorentzVector& momentum) override
  {
    // TODO Not possible atm with MCStepLogger, only secondary info is the number of secondaries produced
    Warning("GetSecondary", "Not yet implemented");
  }

  /// Return the VMC code of the process that has produced the secondary
  /// particles in the current step
  TMCProcess ProdProcess(Int_t isec) const override
  {
    // TODO Not possible atm with MCStepLogger, only secondary info is the number of secondaries produced
    Warning("ProdProcess", "Not yet implemented");
    return TMCProcess();
  }

  /// Return the array of the VMC code of the processes active in the current
  /// step
  Int_t StepProcesses(TArrayI& proc) const override
  {
    // TODO Not possible atm with MCStepLogger, donet have that info
    // Warning("StepProcesses", "Not yet implemented");
    return -1;
  }

  /// Return the information about the transport order needed by the stack
  Bool_t SecondariesAreOrdered() const override
  {
    return kTRUE;
  }

  //
  // ------------------------------------------------
  // Control methods
  // ------------------------------------------------
  //

  /// Initialize MC
  void Init() override;

  /// Initialize MC physics
  void BuildPhysics() override
  {
    Info("BuildPhysics", "Not implemented because not needed at the moment");
  }

  /// Process one event
  // TODO To be implemented
  void ProcessEvent(Int_t eventId) override;

  /// Process one event
  // TODO To be implemented
  void ProcessEvent(Int_t eventId, Bool_t isInterruptible) override
  {
  }

  /// Process one event (backward-compatibility)
  // TODO To be implemented
  void ProcessEvent() override
  {
    // Here we assume that the transport of the next event is desired
    ProcessEvent(mCurrentEvent);
  }

  /// That triggers stopping the transport of the current track without dispatching
  /// to common routines like TVirtualMCApplication::PostTrack() etc.
  void InterruptTrack() override
  {
    Info("InterruptTrack", "Not implemented");
  }

  /// Process one  run and return true if run has finished successfully,
  /// return false in other cases (run aborted by user)
  Bool_t ProcessRun(Int_t nevent) override;

  /// Additional cleanup after a run can be done here (optional)
  // TODO To be implemented
  void TerminateRun() override
  {
    Warning("TerminateRun", "Nothing to be done yet");
  }

  /// Set switches for lego transport
  void InitLego() override
  {
    Warning("InitLego", "Not yet implemented");
  }

  /// (In)Activate collecting TGeo tracks
  void SetCollectTracks(Bool_t collectTracks) override
  {
    Warning("SetCollectTracks", "Not yet implemented");
  }

  /// Return the info if collecting tracks is activated
  Bool_t IsCollectTracks() const override
  {
    Warning("IsCollectTracks", "Not yet implemented");
    return kFALSE;
  }

  /// Return the info if multi-threading is supported/activated
  Bool_t IsMT() const override { return kFALSE; }

  // set filename and treename where to find the steps
  void setStepFilename(const std::string& filename)
  {
    mStepLoggerFilename = filename;
  }
  void setStepTreename(const std::string& treename)
  {
    mStepLoggerTreename = treename;
  }

  // aacept a macro path from where to load a custom function whether or not to keep a step
  void setUserKeepStepMacroPath(const std::string& path)
  {
    mUserKeepStepMacroPath = path;
  }

  void blockSetProcessesCuts(bool value = true)
  {
    mUpdateProcessesCutsBlocked = value;
  }

  void allowStopTrack(bool value = true)
  {
    mAllowStopTrack = value;
  }

  void cutsFromConfig(std::string const& path);

 private:
  // init the run, used to guarantee that both ProcessRun and ProcessEvent
  // work just fine
  bool initRun();

  // map volume IDs assigned by the TGeoManager (replay engine) to those assigned by the engine used in the reference run
  Int_t correctGeoVolID(Int_t geoVolId) const;
  Int_t correctVMCVolID(Int_t vmcVolId) const;

  // That is just a helper function to address some of the
  // interactions with the TGeoManager
  Double_t* makeDoubleArray(Float_t* arrIn, int np) const;

  // Lookup the medium ID for a given volume
  int getMediumId(int volId) const;

  // load cuts and processes for given volume ID
  void loadCurrentCutsAndProcesses(int volId);

  // check whether step is primary
  bool isPrimary(int trackId) const;
  // decide if a step should be kept or not
  bool keepDueToProcesses(const o2::StepInfo& step) const;
  bool keepDueToCuts(const o2::StepInfo& step) const;
  bool keepStep(const o2::StepInfo& step) const;

  // treat a user secondary that was, for instance, pushed to the stack during hit creation
  void transportUserHitSecondary();

  bool startTrack(const o2::StepInfo& step);
  void finishTrack(const o2::StepInfo& nextStep);
  void finishTrack();
  bool stepping();

  // helper function to derive a hash from step properties to seed gRandom
  // On the shoulders of o2::data::Stack
  ULong_t makeHash(const o2::StepInfo& step) const;

  // add process or cut values based on name and value
  template <typename P, typename T, std::size_t N>
  bool insertProcessOrCut(std::vector<std::vector<P>*>& insertInto, const std::array<T, N>& allParamsNames, const std::vector<P>& defaultParams, Int_t mediumId, const char* paramName, P parval, bool forceSet = false)
  {
    auto paramIndex = physics::paramToIndex(allParamsNames, paramName);
    return insertProcessOrCut(insertInto, allParamsNames, defaultParams, mediumId, paramIndex, parval, forceSet);
  }

  // add process or cut values based on name and value
  template <typename P, typename T, std::size_t N>
  bool insertProcessOrCut(std::vector<std::vector<P>*>& insertInto, const std::array<T, N>& allParamsNames, const std::vector<P>& defaultParams, Int_t mediumId, int paramIndex, P parval, bool forceSet = false)
  {
    if (paramIndex < 0) {
      return false;
    }

    if (insertInto.size() <= mediumId) {
      insertInto.resize(mediumId + 1, nullptr);
    }
    auto& currMap = insertInto[mediumId];

    if (!currMap) {
      currMap = new std::vector<P>(defaultParams.begin(), defaultParams.end());
    }
    auto& value = (*currMap)[paramIndex];
    if (value >= 0 && !forceSet) {
      return false;
    }
    value = parval;
    return true;
  }

  template <typename P, typename T, std::size_t N>
  bool insertProcessOrCut(std::vector<P>& insertInto, const std::array<T, N>& allParamsNames, const char* paramName, P parval, bool forceSet = false)
  {
    auto paramIndex = physics::paramToIndex(allParamsNames, paramName);
    return insertProcessOrCut(insertInto, allParamsNames, paramIndex, parval, forceSet);
  }

  template <typename P, typename T, std::size_t N>
  bool insertProcessOrCut(std::vector<P>& insertInto, const std::array<T, N>& allParamsNames, int paramIndex, P parval, bool forceSet = false)
  {
    if (paramIndex < 0) {
      return false;
    }
    auto& value = insertInto[paramIndex];
    if (value >= 0 && !forceSet) {
      return false;
    }
    value = parval;
    return true;
  }

  void printCurrentCuts() const;

 private:
  /// Stack pointers
  TVirtualMCStack* mStack = nullptr;

  /// Run, event, track flags
  /// NOTE These flags could most likely be removed
  bool mIsRunStopped = false;
  bool mIsEventStopped = false;
  bool mIsTrackStopped = false;

  // If we allow for stoopping a track
  bool mAllowStopTrack = true;

  // File and tree name to process
  std::string mStepLoggerFilename;
  std::string mStepLoggerTreename;

  // mapping volume IDs correctly; a reference run might assign other volume IDs than those which would be assigned only by the TGeoManager. Hence, we need to make sure that - whatever engine was used before - we mimic the exact same volume IDs. Otherwise, the logged IDs and those assigned by this engine based on the TGeoManager are not aligned
  std::vector<int> mVolIDMap;
  std::vector<int> mVolIDMapInverse;
  bool mNeedVolIdMapping = false;

  // pointer to opened step file
  TFile* mStepFile = nullptr;

  // branches in the MCStepLogger file
  TBranch* mStepBranch = nullptr;
  TBranch* mLookupBranch = nullptr;

  // flag if initialised, use this flag to either use ProcessRun
  // or to simulate single events with ProcessEvent
  bool mIsInitialised = false;

  // we need to cache the current step and lookups
  // to address all the needs of the TVirtualMC interfaces
  o2::StepLookups* mCurrentLookups = nullptr;
  o2::StepInfo* mCurrentStep = nullptr;

  // keep track of tracks to be skipped
  std::vector<float> mSkipTrack;
  // keep track of track ID assigned by user stack
  std::vector<int> mUserTrackId;

  // current event ID
  int mCurrentEvent = 0;

  // increment the current track length
  float mCurrentTrackLength = 0.;

  // keep track of primary tracks
  int mCurrentPrimaryId = -1;
  // keep track of current track
  int mCurrentTrackId = -1;

  // block in case the using framework should be prevented from setting cuts or processes
  bool mUpdateProcessesCutsBlocked = false;
  // Preliminary process structure
  std::vector<std::vector<Int_t>*> mProcesses;
  // Preliminary cut structure
  std::vector<std::vector<Double_t>*> mCuts;
  // Preliminary process structure for global parameters
  std::vector<Int_t> mProcessesGlobal;
  // Preliminary cut structure for global parameters
  std::vector<Double_t> mCutsGlobal;
  // point to the current map of processes
  std::vector<Int_t>* mCurrentProcesses = nullptr;
  // point to the current map of cuts
  std::vector<Double_t>* mCurrentCuts = nullptr;

  // a user defined cut function, return true if step should be kept, false otherwise
  std::string mUserKeepStepMacroPath;
  user_keep_step_type* mUserKeepStep = nullptr;

  // local pointer to ROOT's geometry manager
  TGeoManager* mGeoManager = nullptr;
  // current TGeoVolume
  TGeoVolume* mCurrentTGeoVolume = nullptr;
  // current TGeoNode
  TGeoNode* mCurrentTGeoNode = nullptr;

  // keep track of number of materials
  int mMaterialCounter = 0;
  // keep track of number of materials
  int mMediumCounter = 0;

  ClassDefOverride(MCReplayEngine, 1);
};
} // end namespace mcreplay

#endif /* MC_REPLAY_ENGINE_H */
