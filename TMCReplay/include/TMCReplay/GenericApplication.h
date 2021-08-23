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

#ifndef TMC_REPLAY_GENERIC_APPLICATION_H
#define TMC_REPLAY_GENERIC_APPLICATION_H

#include <string>

#include <TVirtualMCApplication.h>

class TGeoManager;

namespace tmcreplay
{
class GenericStack;
}

namespace tmcreplay
{
// K is type of user kernel
class GenericApplication : public TVirtualMCApplication
{

 public:
  GenericApplication(const std::string& geoFilename, const std::string& geoKeyname);

  /// For now just default destructor
  virtual ~GenericApplication() = default;

  // ConstructGeometry is the only implementation of a pure virtual method that does something at the moment

  virtual void ConstructGeometry() override;

  // Implement all pure virtual methods, doing nothinig at the moment

  virtual void InitGeometry() override { ; }

  virtual void GeneratePrimaries() override { ; }

  virtual void BeginEvent() override;

  virtual void BeginPrimary() override { ; }

  virtual void PreTrack() override { ; }

  virtual void Stepping() override
  {
    ;
  }

  virtual void PostTrack() override { ; }

  virtual void FinishPrimary() override { ; }

  virtual void FinishEvent() override { ; }

  void setStack(GenericStack* stack)
  {
    fStack = stack;
  }

 private:
  // Filename where geometry can be found
  std::string fGeoFilename;
  // Keyname under which geometry can be found inside the above file
  std::string fGeoKeyname;
  // local pointer to ROOT's geometry manager
  TGeoManager* fGeoManager;
  // stack
  GenericStack* fStack;

  ClassDefOverride(GenericApplication, 1);
};
} // end namespace tmcreplay

#endif /* TMC_REPLAY_GENERIC_APPLICATION_H */
