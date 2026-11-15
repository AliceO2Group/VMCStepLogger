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

// A stand-in for o2::field::MagneticField that returns a fixed field and counts its calls

#include <TVirtualMagField.h>

namespace o2
{
namespace field
{
class MagneticField : public TVirtualMagField
{
 public:
  MagneticField() : TVirtualMagField("fake") {}
  void Field(const double* point, double* bField) override;
};
} // namespace field
} // namespace o2

namespace
{
int nField = 0;
}

void o2::field::MagneticField::Field(const double* point, double* bField)
{
  ++nField;
  bField[0] = 0.;
  bField[1] = 0.;
  bField[2] = 5.;
}

extern "C" TVirtualMagField* makeFakeField() { return new o2::field::MagneticField(); }

extern "C" int fakeFieldCount() { return nField; }
