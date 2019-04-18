// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/* This is a simple high-level analysis on steps. It obtains simple accumulated counts of
 *
 * -> how many steps are done in each volume
 * -> how many steps are done per module
 * -> how many steps are done per pdg
 * -> spatial distribution of steps
 *    -> position
 *    -> steps in r-z plane
 * -> number of secondaries per volume
 * -> number of secondaries per module
 * -> number of secondaries per pdg
 *
 * A user cut can be given from the command line
 *
 */

#ifndef SIMPLESTEP_MCANALYSIS_H_
#define SIMPLESTEP_MCANALYSIS_H_

#include "MCStepLogger/MCAnalysis.h"
class TTree;
class TFile;

namespace o2
{
namespace mcstepanalysis
{

class SimpleStepAnalysis : public MCAnalysis
{
 public:
  SimpleStepAnalysis();

 protected:
  /// custom initialization of histograms
  void initialize() override;
  /// custom event loop
  void analyze(const std::vector<StepInfo>* const steps, const std::vector<MagCallInfo>* const magCalls) override;
  /// custom finalizations of produced histograms
  void finalize() override;

 private:
  // accumulated number of steps per module/region
  TH1I* histNStepsPerMod;
  // accumulated number of steps per volume
  TH1I* histNStepsPerVol;
  // accumulated number of steps per pdg
  TH1I* histNStepsPerPDG;
  TH1I* histNStepsPerVolSorted;

  // track production origins per module and volume
  TH1I* histOriginPerMod;
  TH1I* histOriginPerVol;
  TH1I* histOriginPerVolSorted;

  // accumulated number of secondaries produced per volume
  TH1I* histNSecondariesPerVol;
  // accumulated number of secondaries produces per module
  TH1I* histNSecondariesPerMod;
  // accumulated number of secondaries produces per pdg
  TH1I* histNSecondariesPerPDG;

  // energy spectrum of tracks overall as log10f(E)
  TH1I* histTrackEnergySpectrum;
  TH1I* histTrackPDGSpectrum;
  TH1I* histTrackPDGSpectrumSorted;

  /// the production processes per track
  TH1I* histTrackProdProcess;

  // steps in the r-z plane
  TH2D* histRZ;
  // steps in the x-y plane
  TH2D* histXY;

  // keep steps (under cutting for instance)
  TTree *steptree;
  TFile *stepfile;
  
  // pointing to a user cut function 
  cut_function_type* mUserCutFunction = nullptr; //!

  ClassDefNV(SimpleStepAnalysis, 1);
};

} // namespace mcstepanalysis
} // namespace o2
#endif /* BASIC_MCANALYSIS_H_ */
