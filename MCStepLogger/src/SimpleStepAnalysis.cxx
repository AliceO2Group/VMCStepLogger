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
#include <cmath>

#include "MCStepLogger/SimpleStepAnalysis.h"
#include "MCStepLogger/MCAnalysisUtilities.h"
#include "TInterpreter.h"
#include "TInterpreterValue.h"
#include "TFile.h"
#include "TTree.h"
#include "TDatabasePDG.h"

ClassImp(o2::mcstepanalysis::SimpleStepAnalysis);

using namespace o2::mcstepanalysis;

SimpleStepAnalysis::SimpleStepAnalysis()
  : MCAnalysis("SimpleStepAnalysis")
{
}

void SimpleStepAnalysis::initialize()
{
  // accumulated number of steps per module/region
  histNStepsPerMod = getHistogram<TH1I>("nStepsPerMod", 1, 2., 1.);
  // accumulated number of steps per volume
  histNStepsPerVol = getHistogram<TH1I>("nStepsPerVol", 1, 2., 1.);

  histOriginPerMod = getHistogram<TH1I>("OriginsPerMod", 1, 2., 1.);
  histOriginPerVol = getHistogram<TH1I>("OriginsPerVol", 1, 2., 1.);
  histOriginPerVolSorted = getHistogram<TH1I>("OriginsPerVolSorted", 1, 2., 1.);

  // accumulated number of steps per pdg
  histNStepsPerPDG = getHistogram<TH1I>("nStepsPerPDG", 1, 2., 1.);
  histNStepsPerVolSorted = getHistogram<TH1I>("nStepsPerVolSorted", 1, 2., 1.);

  // accumulated number of secondaries produced per volume
  histNSecondariesPerVol = getHistogram<TH1I>("nSecondariesPerVol", 1, 2., 1.);
  // accumulated number of secondaries produces per module
  histNSecondariesPerMod = getHistogram<TH1I>("nSecondariesPerMod", 1, 2., 1.);
  histNSecondariesPerMod->Sumw2(false);
  // accumulated number of secondaries produces per pdg
  histNSecondariesPerPDG = getHistogram<TH1I>("nSecondariesPerPDG", 1, 2., 1.);

  histTrackEnergySpectrum = getHistogram<TH1I>("trackEnergySpectrum", 400, -10, 4.);
  histTrackPDGSpectrum = getHistogram<TH1I>("trackPDGSpectrum", 1, 2., 1.);
  histTrackPDGSpectrumSorted = getHistogram<TH1I>("trackPDGSpectrumSorted", 1, 2., 1.);
  histTrackProdProcess = getHistogram<TH1I>("trackProdProcess", 1, 2., 1.);

  // steps in the r-z plane
  histRZ = getHistogram<TH2D>("RZOccupancy", 200, -3000., 3000., 200, 0., 3000.);
  // steps in x-y plane
  histXY = getHistogram<TH2D>("XYOccupancy", 200, -3000., 3000., 200, -3000., 3000.);

  // init runtime user cut
  // thanks to discussions with Philippe Canal, Fermilab
  auto cutcondition = getenv("MCSTEPCUT");
  if (cutcondition) {
    std::string expr =
      "#include \"StepInfo.h\"\n #include <cmath>\n bool user_cut(const o2::StepInfo &step, \
                        const std::string &volname,                                       \
                        const std::string &modname,                                       \
                        int pdg, o2::StepLookups* lookup) {";
    expr += std::string(cutcondition);
    expr += std::string("}");
    auto installpath = getenv("MCSTEPLOGGER_ROOT");
    if (installpath) {
      auto includepath = std::string(installpath) + std::string("/include/MCStepLogger");
      std::cout << "Using include path " << includepath << "\n";
      gInterpreter->AddIncludePath(includepath.c_str());
    } else {
      std::cerr << "Could not set path to Steplogger headers; Just-in-time compilation might fail"
                << "\n";
    }
    gInterpreter->Declare(expr.c_str());
    TInterpreterValue* v = gInterpreter->CreateTemporary();
    // std::unique_ptr<TInterpreterValue> v = gInterpreter->MakeInterpreterValue();
    gInterpreter->Evaluate("user_cut", *v);
    mUserCutFunction = (cut_function_type*)v->GetValAddr();
  }

  if (getenv("KEEPSTEPS")) {
    steptree = new TTree("Steps", "Steps");
  }
}

void SimpleStepAnalysis::analyze(const std::vector<StepInfo>* const steps, const std::vector<MagCallInfo>* const magCalls)
{
  static auto pdgdatabase = new TDatabasePDG;
  static StepInfo const* stepptr;
  static TBranch* branch;
  static bool keepsteps = getenv("KEEPSTEPS") != nullptr;
  if (keepsteps && !branch) {
    branch = steptree->Branch("Steps", &stepptr);
    stepfile = new TFile("Steps.root", "RECREATE");
    steptree->SetDirectory(stepfile);
  }

  // to store the volume name
  std::string volName;
  // to store the module name
  std::string modName;

  // total number of steps in this event
  int nSteps = 0;
  // to store particle ID
  int pdgId = 0;
  int nCutSteps = 0;

  int oldTrackID = -1; // to notice when a track changes

  // loop over all steps in an event
  for (const auto& step : *steps) {

    // prepare for PDG ids and volume names
    mAnalysisManager->getLookupPDG(step.trackID, pdgId);
    mAnalysisManager->getLookupVolName(step.volId, volName);
    mAnalysisManager->getLookupModName(step.volId, modName);

    // first increment the total number of steps over all events
    nSteps++;

    // apply user defined cut -- if any
    if (mUserCutFunction && !(*mUserCutFunction)(step, volName, modName, pdgId, mAnalysisManager->getLookups())) {
      nCutSteps++;
      continue;
    }

    int currentTrackID = step.trackID;
    bool newtrack = (currentTrackID != oldTrackID);
    if (newtrack) {
      oldTrackID = currentTrackID;
    }

    if (keepsteps) {
      stepptr = &step;
      steptree->Fill();
    }

    auto pdgparticle = pdgdatabase->GetParticle(pdgId);
    std::string pdgasstring(pdgparticle ? pdgparticle->GetName() : std::to_string(pdgId));

    if (newtrack) {
      histTrackEnergySpectrum->Fill(log10f(step.E));
      histTrackPDGSpectrum->Fill(pdgasstring.c_str(), 1);
      histTrackProdProcess->Fill(step.getProdProcessAsString(), 1);

      auto originid = mAnalysisManager->getLookups()->trackorigin[step.trackID];
      std::string originVolName;
      // to store the module name
      std::string originModName;
      mAnalysisManager->getLookupVolName(originid, originVolName);
      mAnalysisManager->getLookupModName(originid, originModName);

      histOriginPerMod->Fill(originModName.c_str(), 1);
      histOriginPerVol->Fill(originVolName.c_str(), 1);
    }

    // record number of steps per module
    histNStepsPerMod->Fill(modName.c_str(), 1);
    // record number of steps per volume
    histNStepsPerVol->Fill(volName.c_str(), 1);
    // record number of steps per volume
    histNStepsPerPDG->Fill(pdgasstring.c_str(), 1);

    histNSecondariesPerVol->Fill(volName.c_str(), step.nsecondaries);
    histNSecondariesPerMod->Fill(modName.c_str(), step.nsecondaries);
    histNSecondariesPerPDG->Fill(pdgasstring.c_str(), step.nsecondaries);

    histRZ->Fill(step.z, std::sqrt(step.x * step.x + step.y * step.y));
    histXY->Fill(step.x, step.y);
  }
}

void SimpleStepAnalysis::finalize()
{
  *histNStepsPerVolSorted = *histNStepsPerVol;
  histNStepsPerVolSorted->SetName("nStepsPerVolSorted");
  histNStepsPerVolSorted->LabelsOption(">", "X");

  *histOriginPerVolSorted = *histOriginPerVol;
  histOriginPerVolSorted->SetName("OriginPerVolSorted");
  histOriginPerVolSorted->LabelsOption(">", "X");
  histOriginPerVolSorted->SetBins(30, 0, 30);

  std::cerr << "MOD have " << histNStepsPerMod->GetEntries() << " entries \n";

  *histTrackPDGSpectrumSorted = *histTrackPDGSpectrum;
  histTrackPDGSpectrumSorted->SetName("trackPDGSpectrumSorted");
  histTrackPDGSpectrumSorted->LabelsOption(">", "X");
  histTrackPDGSpectrumSorted->SetBins(10, 0, 10);

  // sortit
  // histNStepsPerVolSorted->LabelsOption(">", "X");

  histNStepsPerVolSorted->SetBins(30, 0, 30);
  histNStepsPerMod->LabelsOption(">", "X");
  histNStepsPerMod->SetBins(30, 0, 30);

  histNSecondariesPerMod->LabelsOption(">", "X");
  histNSecondariesPerVol->LabelsOption(">", "X");
  histNSecondariesPerVol->SetBins(30, 0, 30);

  if (getenv("KEEPSTEPS")) {
    std::cout << "Writing step tree\n";
    steptree->Write();
    stepfile->Close();
  }
}
