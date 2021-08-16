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

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>
#include <vector>
#include <string>

#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

#include <MCStepLogger/MCAnalysisFileWrapper.h>
#include <MCStepLogger/MetaInfo.h>
#endif

// Produce overlay histograms for some observables and pass to canvases
void plotHistograms(const char* analysisFilepath, const char* outputDir = "./")
{
  o2::mcstepanalysis::MCAnalysisFileWrapper wrapper;
  wrapper.read(analysisFilepath);

  for (int i = 0; i < wrapper.nHistograms(); i++) {
    TH1& histo = wrapper.getHistogram(i);
    TCanvas* canvas = new TCanvas(histo.GetName(), "", 600, 500);
    TLegend* legend = new TLegend(0.1, 0.7, 0.3, 0.9);

    histo.Draw("hist");
    legend->AddEntry(&histo, wrapper.getAnalysisMetaInfo().label.c_str());
    legend->Draw("same");

    std::string outputDirStr(outputDir);
    if (outputDirStr.back() != '/') {
      outputDirStr += '/';
    }

    canvas->SaveAs(std::string(outputDirStr + histo.GetName() + ".eps").c_str());
    canvas->SaveAs(std::string(outputDirStr + histo.GetName() + ".png").c_str());
    delete canvas;
  }
}

// Produce overlay histograms for some observables and pass to canvases
void plotHistogram(const char* analysisFilepath, const char* histoName,
                   const char* outputDir = "./")
{
  o2::mcstepanalysis::MCAnalysisFileWrapper wrapper;
  wrapper.read(analysisFilepath);

  if (!wrapper.hasHistogram(histoName)) {
    std::cerr << "Histogram " << histoName << " does not exist.\n";
    return;
  }

  TH1& histo = wrapper.getHistogram(histoName);
  TCanvas* canvas = new TCanvas(histo.GetName(), "", 600, 500);
  TLegend* legend = new TLegend(0.1, 0.7, 0.3, 0.9);

  histo.Draw("hist");
  legend->AddEntry(&histo, wrapper.getAnalysisMetaInfo().label.c_str());
  legend->Draw("same");

  std::string outputDirStr(outputDir);
  if (outputDirStr.back() != '/') {
    outputDirStr += '/';
  }

  canvas->SaveAs(std::string(outputDirStr + histoName + ".eps").c_str());
  canvas->SaveAs(std::string(outputDirStr + histoName + ".png").c_str());
  delete canvas;
}
