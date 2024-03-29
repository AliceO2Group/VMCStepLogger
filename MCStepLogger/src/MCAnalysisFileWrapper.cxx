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

#include "TSystem.h" // to check for and create directories

#include "MCStepLogger/MCAnalysisFileWrapper.h"
#include "MCStepLogger/ROOTIOUtilities.h"

ClassImp(o2::mcstepanalysis::MCAnalysisFileWrapper);

using namespace o2::mcstepanalysis;

MCAnalysisFileWrapper::MCAnalysisFileWrapper()
  : mInputFilepath(""), mAnalysisMetaInfo(MCAnalysisMetaInfo()), mHasChanged(false)
{
  mHistograms.clear();
}

bool MCAnalysisFileWrapper::isSane() const
{
  bool sane = true;
  // so far only check number of histograms vs. expected number of histograms
  if (mAnalysisMetaInfo.nHistograms != mHistograms.size()) {
    std::cerr << "ERROR: Histograms are corrupted: found " << mHistograms.size() << " but " << mAnalysisMetaInfo.nHistograms << " expected\n";
    sane = false;
  }
  return sane;
}

bool MCAnalysisFileWrapper::read(const std::string& filepath)
{
  ROOTIOUtilities rootutil(filepath);
  // look for and read meta info of step logger and analysis run
  if (!rootutil.hasObject(defaults::mcAnalysisMetaInfoName)) {
    rootutil.close();
    return false;
  }
  rootutil.readObject(mAnalysisMetaInfo, defaults::mcAnalysisMetaInfoName);

  // try to recover histograms
  TH1* histoRecover = nullptr;

  if (!rootutil.changeToTDirectory(defaults::mcAnalysisObjectsDirName)) {
    rootutil.close();
    return false;
  }
  while (true) {
    rootutil.readObject(histoRecover);
    // assuming that all objects have been read
    if (!histoRecover) {
      break;
    }
    TH1* histo = dynamic_cast<TH1*>(histoRecover->Clone());
    histo->SetDirectory(0);
    mHistograms.push_back(std::shared_ptr<TH1>(histo));
  }
  rootutil.close();
  isSane();
  return true;
}

void MCAnalysisFileWrapper::write(const std::string& filedir) const
{
  if (!isSane()) {
    std::cerr << "ERROR: Analysis file cannot be written, see above.\n";
    return;
  }
  const std::string outputDir = filedir + "/" + mAnalysisMetaInfo.analysisName;
  if (!createDirectory(outputDir)) {
    std::cerr << "ERROR: Directory " << outputDir << " could not be created for analysis " << mAnalysisMetaInfo.analysisName << ". Skip...\n";
    return;
  }
  ROOTIOUtilities rootutil(outputDir + "/Analysis.root", ETFileMode::kRECREATE);

  std::cerr << "INFO: Save histograms of analysis " << mAnalysisMetaInfo.analysisName << " at " << filedir << std::endl;
  rootutil.writeObject(&mAnalysisMetaInfo, defaults::mcAnalysisMetaInfoName);
  rootutil.changeToTDirectory(defaults::mcAnalysisObjectsDirName);
  for (const auto& h : mHistograms) {
    rootutil.writeObject(h.get());
  }
  rootutil.close();
}

TH1* MCAnalysisFileWrapper::findHistogram(const std::string& name)
{
  for (auto& h : mHistograms) {
    if (name.compare(h->GetName()) == 0) {
      return h.get();
    }
  }
  return nullptr;
}

bool MCAnalysisFileWrapper::hasHistogram(const std::string& name)
{
  return (findHistogram(name) != nullptr);
}

bool MCAnalysisFileWrapper::createDirectory(const std::string& dir)
{
  gSystem->mkdir(dir.c_str(), true);
  // according to documentation returns false if possible to access
  return (gSystem->AccessPathName(dir.c_str()) == 0);
}

MCAnalysisMetaInfo& MCAnalysisFileWrapper::getAnalysisMetaInfo()
{
  return mAnalysisMetaInfo;
}

int MCAnalysisFileWrapper::nHistograms() const
{
  return mHistograms.size();
}

void MCAnalysisFileWrapper::printAnalysisMetaInfo() const
{
  std::cerr << "INFO: Meta info of analysis file\n";
  mAnalysisMetaInfo.print();
}

void MCAnalysisFileWrapper::printHistogramInfo(const std::string& option) const
{
  if (mHistograms.empty()) {
    return;
  }
  std::cerr << "INFO: Histograms of analysis file\n";
  for (auto& h : mHistograms) {
    std::cerr << "Histogram of class " << h->ClassName() << std::endl;
    h->Print(option.c_str());
  }
}
