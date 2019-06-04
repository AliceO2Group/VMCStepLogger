// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


#include "MCStepLogger/StepLoggerUtilities.h"
#include "MCStepLogger/StepInfo.h"

namespace o2
{
namespace mcsteploggerutilities
{

const char* getLogFileName()
{
  if (const char* f = std::getenv("MCSTEPLOG_OUTFILE")) {
    return f;
  } else {
    return "MCStepLoggerOutput.root";
  }
}

const char* getVolMapFile()
{
  if (const char* f = std::getenv("MCSTEPLOG_VOLMAPFILE")) {
    return f;
  } else {
    return "MCStepLoggerVolMap.dat";
  }
}

const char* getSensitiveVolFile()
{
  if (const char* f = std::getenv("MCSTEPLOG_SENSVOLFILE")) {
    return f;
  } else {
    return "MCStepLoggerSenVol.dat";
  }
}

// initializes a mapping from volumename to detector
// used for step resolution to detectors
void initVolumeMap()
{
  auto volmap = new std::map<std::string, std::string>;
  // open for reading or fail
  std::ifstream ifs;
  auto f = getVolMapFile();
  std::cerr << "[MCLOGGER:] TRYING TO READ VOLUMEMAPS FROM " << f << "\n";
  ifs.open(f);
  if (ifs.is_open()) {
    std::string line;
    while (std::getline(ifs, line)) {
      std::istringstream ss(line);
      std::string token;
      // split the line into key + value
      int counter = 0;
      std::string keyvalue[2] = { "NULL", "NULL" };
      while (counter < 2 && std::getline(ss, token, ':')) {
        if (!token.empty()) {
          keyvalue[counter] = token;
          counter++;
        }
      }
      // put into map
      volmap->insert({ keyvalue[0], keyvalue[1] });
    }
    ifs.close();
    StepInfo::volnametomodulemap = volmap;
  } else {
    std::cerr << "[MCLOGGER:] VOLUMEMAPSFILE NOT FILE\n";
    StepInfo::volnametomodulemap = nullptr;
  }
}

void initTFile()
{
  if (!std::getenv("MCSTEPLOG_TTREE")) {
    return;
  }
  TFile* f = new TFile(getLogFileName(), "RECREATE");
  f->Close();
  delete f;
}

} // namespace mcsteploggerutilities
} // namespace o2
