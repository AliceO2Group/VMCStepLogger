// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//  @file   MCStepLoggerImpl.cxx
//  @author Sandro Wenzel
//  @since  2017-06-29
//  @brief  A logging service for MCSteps (hooking into Stepping of TVirtualMCApplication's)

#include "MCStepLogger/StepInfo.h"
#include "MCStepLogger/MetaInfo.h"
#include "MCStepLogger/StepLogger.h"
#include "MCStepLogger/FieldLogger.h"
#include "MCStepLogger/StepLoggerUtilities.h"
#include <TBranch.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <TVirtualMCApplication.h>
#include <TVirtualMagField.h>
#include <sstream>

#include <dlfcn.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

// the global logging instances (in anonymous namespace)
// pointers to dissallow construction at each library load
// NOTE in case something breaks
namespace o2 {
StepLogger* logger;
FieldLogger* fieldlogger;
}

// a generic function that can dispatch to the original method of a TVirtualMCApplication
extern "C" void dispatchOriginal(TVirtualMCApplication* app, char const* libname, char const* origFunctionName)
{
  typedef void (TVirtualMCApplication::*StepMethodType)();
  o2::mcsteploggerutilities::dispatchOriginalKernel<TVirtualMCApplication, StepMethodType>(app, libname, origFunctionName);
}

// a generic function that can dispatch to the original method of a TVirtualMagField
extern "C" void dispatchOriginalField(TVirtualMagField* field, char const* libname, char const* origFunctionName,
                                      const double x[3], double* B)
{
  typedef void (TVirtualMagField::*MethodType)(const double[3], double*);
  o2::mcsteploggerutilities::dispatchOriginalKernel<TVirtualMagField, MethodType>(field, libname, origFunctionName, x, B);
}

extern "C" void performLogging(TVirtualMCApplication* app)
{
  static TVirtualMC* mc = TVirtualMC::GetMC();
  o2::logger->addStep(mc);
}

extern "C" void logField(double const* p, double const* b)
{
  static TVirtualMC* mc = TVirtualMC::GetMC();
  o2::fieldlogger->addStep(mc, p, b);
}

extern "C" void initLogger()
{
  // init TFile for logging output
  o2::mcsteploggerutilities::initTFile();
  // initializes the logging instances
  o2::logger = &o2::StepLogger::Instance();
  o2::fieldlogger = &o2::FieldLogger::Instance();
}

extern "C" void flushLog()
{
  std::cerr << "[MCLOGGER:] START FLUSHING ----\n";
  o2::logger->flush();
  o2::fieldlogger->flush();
  std::cerr << "[MCLOGGER:] END FLUSHING ----\n";
}
