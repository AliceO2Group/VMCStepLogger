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

// Drives the fake application and field under the preloaded interceptor and checks the originals are reached

#include <TFile.h>
#include <TTree.h>
#include <TVirtualMCApplication.h>
#include <TVirtualMagField.h>
#include <dlfcn.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#ifdef FAKE_LINKED
extern "C" TVirtualMCApplication* makeFakeApplication();
extern "C" void fakeApplicationCounts(int*, int*, int*);
extern "C" TVirtualMagField* makeFakeField();
extern "C" int fakeFieldCount();
#endif

namespace
{
int failures = 0;

void check(bool ok, std::string const& what)
{
  std::cout << (ok ? "OK   " : "FAIL ") << what << "\n";
  failures += ok ? 0 : 1;
}

// the global definition of symbol must be the interceptor's
void checkInterposed(char const* symbol)
{
  Dl_info info{};
  void* address = dlsym(RTLD_DEFAULT, symbol);
  bool ok = address && dladdr(address, &info) && info.dli_fname && std::strstr(info.dli_fname, "MCStepLoggerInterceptSteps");
  check(ok, std::string("interceptor defines ") + symbol);
}

template <typename F>
F lookup(void* handle, char const* name)
{
  void* address = dlsym(handle, name);
  if (!address) {
    std::cerr << "cannot find " << name << ": " << dlerror() << "\n";
    std::exit(2);
  }
  return reinterpret_cast<F>(address);
}

void* open(char const* path, int scope)
{
  void* handle = dlopen(path, RTLD_NOW | scope);
  if (!handle) {
    std::cerr << "cannot open " << path << ": " << dlerror() << "\n";
    std::exit(2);
  }
  return handle;
}
} // namespace

int main(int argc, char* argv[])
{
  char const* outfile = std::getenv("MCSTEPLOG_OUTFILE");
  if (!outfile || !std::getenv("MCSTEPLOG_TTREE")) {
    std::cerr << "MCSTEPLOG_OUTFILE and MCSTEPLOG_TTREE must be set\n";
    return 2;
  }
  std::remove(outfile);

  using MakeApp = TVirtualMCApplication* (*)();
  using AppCounts = void (*)(int*, int*, int*);
  using MakeField = TVirtualMagField* (*)();
  using FieldCount = int (*)();
  MakeApp makeApp;
  AppCounts appCounts;
  MakeField makeField;
  FieldCount fieldCount;
#ifdef FAKE_LINKED
  makeApp = makeFakeApplication;
  appCounts = fakeApplicationCounts;
  makeField = makeFakeField;
  fieldCount = fakeFieldCount;
#else
  // usage: driver global|local <application library> <field library>
  if (argc != 4) {
    std::cerr << "usage: " << argv[0] << " global|local <application library> <field library>\n";
    return 2;
  }
  int scope = std::strcmp(argv[1], "local") == 0 ? RTLD_LOCAL : RTLD_GLOBAL;
  void* appLib = open(argv[2], scope);
  void* fieldLib = open(argv[3], scope);
  makeApp = lookup<MakeApp>(appLib, "makeFakeApplication");
  appCounts = lookup<AppCounts>(appLib, "fakeApplicationCounts");
  makeField = lookup<MakeField>(fieldLib, "makeFakeField");
  fieldCount = lookup<FieldCount>(fieldLib, "fakeFieldCount");
#endif

  checkInterposed("_ZN17FairMCApplication8SteppingEv");
  checkInterposed("_ZN17FairMCApplication11FinishEventEv");
  checkInterposed("_ZN17FairMCApplication17ConstructGeometryEv");
  checkInterposed("_ZN2o25field13MagneticField5FieldEPKdPd");

  TVirtualMCApplication* app = makeApp();
  TVirtualMagField* field = makeField();
  app->ConstructGeometry();
  double point[3] = {1., 2., 3.};
  double b[3] = {0., 0., 0.};
  field->Field(point, b);
  field->Field(point, b);
  for (int i = 0; i < 3; ++i) {
    app->Stepping();
  }
  app->FinishEvent();

  int nConstruct = 0, nStepping = 0, nFinish = 0;
  appCounts(&nConstruct, &nStepping, &nFinish);
  check(nConstruct == 1, "original ConstructGeometry reached once");
  check(nStepping == 3, "original Stepping reached three times");
  check(nFinish == 1, "original FinishEvent reached once");
  check(fieldCount() == 2 && b[2] == 5., "original Field reached twice and its value returned");

  // without a VMC engine no step is recorded, but the logger writes this file only when its
  // ConstructGeometry and FinishEvent hooks ran
  TFile f(outfile, "READ");
  auto tree = f.IsZombie() ? nullptr : f.Get<TTree>("StepLoggerTree");
  check(tree && tree->GetBranch("Steps") && tree->GetBranch("Calls"), "logger set up and flushed its output");

  return failures == 0 ? 0 : 1;
}
