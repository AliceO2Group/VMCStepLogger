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

//  @file   MCStepInterceptor.cxx
//  @author Sandro Wenzel
//  @since  2017-06-29
//  @brief  A LD_PRELOAD logger hooking into Stepping of TVirtualMCApplication's

class TVirtualMCApplication;
class TVirtualMagField;

// (re)declare symbols to be able to hook into them
#define DECLARE_INTERCEPT_SYMBOLS(APP) \
  class APP                            \
  {                                    \
   public:                             \
    void Stepping();                   \
    void FinishEvent();                \
    void ConstructGeometry();          \
  };

DECLARE_INTERCEPT_SYMBOLS(FairMCApplication)
DECLARE_INTERCEPT_SYMBOLS(AliMC)
namespace mcreplay
{
DECLARE_INTERCEPT_SYMBOLS(MCReplayGenericApplication)
}

// same for field
#define DECLARE_INTERCEPT_FIELD_SYMBOLS(FIELD)       \
  class FIELD                                        \
  {                                                  \
   public:                                           \
    void Field(const double* point, double* bField); \
  };

namespace o2
{
namespace field
{
DECLARE_INTERCEPT_FIELD_SYMBOLS(MagneticField);
}
} // namespace o2
DECLARE_INTERCEPT_FIELD_SYMBOLS(AliMagF);

extern "C" void performLogging(TVirtualMCApplication*);
extern "C" void logField(const double*, const double*);
extern "C" void dispatchOriginal(TVirtualMCApplication*, char const* libname, char const*);
extern "C" void dispatchOriginalField(TVirtualMagField*, char const* libname, char const*, const double x[3],
                                      double* B);
extern "C" void flushLog();
extern "C" void initLogger();

#define INTERCEPT_STEPPING(APP, LIB, SYMBOL)                       \
  void APP::Stepping()                                             \
  {                                                                \
    auto baseptr = reinterpret_cast<TVirtualMCApplication*>(this); \
    performLogging(baseptr);                                       \
    dispatchOriginal(baseptr, LIB, SYMBOL);                        \
  }

#define INTERCEPT_FINISHEVENT(APP, LIB, SYMBOL)                    \
  void APP::FinishEvent()                                          \
  {                                                                \
    auto baseptr = reinterpret_cast<TVirtualMCApplication*>(this); \
    flushLog();                                                    \
    dispatchOriginal(baseptr, LIB, SYMBOL);                        \
  }

// we use the ConstructGeometry hook to setup the logger
#define INTERCEPT_GEOMETRYINIT(APP, LIB, SYMBOL)                   \
  void APP::ConstructGeometry()                                    \
  {                                                                \
    auto baseptr = reinterpret_cast<TVirtualMCApplication*>(this); \
    dispatchOriginal(baseptr, LIB, SYMBOL);                        \
    initLogger();                                                  \
  }

// the runtime will now dispatch to these functions due to LD_PRELOAD
INTERCEPT_STEPPING(FairMCApplication, "libBase", "_ZN17FairMCApplication8SteppingEv")
INTERCEPT_STEPPING(AliMC, "libSTEER", "_ZN5AliMC8SteppingEv")
INTERCEPT_STEPPING(mcreplay::MCReplayGenericApplication, "libMCReplayCore", "_ZN8mcreplay26MCReplayGenericApplication8SteppingEv")

INTERCEPT_FINISHEVENT(FairMCApplication, "libBase", "_ZN17FairMCApplication11FinishEventEv")
INTERCEPT_FINISHEVENT(AliMC, "libSTEER", "_ZN5AliMC11FinishEventEv")
INTERCEPT_FINISHEVENT(mcreplay::MCReplayGenericApplication, "libMCReplayCore", "_ZN8mcreplay26MCReplayGenericApplication11FinishEventEv")

INTERCEPT_GEOMETRYINIT(FairMCApplication, "libBase", "_ZN17FairMCApplication17ConstructGeometryEv")
INTERCEPT_GEOMETRYINIT(AliMC, "libSTEER", "_ZN5AliMC17ConstructGeometryEv")
INTERCEPT_GEOMETRYINIT(mcreplay::MCReplayGenericApplication, "libMCReplayCore", "_ZN8mcreplay26MCReplayGenericApplication17ConstructGeometryEv")

#define INTERCEPT_FIELD(FIELD, LIB, SYMBOL)                     \
  void FIELD::Field(const double* point, double* bField)        \
  {                                                             \
    auto baseptr = reinterpret_cast<TVirtualMagField*>(this);   \
    dispatchOriginalField(baseptr, LIB, SYMBOL, point, bField); \
    logField(point, bField);                                    \
  }

namespace o2
{
namespace field
{
INTERCEPT_FIELD(MagneticField, "libO2Field", "_ZN2o25field13MagneticField5FieldEPKdPd");
}
} // namespace o2
// for AliRoot
INTERCEPT_FIELD(AliMagF, "libSTEERBase", "_ZN7AliMagF5FieldEPKdPd")
