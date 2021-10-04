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

#ifndef O2_STEPINFO
#define O2_STEPINFO

#include <Rtypes.h>
#include <chrono>
#include <iostream>
#include <map>
#include <vector>

class TVirtualMC;
class TGeoVolume;
class TGeoMedium;

namespace o2
{

struct VolInfoContainer {
  VolInfoContainer() = default;

  // keeps info about volumes (might exist somewhere else already)
  // essentially a mapping from volumeID x copyNo to TGeoVolumes
  std::vector<std::vector<TGeoVolume const*>*> volumes; // sparse container

  void insert(int id, int copyNo, TGeoVolume const* vol)
  {
    if (volumes.size() <= id) {
      volumes.resize(id + 1, nullptr);
    }
    if (volumes[id] == nullptr) {
      volumes[id] = new std::vector<TGeoVolume const*>;
    }
    if (volumes[id]->size() <= copyNo) {
      volumes[id]->resize(copyNo + 1, nullptr);
    }
    (*volumes[id])[copyNo] = vol;
  }

  TGeoVolume const* get(int id, int copy) const { return (*volumes[id])[copy]; }
  ClassDefNV(VolInfoContainer, 1);
};

// LookupStructures to translate some step information such
// as volumeid, to readable names
struct StepLookups {
  // using pointers to allow "nullptr==unknown" and faster query
  std::vector<std::string*> volidtovolname;
  std::vector<std::string*> volidtomodule;
  std::vector<std::string*> volidtomedium;
  std::vector<int> geovolidtovmcvolid;
  std::vector<bool> volidtoissensitive; // keep track of which volume is sensitive
  std::vector<int> tracktopdg;
  std::vector<float> tracktocharge;
  std::vector<float> tracktomass;
  std::vector<int> tracktoparent; // when parent is -1 we mean primary
  std::vector<int> stepcounterpertrack;

  std::vector<int> trackorigin;        // the volumeID where the track originates
  std::vector<float> tracktoenergy;    // mapping of trackID to the track energy
  std::vector<bool> crossedboundary;   // if track every crossed a geometry boundary
  std::vector<bool> producedsecondary; // if track ever produced a secondary

  void insertVolName(int index, std::string const& s) { insertValueAt(index, s, volidtovolname); }
  void insertModuleName(int index, std::string const& s) { insertValueAt(index, s, volidtomodule); }
  void insertMediumName(int index, std::string const& s) { insertValueAt(index, s, volidtomedium); }
  std::string* getModuleAt(int index) const
  {
    if (index >= volidtomodule.size())
      return nullptr;
    return volidtomodule[index];
  }

  void insertPDG(int trackindex, int pdg)
  {
    constexpr int INVALIDPDG = 0;
    if (trackindex >= tracktopdg.size()) {
      tracktopdg.resize(trackindex + 1, INVALIDPDG);
    }
    auto prev = tracktopdg[trackindex];
    if (prev != INVALIDPDG && prev != pdg) {
      std::cerr << "Warning: Seeing more than one pdg (prev: " << prev << " curr: " << pdg << ") for same trackID " << trackindex << "\n";
    }
    tracktopdg[trackindex] = pdg;
  }

  void setTrackCharge(int trackindex, float c)
  {
    if (trackindex >= tracktocharge.size()) {
      tracktocharge.resize(trackindex + 1, 0.);
    }
    tracktocharge[trackindex] = c;
  }

  void setTrackMass(int trackindex, float m)
  {
    if (trackindex >= tracktomass.size()) {
      tracktomass.resize(trackindex + 1, 0.);
    }
    tracktomass[trackindex] = m;
  }

  void incStepCount(int trackindex)
  {
    if (trackindex >= stepcounterpertrack.size()) {
      stepcounterpertrack.resize(trackindex + 1, 0);
    }
    stepcounterpertrack[trackindex]++;
  }

  void setProducedSecondary(int trackindex, bool b)
  {
    if (trackindex >= producedsecondary.size()) {
      producedsecondary.resize(trackindex + 1, false);
    }
    producedsecondary[trackindex] = b;
  }

  void setCrossedBoundary(int trackindex, bool b)
  {
    if (trackindex >= crossedboundary.size()) {
      crossedboundary.resize(trackindex + 1, false);
    }
    crossedboundary[trackindex] = b;
  }

  void setTrackEnergy(int trackindex, float e)
  {
    if (trackindex >= tracktoenergy.size()) {
      tracktoenergy.resize(trackindex + 1, -1.);
    }
    // only take the starting energy
    if (tracktoenergy[trackindex] == -1.) {
      tracktoenergy[trackindex] = e;
    }
  }

  void setTrackOrigin(int trackindex, int volID)
  {
    if (trackindex >= trackorigin.size()) {
      trackorigin.resize(trackindex + 1, -1);
    }
    trackorigin[trackindex] = volID;
  }

  bool initSensitiveVolLookup(std::string const& filename);
  bool initGeoVolIdToVMCVolid(TVirtualMC* mc);

  void insertParent(int trackindex, int parent)
  {
    constexpr int PRIMARY = -1;
    if (trackindex >= tracktoparent.size()) {
      tracktoparent.resize(trackindex + 1, PRIMARY);
    }
    auto prev = tracktoparent[trackindex];
    if (prev != PRIMARY && prev != parent) {
      std::cerr << "Warning: Seeing more than one parent for same trackID\n";
    }
    tracktoparent[trackindex] = parent;
  }

  void clearTrackLookups()
  {
    tracktoparent.clear();
    producedsecondary.clear();
    crossedboundary.clear();
    tracktopdg.clear();
    stepcounterpertrack.clear();
    tracktoenergy.clear();
    tracktocharge.clear();
    tracktomass.clear();
  }

 private:
  void insertValueAt(int index, std::string const& s, std::vector<std::string*>& container)
  {
    if (index >= container.size()) {
      container.resize(index + 1, nullptr);
    }
    //#ifdef CHECKMODE
    // check that if a value exists at some index it is the same that we want to write
    if (container[index] != nullptr) {
      auto previous = *(container[index]);
      if (s.compare(previous) != 0) {
        std::cerr << "trying to override " << previous << " with " << s << "\n";
      }
    }
    //#endif
    // can we use unique pointers??
    container[index] = new std::string(s);
  }

  ClassDefNV(StepLookups, 1);
};

// struct collecting info about one MC step done
struct StepInfo {
  StepInfo() = default;
  // construct directly using virtual mc
  StepInfo(TVirtualMC* mc);

  // long cputimestamp;
  int stepid = -1; // serves as primary key
  int volId = -1;  // keep another branch somewhere mapping this to name, medium, etc.
  int copyNo = -1;
  int trackID = -1;
  float t = 0.;
  float x = 0.;
  float y = 0.;
  float z = 0.;
  float E = 0.;
  float px = 0.;
  float py = 0.;
  float pz = 0.;
  float edep = 0.;
  float mass = 0.;
  float step = 0.;
  float maxstep = 0.;
  int nsecondaries = 0;
  int prodprocess = -1;           // prod process of current track
  int nprocessesactive = 0;       // number of active processes
  bool disappeared = false;       // if track disappeared
  bool stopped = false;           // if track was stopped during last step
  bool entered = false;           // if track entered volume during last step
  bool exited = false;            // if track exited volume during last step
  bool inside = false;            // if track is inside world
  bool outside = false;           // if track is outside world
  bool newtrack = false;          // if track is new
  bool insensitiveRegion = false; // whether step done in sensitive region
  std::string* geopath = nullptr;

  const char* getProdProcessAsString() const;

  static int stepcounter;           //!
  static StepInfo* currentinstance; //!
  static std::chrono::time_point<std::chrono::high_resolution_clock> starttime;
  static void resetCounter() { stepcounter = -1; }
  static std::map<std::string, std::string>* volnametomodulemap;
  static std::vector<std::string*> volidtomodulevector;

  static StepLookups lookupstructures;
  ClassDefNV(StepInfo, 3);
};

struct MagCallInfo {
  MagCallInfo() = default;
  MagCallInfo(TVirtualMC* mc, float x, float y, float z, float Bx, float By, float Bz);

  long id = -1;
  long stepid = -1; // cross-reference to current MC stepid (if any??)
  float x = 0.;
  float y = 0.;
  float z = 0.;
  float B = 0.; // absolute value of the B field

  static int stepcounter;
  ClassDefNV(MagCallInfo, 1);
};
} // namespace o2
#endif
