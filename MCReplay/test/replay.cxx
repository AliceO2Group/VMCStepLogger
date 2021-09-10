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

// Test basic replay functionality

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MCReplayCoreTest

#include <boost/test/unit_test.hpp>
#include <boost/test/results_collector.hpp>

#include "MCReplay/MCReplayGenericApplication.h"
#include "MCReplay/MCReplayGenericStack.h"
#include "MCReplay/MCReplayEngine.h"

BOOST_AUTO_TEST_CASE(testReplay)
{

  // require
  // 1: Path to MCStepLogger output ROOT file
  // 2. Name of TTree were steps have been logged
  // 3. Path to ROOT geometry file
  // 4. Key name were to find geometry therein
  BOOST_REQUIRE(boost::unit_test::framework::master_test_suite().argc == 5);

  mcreplay::MCReplayGenericStack stack;
  mcreplay::MCReplayGenericApplication app{boost::unit_test::framework::master_test_suite().argv[3], boost::unit_test::framework::master_test_suite().argv[4], boost::unit_test::framework::master_test_suite().argv[1], boost::unit_test::framework::master_test_suite().argv[2]};
  mcreplay::MCReplayEngine mc{boost::unit_test::framework::master_test_suite().argv[1], boost::unit_test::framework::master_test_suite().argv[2]};

  mc.SetStack(&stack);
  app.setStack(&stack);
  mc.Init();

  // replay all events
  BOOST_REQUIRE(mc.ProcessRun(-1));
}
