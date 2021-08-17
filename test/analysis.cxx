// Test analysis functionality

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MCStepLoggerAnalysisCoreTest

#include <fstream>

#include <boost/test/unit_test.hpp>
#include <boost/test/results_collector.hpp>

#include "MCStepLogger/MCAnalysisManager.h"
#include "MCStepLogger/MCAnalysisFileWrapper.h"
#include "MCStepLogger/SimpleStepAnalysis.h"
#include "MCStepLogger/MCAnalysisFileWrapper.h"

using namespace o2::mcstepanalysis;

// checking meta info and whether there is a value container the analysis can be conpared with
BOOST_AUTO_TEST_CASE(testAnalysis)
{

  // require input and output path
  BOOST_REQUIRE(boost::unit_test::framework::master_test_suite().argc == 3);

  new SimpleStepAnalysis();
  auto& anamgr = MCAnalysisManager::Instance();
  anamgr.setLabel("testLabel");
  anamgr.setInputFilepath(boost::unit_test::framework::master_test_suite().argv[1]);

  // everything loaded and ready to run?
  BOOST_REQUIRE(anamgr.checkReadiness());

  // runs all events found and should produce an output file
  anamgr.run();

  // write output
  std::string outputDir(boost::unit_test::framework::master_test_suite().argv[2]);
  std::string analysisFilepath = outputDir + "/SimpleStepAnalysis/Analysis.root";
  anamgr.write(outputDir);

  // Now check if the analysis file was produced properly
  std::ifstream file(analysisFilepath);
  BOOST_REQUIRE(file);

  MCAnalysisFileWrapper anaFile;
  anaFile.read(analysisFilepath);
  BOOST_REQUIRE(anaFile.isSane());

  // one last check for meta information
  auto& anaMetaInfo = anaFile.getAnalysisMetaInfo();
  BOOST_REQUIRE(anaMetaInfo.analysisName.compare("SimpleStepAnalysis") == 0);
}
