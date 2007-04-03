// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Analysis abstract base class.
//

#include "Rivet/RivetHandler.hh"
#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Analysis/TestAnalysis.hh"
#include "Rivet/Analysis/HZ95108.hh"
#include "Rivet/Analysis/PL273B181.hh"
#include "Rivet/Analysis/HepEx0112029.hh"
#include "Rivet/Analysis/HepEx0409040.hh"
#include "Rivet/Tools/Logging.hh"
using namespace Rivet;

#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"
using namespace AIDA;

#include <stdexcept>


//////////////////////////////////////////////////////////////


/// If you write a new analysis, add the ID for it here.
Analysis& Analysis::getAnalysis(const AnalysisName atype) {
  switch (atype) {
  case ANALYSIS_TEST:
    return *(new TestAnalysis());
  case ANALYSIS_HZ95108:
    return *(new HZ95108());
  case ANALYSIS_PL273B181:
    return *(new PL273B181());
  case ANALYSIS_HEPEX0112029:
    return *(new HepEx0112029());
  case ANALYSIS_HEPEX0409040:
    return *(new HepEx0409040());  
  }
  throw runtime_error("Tried to get an analysis not known in the Rivet::AnalysisName enum.");
}


//////////////////////////////////////////////////////////////


Analysis::~Analysis() {}

void Analysis::init() {}


void Analysis::analyze(const Event &) {}


void Analysis::finalize() {}


RivetInfo Analysis::getInfo() const {
  return info;
}

//////////////////////////////////////////////////////////////


IAnalysisFactory& Analysis::analysisFactory() {
  return handler().analysisFactory();
}


ITree& Analysis::tree() {
  return handler().tree();
}


IHistogramFactory& Analysis::histogramFactory() {
  return handler().histogramFactory();
}


Log& Analysis::getLog() {
  string logname = "Rivet.Analysis." + name();
  return Log::getLog(logname);
}


IHistogram1D* bookHistogram1D(const unsigned int paperId, const unsigned int datasetId, 
                              const unsigned int axisId, const std::string& title) {
  throw runtime_error("Analysis::bookHistogram1D(int paperId, int datasetId, int axisId, string title) is not yet implemented.");
}

IHistogram1D* bookHistogram1D(const string hdcode, const std::string& title) {
  throw runtime_error("Analysis::bookHistogram1D(string hdcode, string title) is not yet implemented.");
}


IHistogram1D* Analysis::bookHistogram1D(const std::string& name, const std::string& title, 
                                        const int nbins, const double lower, const double upper) {
  makeHistoDir();
  return histogramFactory().createHistogram1D(histoDir() + "/" + name, title, nbins, lower, upper);
}


void Analysis::makeHistoDir() {
  if (!_madeHistoDir) {
    if (!name().empty()) {
      tree().mkdir(histoDir());
    }
    _madeHistoDir = true;
  }
}
