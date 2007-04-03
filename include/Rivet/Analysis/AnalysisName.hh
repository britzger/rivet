// -*- C++ -*-
//
// This is the forward declaration of the Analysis class.
//
#ifndef RIVET_AnalysisName_HH
#define RIVET_AnalysisName_HH

#include "Rivet/Rivet.hh"

namespace Rivet {


  /// List of known available analyses
  enum AnalysisName { 
    ANALYSIS_TEST, 
    ANALYSIS_HZ95108, 
    ANALYSIS_PL273B181, 
    ANALYSIS_HEPEX0112029,
    ANALYSIS_HEPEX0409040
  };


  /// Typedef for a map of analysis name enums to strings.
  typedef std::map<AnalysisName, std::string> AnalysisMap;


  /// Typedef for a map of analysis name strings to enums.
  typedef std::map<std::string, AnalysisName> AnalysisMapR;


  /// Function which returns a map from analysis enums to the corresponding name strings.
  inline AnalysisMap getKnownAnalyses() {
    AnalysisMap amap;
    amap[ANALYSIS_TEST] = "TEST";
    amap[ANALYSIS_HZ95108] = "HZ95108";
    amap[ANALYSIS_PL273B181] = "PL273B181";
    amap[ANALYSIS_HEPEX0112029] = "HEPEX0112029";
    amap[ANALYSIS_HEPEX0409040] = "HEPEX0409040";
    return amap;
  }

  /// Function which returns a map from analysis name strings to the corresponding enums.
  inline AnalysisMapR getKnownAnalysesR() {
    AnalysisMap amap = getKnownAnalyses();
    AnalysisMapR amapr;
    for (AnalysisMap::const_iterator a = amap.begin(); a != amap.end(); ++a) {
      amapr[a->second] = a->first;
    }
    return amapr;
  }


  /// Typedef for a collection of analysis name enums.
  typedef std::set<AnalysisName> AnalysisList;
  typedef const std::set<AnalysisName> cAnalysisList;


  /// Function which returns a vector of all the analysis values in 
  /// the AnalysisName enum.
  inline AnalysisList getKnownAnalysisEnums() {
    AnalysisList names;
    AnalysisMap amap = getKnownAnalyses();
    for (AnalysisMap::const_iterator a = amap.begin(); a != amap.end(); ++a) {
      names.insert(a->first);
    }
    return names;
  }


  /// Function which returns a vector of all the analysis name strings.
  inline std::vector<std::string> getKnownAnalysisNames() {
    std::vector<std::string> names;
    AnalysisMap amap = getKnownAnalyses();
    for (AnalysisMap::const_iterator a = amap.begin(); a != amap.end(); ++a) {
      names.push_back(a->second);
    }
    return names;
  }


}


#endif
