// -*- C++ -*-
#include "Rivet/Analysis.hh"

namespace Rivet {
  
  class HeavyIonAnalysis : public Analysis {
    
  public:
    HeavyIonAnalysis(const std::string &name = "");
    
    // methods specific to heavy-ion analyses
    enum CentralityMethod {
      ImpactParameter,
      Multiplicity,
      Undefined
    };

    double quantile(double value, const Histo1DPtr hist) const;

    /// centrality implementation
    void setCentralityMethod(const CentralityMethod method,
			     const size_t numEventsRequired,
			     const FinalState *fs = 0x0);

    double centrality(const Event& event) const;

    void centmethod_not_found() const;

    float calculateObservable(const Event &event) const;

    void checkObservable(const float observable) const;

  private:
    string trimString(const string& str) const;

    /// Histogram for centrality calibration
    Histo1DPtr _histCentralityCalibration;
    /// Histogram for centrality control. It may be used to compare distribution 
    /// in the current run to the one provided in calibration histogram.
    Histo1DPtr _histCentralityControl;

    /// String with the cuts
    string _cutString;
    
    /// Method of the centrality calibration
    CentralityMethod _centmethod;
    
    /// Number of events required for a selected method of centrality calibration
    size_t _numEventsRequired;
    
    /// Name of the centrality calibration method
    string _methodName;
  };
}
