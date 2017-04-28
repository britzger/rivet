// -*- C++ -*-
#ifndef RIVET_Centrality_HH
#define RIVET_Centrality_HH

#include "Rivet/Projection.hh"
#include "YODA/WriterYODA.h"

namespace Rivet {


  class Centrality : public Projection {
  public:

    enum Method {
      ImpactParameter,
      Multiplicity
    };

    Centrality (const FinalState &fs, const Method method, const size_t numEventsRequired = 10000);
    
    /// Clone on the heap
    DEFAULT_RIVET_PROJ_CLONE(Centrality);

    /// Is the centrality valid?
    bool isValid() const { return (_centrality >= 0.); }
    
    /// Add calibration file
    void addCalibrationFile(const string& filename);
    
    /// Get centrality of the event
    double getCentrality() const { return _centrality; }
    
    /// Finalize the centrality calculation; save histograms
    void finalize() const;

    /// Add cut to the projection
    void addCut(const Cut& cut);
    
  protected:
    
    /// Apply the projection to the event
    virtual void project(const Event& e);
    
    /// Compare projections
    virtual int compare(const Projection& p) const;
    
    /// Calculate the value of the observable using the selected method
    float calculateObservable(const Event& e) const;
    
    /// Check if observable is correct
    void checkObservable(const float observable) const;
    
    /// Calculate centrality percentile
    void calculateCentralityPercentile(const float observable);
    
    /// Save 1D histogram
    void saveHistogram(const string& filename, const Histo1DPtr hist) const;
    
    /// Trim string by removing some special characters
    void trimString(string& str) const;
    
    /// Add cut string to the projection
    void addCutString(string cutString);
    
    /// When there is no implementation for requested method
    void method_not_found() const;
    
  private:
    
    /// Histogram for centrality calibration
    Histo1DPtr _histCentralityCalibration;
    /// Histogram for centrality control. It may be used to compare distribution 
    /// in the current run to the one provided in calibration histogram.
    Histo1DPtr _histCentralityControl;
    
    /// String with the cuts
    string _cutString;
    
    /// Method of the centrality calibration
    Method _method;
    
    /// Number of events required for a selected method of centrality calibration
    size_t _numEventsRequired;
    
    /// Centrality value
    float _centrality;
  
    /// Name of the centrality calibration method
    string _methodName;
    
  };

}


#endif
