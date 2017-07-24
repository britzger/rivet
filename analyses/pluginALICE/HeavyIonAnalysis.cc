// -*- C++ -*-
#include "pluginALICE/HeavyIonAnalysis.hh"

namespace Rivet {
  HeavyIonAnalysis::HeavyIonAnalysis(const std::string &name) : 
      Analysis(name),
      _centmethod(Undefined)
    { }

  void HeavyIonAnalysis::setCentralityMethod(const CentralityMethod method,
					     const size_t numEventsRequired,
					     const FinalState *fs) {
    switch (method) {

    case ImpactParameter:
      std::cout << "Method: impact parameter" << std::endl;
      _methodName = "ImpactParameter";
      _histCentralityCalibration = bookHisto1D("calib_impactpar", 200, 0, 20.0, "Calibration histogram", "xlabel", "ylabel");
      std::cout << "Calibration histogram entries: " << _histCentralityCalibration->numEntries() << std::endl;
      _histCentralityControl = bookHisto1D("control_impactpar", 200, 0, 20.0, "Control histogram", "xlabel", "ylabel");
      break;
      
    case Multiplicity:
      if (!fs)
	throw Exception("something bad happened ...");
      addProjection(*fs, "FS");
      _methodName = "Multiplicity";
      _histCentralityCalibration = bookHisto1D("calib_multiplicity_" + trimString(fs->getCuts()->description()), 10000, 0, 10000.0, "Calibration histogram", "xlabel", "ylabel");
      _histCentralityControl = bookHisto1D("control_multiplicity_" + trimString(fs->getCuts()->description()), 10000, 0, 10000.0, "Control histogram", "xlabel", "ylabel");
      break;
      
    default:
      centmethod_not_found();
    }
    
    _numEventsRequired = numEventsRequired;
    _centmethod = method;
  }

  void HeavyIonAnalysis::centmethod_not_found() const {
    throw Exception("Unimplemented method of centrality calibration.");
  }

  double HeavyIonAnalysis::centrality(const Event& event) const {
    // Initialize value of the event observable to be filled to the histograms
    if (_centmethod == Undefined)
      throw Exception("no centrality method defined");

    const float observable = calculateObservable(event);
    float centrality = -1.;

    // Check if there are enough events in the calibration histogram
    if (_histCentralityCalibration->numEntries() >= _numEventsRequired) {  
      // Calculate centrality as percentile using observable calculated with the selected method
      centrality = quantile(observable, _histCentralityCalibration) * 100.;
      if (_centmethod == Multiplicity)
	centrality = 100. - centrality;
      // Fill the control histogram with the impact parameter value and event weight
      MSG_INFO("Adding control point nr " << _histCentralityControl->numEntries() << ": (" << observable << ", " << event.weight() << ")");
      _histCentralityControl->fill(observable, event.weight());
    }
    // Otherwise, if there are not enough events in the calibration histogram
    else {
      // Fill the calibration histogram with the impact parameter value and event weight
      MSG_INFO("Adding calibration point nr " << _histCentralityCalibration->numEntries() << ": (" << observable << ", " << event.weight() << ")");
      _histCentralityCalibration->fill(observable, event.weight());
    }
    return centrality;
  }

  float HeavyIonAnalysis::calculateObservable(const Event& e) const {
    
    // Initialize observable
    float observable = -1.;
    
    // Calculate its value according to the selected method
    switch (_centmethod) {
      
    case ImpactParameter:
	// Get impact parameter of the event
	observable = e.genEvent()->heavy_ion() ? e.genEvent()->heavy_ion()->impact_parameter() : -1.;
	break;
	
    case Multiplicity:
      {
	// Get multiplicity of the event
	const FinalState& fs = applyProjection<FinalState>(e, "FS");
	observable = e.genEvent()->heavy_ion() ? fs.particles().size() : -1.;
	break;
      }
	
    default:
      centmethod_not_found();
    }
    
    // Check if observable is correct. If not, skip this event
    checkObservable(observable);
    
    return observable;
  }
  
  
  void HeavyIonAnalysis::checkObservable(const float observable) const {
    
    try {
      // Check if the observable is greater than 0
      if (observable < 0.)
	throw UserError("Calculated observable is lower than 0!");
      // Check if the observable is inside histogram ranges
      if ((observable < _histCentralityCalibration->xMin()) ||
	  (observable > _histCentralityCalibration->xMax()) ||
	  (_histCentralityCalibration->binIndexAt(observable) < 0))
	throw RangeError("Calculated observable is out of bounds!");
    } catch (...) {
      MSG_WARNING("Skipping this event...");
    }
    
  }

  string HeavyIonAnalysis::trimString(const string& str) const {
    // @note Which characters should be removed?
    string trimmedString = str;

    const char characters[] = "\\ ";
    for (const auto character : characters) {
      trimmedString.erase (std::remove (trimmedString.begin(), trimmedString.end(), character), trimmedString.end());
    }
    return trimmedString;
  }
  
  // Function for calculating the quantile from histogram at selected value
  // (below this value, including overflow)
  // @note this function is here temporarily, it will be moved and modified 
  // to create a general function for calculating the quantile of a histogram
  double HeavyIonAnalysis::quantile(double value, const Histo1DPtr hist) const {
    
    // Check if there are entries in the histogram
    if (hist->numEntries() == 0) {
      throw WeightError("There are no entires in the histogram!");
    }
    
    // Integration ranges
    size_t upperBin = hist->binIndexAt(value);

    // Calculate centrality as percentile
    return hist->integralTo(upperBin, true) / hist->integral(true);
    
  }
}
