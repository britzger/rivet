// -*- C++ -*-
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Analysis.hh"
#include "YODA/WriterYODA.h"
#include "YODA/ReaderYODA.h"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Centrality.hh"

#include <algorithm>

namespace Rivet {
  
  Centrality::Centrality(const FinalState &fs, const Method method, const size_t numEventsRequired) :
    _cutString("true"),
    _method(method),
    _numEventsRequired(numEventsRequired),
    _centrality(-1.)
  {
    addProjection(fs, "FS");
    
    switch (_method) {

    case ImpactParameter:
      _methodName = "ImpactParameter";
      _histCentralityCalibration = make_shared<Histo1D>(200, 0.0, 20.0);
      _histCentralityControl = make_shared<Histo1D>(200, 0.0, 20.0);
      break;
      
    case Multiplicity:
      _methodName = "Multiplicity";
      _histCentralityCalibration = make_shared<Histo1D>(10000, 0.0, 10000.0);
      _histCentralityControl = make_shared<Histo1D>(10000, 0.0, 10000.0);
      break;
      
    default:
      method_not_found();
      
    }

    addCut(fs.getCuts());
    MSG_INFO("Using centrality projection for cuts " + _cutString + ", calibration method " + _methodName + ", and " << _numEventsRequired << " required events");
  }
  
  
  int Centrality::compare(const Projection& p) const {
    // @note Is it okay to cast projetion to centrality? What if it fails?
    const Centrality& other = dynamic_cast<const Centrality&>(p);
    const FinalState& fs = getProjection<FinalState>("FS");
    const FinalState& otherfs = other.getProjection<FinalState>("FS");
    return 					\
      (cmp(fs.getCuts(), otherfs.getCuts()) || 
       cmp(_method, other._method) || 
       cmp(_numEventsRequired, other._numEventsRequired) || 
       mkNamedPCmp(p, "FS"));
  }
  

  void Centrality::project(const Event& e) {
    
    // Set the default value of centrality
    _centrality = -1.;

    // Initialize value of the event observable to be filled to the histograms
    const float observable = calculateObservable(e);
    
    // Check if there are enough events in the calibration histogram
    if (_histCentralityCalibration->numEntries() >= _numEventsRequired) {  
      // Calculate centrality as percentile using observable calculated with the selected method
      calculateCentralityPercentile(observable);
      // Fill the control histogram with the impact parameter value and event weight
      MSG_DEBUG("Adding control point nr " << _histCentralityControl->numEntries() << ": (" << observable << ", " << e.weight() << ")");
      _histCentralityControl->fill(observable, e.weight());
    }
    // Otherwise, if there are not enough events in the calibration histogram
    else {
      // Fill the calibration histogram with the impact parameter value and event weight
      MSG_DEBUG("Adding calibration point nr " << _histCentralityCalibration->numEntries() << ": (" << observable << ", " << e.weight() << ")");
      _histCentralityCalibration->fill(observable, e.weight());
    }
    
  }

  
  float Centrality::calculateObservable(const Event& e) const {
    
    // Initialize observable
    float observable = -1.;
    
    // Calculate its value according to the selected method
    switch (_method) {
      
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
      method_not_found();
    }
    
    // Check if observable is correct. If not, skip this event
    checkObservable(observable);
    
    return observable;
  }
  
  
  void Centrality::checkObservable(const float observable) const {
    
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
  
  
  void Centrality::calculateCentralityPercentile(const float observable) {
    
    // Default integration ranges
    size_t bin1 = 0;
    size_t bin2 = _histCentralityCalibration->numBins() - 1;
    
    // Select proper integration ranges according to the selected method
    switch (_method) {
    
    case ImpactParameter:
	// Calculate centrality using impact parameter distribution
	bin2 = _histCentralityCalibration->binIndexAt(observable);
	break;
    
    case Multiplicity:
	// Calculate centrality using multiplicity distribution
	bin1 = _histCentralityCalibration->binIndexAt(observable);
	break;
      
    default:
      method_not_found();
    }
    
    // Calculate centrality as percentile
    if ((0 <= bin1) && (bin1 <= bin2) && (bin2 < _histCentralityCalibration->numBins())) {
      _centrality = _histCentralityCalibration->integralRange(bin1, bin2) * 100.0 / _histCentralityCalibration->integral();
    }
    
  }
  
  
  void Centrality::saveHistogram(const string& filename, const Histo1DPtr hist) const {
    Histo1DPtr histCopy = make_shared<Histo1D>(*hist);
    histCopy->scaleW(1. / hist->numEntries());
    try {
      MSG_INFO("Saving the centrality histogram to " << filename << "...");
      YODA::WriterYODA::write(filename, *histCopy);
    } catch (...) {
      throw UserError("Unexpected error in writing file to: " + filename);
    }
    histCopy.reset();
  }
  
  
  void Centrality::finalize() const {
    
    // Save histograms
    saveHistogram(name() + "_calibration.yoda", _histCentralityCalibration);
    if (_histCentralityControl->numEntries() > 0)
      saveHistogram(name() + "_control.yoda", _histCentralityControl);
    
    // @note Is there anything else that should be done here?
    
  }
  
  
  void Centrality::addCalibrationFile(const string& filename) {
    
    // Open the calibration file and read all objects
    std::vector<AnalysisObject*> aos;
    YODA::Reader& reader = YODA::mkReader(filename);
    std::ifstream istring(filename);
    try {
      MSG_INFO("Reading the calibration file: " << filename << "...");
      reader.read(istring, aos);
    } catch (const YODA::ReadError& err) {
      throw UserError("Unexpected error in reading the file " + filename + ": " + err.what());
    }
    
    // Create a histogram path template to look for
    const string histPath = "/Centrality/" + _methodName + "/" + _cutString + "/calibration";
    MSG_INFO("Searching for: " << histPath);
    
    // Extract the proper calibration histogram from the file
    for (const auto ao : aos) {
      AnalysisObjectPtr aoPtr(ao);
      Histo1DPtr hist = dynamic_pointer_cast<Histo1D>(aoPtr);
      if (hist != nullptr) {
	MSG_DEBUG(hist->path() << " ? " << histPath);
	if (hist->path().compare(histPath) == 0) {
	  // Check if calibration histogram already exists
	  if (_histCentralityCalibration && (_histCentralityCalibration->numEntries() > 0)) {
	    MSG_WARNING("Calibration histogram already exists! The histogram found in " << filename << " file will be omitted!");
	    return;
	  }
	  MSG_INFO("Histogram found for centrality projection. Setting the new calibration histogram...");
	  _histCentralityCalibration = make_shared<Histo1D>(*hist);
	}
      }
      else {
	MSG_ERROR("Unsupported type of histogram");
      }
    }
    
  }
  
  void Centrality::trimString(string& str) const {
    // @note Which characters should be removed?
    const char characters[] = "\\ ";
    for (const auto character : characters) {
      str.erase (std::remove (str.begin(), str.end(), character), str.end());
    }
  }
  
  void Centrality::addCut(const Cut& cut) {
    addCutString(cut->description());
  }
  
  void Centrality::addCutString(string cutString) {
    trimString(cutString);
    _cutString = cutString;
    _histCentralityCalibration->setPath("/Centrality/" + _methodName + "/" + _cutString + "/calibration");
    _histCentralityCalibration->setTitle("/Centrality/" + _methodName + "/" + _cutString + "/calibration");
    _histCentralityControl->setPath("/Centrality/" + _methodName + "/" + _cutString + "/control");
    _histCentralityControl->setTitle("/Centrality/" + _methodName + "/" + _cutString + "/control");
    setName("Centrality_" + _methodName + "_" + _cutString);
  }
  
  void Centrality::method_not_found() const {
    throw Exception("Unimplemented method of centrality calibration.");
  }
  
}
