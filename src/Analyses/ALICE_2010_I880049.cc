// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet {
  
  double quantile(double value, const Histo1DPtr hist);
  
  class ALICE_2010_I880049 : public Analysis {
  
  public:
    
    ALICE_2010_I880049() : 
      Analysis("ALICE_2010_I880049")
    {  }


    void init() {

      // Charged final states with |eta| < 0.5 and pT > 50 MeV
      const Cut& cut = Cuts::abseta < 0.5 && Cuts::pT > 50*MeV;
      const ChargedFinalState cfs(cut);
      addProjection(cfs,"CFS");
      
      // Multiplicity distribution histogram for centrality calculation. If the histogram with the same name
      // is loaded in preloadData() function (which will be called before init(), to be written by Andy)
      // it will bind such preloaded histogram to the following histogram. Otherwise, empty histogram 
      // will be created
      // @note the name was not discussed yet, this is only an example: calib_mult...
      _histCalibration = bookHisto1D("calib_mult_(eta<0.5)&&(pT>50*MeV)", 2000, 0, 2000, "Calibration histogram", "xlabel", "ylabel");
      _histControl = bookHisto1D("control_mult_(eta<0.5)&&(pT>50*MeV)", 2000, 0, 2000, "Control histogram", "xlabel", "ylabel");
      
      // Histograms and variables initialization. Do it for each centrality range
      _histNchVsCentr = bookProfile1D(1, 1, 1);
      
    }

    void analyze(const Event& event) {

      const ChargedFinalState& charged = applyProjection<ChargedFinalState>(event, "CFS");
      Particles chargedParticles = charged.particlesByPt();
      
      std::cout << "num entries: " << _histCalibration->numEntries() << std::endl;
      //MSG_INFO("num entries: " + _histCalibration->numEntries());
      // Filling calibration, control, and analysis histograms
      if (_histCalibration->numEntries() >= 50) {
	_histControl->fill(charged.particles().size(), 1);
	MSG_INFO("Filling the control histogram");
	double c = quantile(charged.particles().size(), _histCalibration) * 100.0;
	_histNchVsCentr->fill(c, charged.particles().size(), event.weight());
	MSG_INFO("Filling the analysis plots");
      } else {
	_histCalibration->fill(charged.particles().size(), 1);
	MSG_INFO("Filling the calibration histogram");
      }
      
    }

    void finalize() { }
    
  private:
    
    Profile1DPtr _histNchVsCentr;
    Histo1DPtr _histCalibration;
    Histo1DPtr _histControl;
    
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2010_I880049);
  
  // Function for calculating the quantile from histogram at selected value
  // (below this value, including overflow)
  // @note this function is here temporarily, it will be moved and modified 
  // to create a general function for calculating the quantile of a histogram
  double quantile(double value, const Histo1DPtr hist) {
    
    // Check if there are entries in the histogram
    try {
      if (hist->numEntries() == 0) {
	throw WeightError("There are no entires in the histogram!");
	//throw YODA::LowStatsError("Insufficient statistics in the histogram!");
      }
    } catch (...) {
      //MSG_ERROR("Unable to calculate quantile of a histogram");
      
    }
    
    // Integration ranges
    size_t upperBin = hist->binIndexAt(value);

    // Calculate centrality as percentile
    return hist->integralTo(upperBin, true) / hist->integral(true);
    
  }

}
