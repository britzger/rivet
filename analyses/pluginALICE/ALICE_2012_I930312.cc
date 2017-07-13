// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

#define PT_BINS 4
#define EVENT_TYPES 3

namespace Rivet {
  
  double quantile(double value, const Histo1DPtr hist);
  
  class ALICE_2012_I930312 : public Analysis {
  
  public:
    
    ALICE_2012_I930312() : 
      Analysis("ALICE_2012_I930312")
    {  }


    void init() {

      // Charged final states with |eta| < 1.0 and 8 < pT < 15 GeV/c for trigger particles
      const Cut& cutTrigger = Cuts::abseta < 1.0 && Cuts::pT > 8*GeV && Cuts::pT < 15*GeV;
      const ChargedFinalState cfsTrigger(cutTrigger);
      addProjection(cfsTrigger,"CFSTrigger");
      
      // Chaged final states with |eta| < 1.0 and 3 < pT < 15 GeV/c for associated particles
      for (int ipt = 0; ipt < PT_BINS; ipt++)
	num_trig[ipt] = 0;
      
      pt_limits[0] = 3;
      pt_limits[1] = 4;
      pt_limits[2] = 6;
      pt_limits[3] = 8;
      pt_limits[4] = 10;
      
      // Cut cutAssoc[4];
      // ChargedFinalState cfsAssoc[4];
      for (int ipt = 0; ipt < PT_BINS; ipt++) {
	
	Cut mycut = Cuts::abseta < 1.0 && Cuts::pT > pt_limits[ipt]*GeV && Cuts::pT < pt_limits[ipt + 1]*GeV;
	addProjection(ChargedFinalState(mycut), "CFSAssoc" + std::to_string(ipt));

	//ChargedFinalState mycfs = ChargedFinalState(mycut);
	// const ChargedFinalState &mycfs(mycut);
	//addProjection(mycfs, "CFSAssoc" + std::to_string(ipt));
      }
      
      // Impact parameter distribution histogram for centrality calculation. If the histogram with the same name
      // is loaded in readData() function, it will bind such preloaded histogram to the following histogram.
      // Otherwise, empty histogram will be created
      // @note the name was not discussed yet, this is only an example: calib_mult...
      std::cout << "Initialization..." << std::endl;
      _histCalibration = bookHisto1D("calib_impactpar_(eta<0.5)&&(pT>50*MeV)", 200, 0, 20.0, "Calibration histogram", "xlabel", "ylabel");
      _histControl = bookHisto1D("control_impactpar_(eta<0.5)&&(pT>50*MeV)", 200, 0, 20.0, "Control histogram", "xlabel", "ylabel");
      
      // Histograms and variables initialization. Do it for each centrality range
      for (int ipt = 0; ipt < PT_BINS; ipt++) {
	_histYield[0][ipt] = bookHisto1D("Yield_pp_" + std::to_string(ipt), 36, -0.5 * M_PI, 1.5 * M_PI,
							   "Associated particle per trigger particle yield", "#Delta#eta (rad)", "1 / N_trig dN_assoc / d#Delta#eta (rad^-1)");
	_histYield[1][ipt] = bookHisto1D("Yield_central_" + std::to_string(ipt), 36, -0.5 * M_PI, 1.5 * M_PI,
							   "Associated particle per trigger particle yield", "#Delta#eta (rad)", "1 / N_trig dN_assoc / d#Delta#eta (rad^-1)");
	_histYield[2][ipt] = bookHisto1D("Yield_peripheral_" + std::to_string(ipt), 36, -0.5 * M_PI, 1.5 * M_PI,
							   "Associated particle per trigger particle yield", "#Delta#eta (rad)", "1 / N_trig dN_assoc / d#Delta#eta (rad^-1)");
	_histYieldBkgRemoved[0][ipt] = bookHisto1D("Yield_pp_nobkg_" + std::to_string(ipt), 36, -0.5 * M_PI, 1.5 * M_PI,
							   "Associated particle per trigger particle yield no bkg", "#Delta#eta (rad)", "1 / N_trig dN_assoc / d#Delta#eta (rad^-1)");
	_histYieldBkgRemoved[1][ipt] = bookHisto1D("Yield_central_nobkg_" + std::to_string(ipt), 36, -0.5 * M_PI, 1.5 * M_PI,
							   "Associated particle per trigger particle yield no bkg", "#Delta#eta (rad)", "1 / N_trig dN_assoc / d#Delta#eta (rad^-1)");
	_histYieldBkgRemoved[2][ipt] = bookHisto1D("Yield_peripheral_nobkg_" + std::to_string(ipt), 36, -0.5 * M_PI, 1.5 * M_PI,
							   "Associated particle per trigger particle yield no bkg", "#Delta#eta (rad)", "1 / N_trig dN_assoc / d#Delta#eta (rad^-1)");
	
      }
      _histTriggerCounter[0] = bookHisto1D("Trigger_pp", PT_BINS, 0.0, PT_BINS, "Trigger counter pp", "pt bin", "N");
      _histTriggerCounter[1] = bookHisto1D("Trigger_central", PT_BINS, 0.0, PT_BINS, "Trigger counter central", "pt bin", "N");
      _histTriggerCounter[2] = bookHisto1D("Trigger_peripheral", PT_BINS, 0.0, PT_BINS, "Trigger counter peripheral", "pt bin", "N");
      
      _histIAA[0] = bookHisto1D(1, 1, 1);
      _histIAA[1] = bookHisto1D(4, 1, 1);
      _histIAA[2] = bookHisto1D(5, 1, 1);
      _histIAA[3] = bookHisto1D(6, 1, 1);
      
    }

    void analyze(const Event& event) {
      
      if (event.genEvent()->heavy_ion()) {
	std::cout << "HI EVENT" << std::endl;
      } else {
	std::cout << "PP EVENT" << std::endl;
      }
      std::cout << event.genEvent()->heavy_ion() << std::endl;
      
      const ChargedFinalState& triggerFinalState = applyProjection<ChargedFinalState>(event, "CFSTrigger");
      Particles triggerParticles = triggerFinalState.particlesByPt();
      std::cout << "Trigger particles: " << triggerParticles.size() << std::endl;
      
      ChargedFinalState associatedFinalState[PT_BINS];
      Particles associatedParticles[PT_BINS];
      std::cout << "Associated particles: ";
      for (int ipt = 0; ipt < PT_BINS; ipt++) {
	associatedFinalState[ipt] = applyProjection<ChargedFinalState>(event, "CFSAssoc" + std::to_string(ipt));
	associatedParticles[ipt] = associatedFinalState[ipt].particlesByPt();
	std::cout << associatedParticles[ipt].size() << " ";
      }
      std::cout << std::endl;
      
      // Check event type
      if (event.genEvent()->heavy_ion()) {
	// For PbPb events fill calibration and control histogram
	if (_histCalibration->numEntries() >= 10000) {
	  _histControl->fill(event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->impact_parameter() : -1.0, 1);
	} else {
	  _histCalibration->fill(event.genEvent()->heavy_ion() ? event.genEvent()->heavy_ion()->impact_parameter() : -1.0, 1);
	  // Veto event in case there are not enough events in calibration histogram
	  vetoEvent;
	}
	double centr = quantile(event.genEvent()->heavy_ion()->impact_parameter(), _histCalibration) * 100.0;
	if (centr > 0.0 && centr < 5.0)
	  event_type = 1; // PbPb, central
	else if (centr > 60.0 && centr < 90.0)
	  event_type = 2; // PbPb, peripherial
	else
	  event_type = 3; // PbPb, other
      } else {
	event_type = 0; // pp
      }
      std::cout << "Event type: " << event_type << std::endl;
      
      if (event_type == 3)
	vetoEvent;
      
      // Loop over trigger particles
      foreach (const Particle& triggerParticle, triggerParticles) {
	
	for (int ipt = 0; ipt < PT_BINS; ipt++) {	  
	  if (associatedParticles[ipt].size() > 1 || 
	      (associatedParticles[ipt].size() == 1 && associatedParticles[ipt].at(0) != triggerParticle)) {
	    num_trig[ipt] += 1;
	    _histTriggerCounter[event_type]->fill(ipt, 1);
	    std::cout << "Found some associated particles" << std::endl;
	    foreach (const Particle& associatedParticle, associatedParticles[ipt]) {
	      if (triggerParticle.pt() > associatedParticle.pt()) {
		double dPhi = triggerParticle.phi() - associatedParticle.phi();
		while (dPhi > 1.5 * M_PI)  { dPhi -= 2 * M_PI; }
		while (dPhi < -0.5 * M_PI) { dPhi += 2 * M_PI; }
		std::cout << "Filling yield histogram for pt = " << pt_limits[ipt] << "-" << pt_limits[ipt + 1] << " GeV/c for delta phi = " << dPhi << std::endl;
		_histYield[event_type][ipt]->fill(dPhi, 1);
	      }
	    }
	  }
	  else {
	    std::cout << "No associated particles with pt = " << pt_limits[ipt] << "-" << pt_limits[ipt + 1] << " GeV/c found for that trigger particle" << std::endl;
	  }
	}
	
      }
      
    }

    void finalize() { 
      
      std::cout << "Finalising..." << std::endl;
      // Normalize all histograms by number of trigger particles used to create each one of them
      /*for (int ipt = 0; ipt < PT_BINS; ipt++) {
	std::cout << "finalize(): number of trigger particles: " << num_trig[ipt] << std::endl;
	if (num_trig[ipt] != 0) {
	  scale(_histYield[event_type][ipt], (1. / num_trig[ipt]));
	}
	}*/
      
    }
    
    void post() {
      
      double nearSide[EVENT_TYPES][PT_BINS] = { {0.0} };
      double awaySide[EVENT_TYPES][PT_BINS] = { {0.0} };
      
      std::cout << "Calibration entries: " << _histCalibration->numEntries() << std::endl;
      std::cout << "Control entries: " << _histControl->numEntries() << std::endl;
      std::cout << "pp trigger counter histogram entries: " <<_histTriggerCounter[0]->numEntries() << std::endl;
      std::cout << "Central trigger counter histogram entries: " << _histTriggerCounter[1]->numEntries() << std::endl;
      std::cout << "Peripheral trigger counter histogram entries: " << _histTriggerCounter[2]->numEntries() << std::endl;
      std::cout << "pp yield pt bin 0, entries: " << _histYield[0][0]->numEntries() << std::endl;
      std::cout << "pp yield pt bin 1, entries: " << _histYield[0][1]->numEntries() << std::endl;
      std::cout << "pp yield pt bin 2, entries: " << _histYield[0][2]->numEntries() << std::endl;
      std::cout << "pp yield pt bin 3, entries: " << _histYield[0][3]->numEntries() << std::endl;
      std::cout << "Central yield pt bin 0, entries: " << _histYield[1][0]->numEntries() << std::endl;
      std::cout << "Central yield pt bin 1, entries: " << _histYield[1][1]->numEntries() << std::endl;
      std::cout << "Central yield pt bin 2, entries: " << _histYield[1][2]->numEntries() << std::endl;
      std::cout << "Central yield pt bin 3, entries: " << _histYield[1][3]->numEntries() << std::endl;
      std::cout << "Peripheral yield pt bin 0, entries: " << _histYield[2][0]->numEntries() << std::endl;
      std::cout << "Peripheral yield pt bin 1, entries: " << _histYield[2][1]->numEntries() << std::endl;
      std::cout << "Peripheral yield pt bin 2, entries: " << _histYield[2][2]->numEntries() << std::endl;
      std::cout << "Peripheral yield pt bin 3, entries: " << _histYield[2][3]->numEntries() << std::endl;
      
      // For each event type
      for (int itype = 0; itype < EVENT_TYPES; itype++) {
	// For each pT range
	for (int ipt = 0; ipt < PT_BINS; ipt++) {
	  
	  // Check if histograms are fine
	  if (_histTriggerCounter[itype]->numEntries() == 0 || _histYield[itype][ipt]->numEntries() == 0) {
	    std::cout << _histTriggerCounter[itype]->numEntries() << std::endl;
	    std::cout << _histYield[itype][ipt]->numEntries() << std::endl;
	    std::cout << "There are no entries in one of the histograms" << std::endl;
	    continue;
	  }
	  
	  // Scale yield histogram
	  std::cout << "Scaling yield histograms..." << std::endl;
	  if (_histTriggerCounter[itype]->bin(ipt).sumW() != 0) {
	    std::cout << "Scaling histogram type " << itype << " pt bin " << ipt << "..." << std::endl;
	    scale(_histYield[itype][ipt], (1. / _histTriggerCounter[itype]->bin(ipt).sumW()));
	  }
	  
	  // Calculate background
	  std::cout << "Calculating background" << std::endl;
	  double sum = 0.0;
	  int nbins = 0;
	  for (unsigned int ibin = 0; ibin < _histYield[itype][ipt]->numBins(); ibin++) {
	    std::cout << "Bin " << ibin << std::endl;
	    if (_histYield[itype][ipt]->bin(ibin).xMid() > (0.5 * M_PI - 0.4) &&
		_histYield[itype][ipt]->bin(ibin).xMid() < (0.5 * M_PI + 0.4)) {
	      std::cout << "Adding " << _histYield[itype][ipt]->bin(ibin).sumW() << std::endl;
	      sum += _histYield[itype][ipt]->bin(ibin).sumW();
	      nbins++;
	    }
	  }
	  if (nbins == 0) {
	    std::cout << "Failed to estimate background!" << std::endl;
	    continue;
	  }
	  std::cout << "Dividing " << sum << " / " << nbins << std::endl;
	  double bkg = sum / nbins;
	  
	  // Fill histograms with removed background
	  std::cout << "Filling histogram with removed background..." << std::endl;
	  for (unsigned int ibin = 0; ibin < _histYield[itype][ipt]->numBins(); ibin++) {
	    std::cout << "Bin " << ibin << ", value " << _histYield[itype][ipt]->bin(ibin).sumW() - bkg << std::endl;
	    _histYieldBkgRemoved[itype][ipt]->fillBin(ibin, _histYield[itype][ipt]->bin(ibin).sumW() - bkg);
	  }
	  
	  std::cout << "Integrating the whole histogram: " << _histYieldBkgRemoved[itype][ipt]->integral() << std::endl;
	  
	  // Integrate near-side yield
	  std::cout << "Integrating near-side yield..." << std::endl;
	  unsigned int lowerBin = _histYieldBkgRemoved[itype][ipt]->binIndexAt(-0.7);
	  unsigned int upperBin = _histYieldBkgRemoved[itype][ipt]->binIndexAt(0.7);
	  nearSide[itype][ipt] = _histYieldBkgRemoved[itype][ipt]->integralRange(lowerBin, upperBin);
	  std::cout << "Integrating bins (" << lowerBin << ", " << upperBin << "): " << nearSide[itype][ipt] << std::endl;
	  // Integrate away-side yield
	  std::cout << "Integrating away-side yield..." << std::endl;
	  lowerBin = _histYieldBkgRemoved[itype][ipt]->binIndexAt(M_PI - 0.7);
	  upperBin = _histYieldBkgRemoved[itype][ipt]->binIndexAt(M_PI + 0.7);
	  awaySide[itype][ipt] = _histYieldBkgRemoved[itype][ipt]->integralRange(lowerBin, upperBin);
	  std::cout << "Integrating bins (" << lowerBin << ", " << upperBin << "): " << awaySide[itype][ipt] << std::endl;
	  
	}
      }
      
      // Print I_CP results
      std::cout << "Near side" << std::endl;
      for (int ipt = 0; ipt < PT_BINS; ipt++) {
	std::cout << "I_AA, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << nearSide[0][ipt] / nearSide[1][ipt] << std::endl;
	_histIAA[0]->fillBin(ipt, nearSide[0][ipt] / nearSide[1][ipt]);
	std::cout << "I_CP, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << nearSide[1][ipt] / nearSide[2][ipt] << std::endl;
	_histIAA[1]->fillBin(ipt, nearSide[1][ipt] / nearSide[2][ipt]);
      }
      std::cout << "Away side" << std::endl;
      for (int ipt = 0; ipt < PT_BINS; ipt++) {	
	std::cout << "I_AA, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << awaySide[0][ipt] / awaySide[1][ipt] << std::endl;
	_histIAA[2]->fillBin(ipt, awaySide[0][ipt] / awaySide[1][ipt]);
	std::cout << "I_CP, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << awaySide[1][ipt] / awaySide[2][ipt] << std::endl;
	_histIAA[3]->fillBin(ipt, awaySide[1][ipt] / awaySide[2][ipt]);
      }
      
    }
    
  private:
    
    Histo1DPtr _histCalibration;
    Histo1DPtr _histControl;
    Histo1DPtr _histYield[EVENT_TYPES][PT_BINS];
    Histo1DPtr _histYieldBkgRemoved[EVENT_TYPES][PT_BINS];
    Histo1DPtr _histTriggerCounter[EVENT_TYPES];
    Histo1DPtr _histIAA[4];
    int num_trig[PT_BINS];
    int pt_limits[5];
    int event_type;
    
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I930312);
  
  // Function for calculating the quantile from histogram at selected value
  // (below this value, including overflow)
  // @note this function is here temporarily, it will be moved and modified 
  // to create a general function for calculating the quantile of a histogram
  double quantile(double value, const Histo1DPtr hist) {
    
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
