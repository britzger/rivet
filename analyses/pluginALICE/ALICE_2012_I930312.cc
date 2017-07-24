// -*- C++ -*-
#include "pluginALICE/HeavyIonAnalysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet {
  
  class ALICE_2012_I930312 : public HeavyIonAnalysis {
  
  public:
    
    ALICE_2012_I930312() : 
      HeavyIonAnalysis("ALICE_2012_I930312")
    {  }


    void init() {
      HeavyIonAnalysis::init();

      setCentralityMethod(HeavyIonAnalysis::ImpactParameter, 10000);

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
      _histTriggerCounter = bookHisto1D("Trigger", EVENT_TYPES, 0.0, EVENT_TYPES, "Trigger counter", "event type", "N");
      //_histTriggerCounter[0] = bookHisto1D("Trigger_pp", PT_BINS, 0.0, PT_BINS, "Trigger counter pp", "pt bin", "N");
      //_histTriggerCounter[1] = bookHisto1D("Trigger_central", PT_BINS, 0.0, PT_BINS, "Trigger counter central", "pt bin", "N");
      //_histTriggerCounter[2] = bookHisto1D("Trigger_peripheral", PT_BINS, 0.0, PT_BINS, "Trigger counter peripheral", "pt bin", "N");
      
      _histIAA[0] = bookHisto1D(1, 1, 1);
      _histIAA[1] = bookHisto1D(2, 1, 1);
      _histIAA[2] = bookHisto1D(5, 1, 1);

      _histIAA[3] = bookHisto1D(3, 1, 1);
      _histIAA[4] = bookHisto1D(4, 1, 1);
      _histIAA[5] = bookHisto1D(6, 1, 1);
    }

    void analyze(const Event& event) {
      
      std::pair<HepMC::GenParticle*,HepMC::GenParticle*> beam_pair = event.genEvent()->beam_particles();
      std::cout << "PDG ID 1: " << (std::get<0>(beam_pair))->pdg_id() << std::endl;
      std::cout << "PDG ID 2: " << (std::get<1>(beam_pair))->pdg_id() << std::endl;
      bool is_heavy_ion = ((std::get<0>(beam_pair))->pdg_id() == 1000822080) && ((std::get<1>(beam_pair))->pdg_id() == 1000822080);
      if (is_heavy_ion) {
	std::cout << "HI EVENT ";
      } else {
	std::cout << "PP EVENT ";
      }
      std::cout << event.genEvent()->heavy_ion();
      
      const ChargedFinalState& triggerFinalState = applyProjection<ChargedFinalState>(event, "CFSTrigger");
      Particles triggerParticles = triggerFinalState.particlesByPt();
      std::cout << "trig: " << triggerParticles.size();
      
      ChargedFinalState associatedFinalState[PT_BINS];
      Particles associatedParticles[PT_BINS];
      std::cout << ", assoc: ";
      for (int ipt = 0; ipt < PT_BINS; ipt++) {
	associatedFinalState[ipt] = applyProjection<ChargedFinalState>(event, "CFSAssoc" + std::to_string(ipt));
	associatedParticles[ipt] = associatedFinalState[ipt].particlesByPt();
	std::cout << associatedParticles[ipt].size() << " ";
      }
      std::cout << std::endl;
      
      // Check event type
      string event_string = "";
      if (is_heavy_ion) {
	double centr = centrality(event);
	std::cout << "centrality: " << centr << std::endl;
	if (centr > 0.0 && centr < 5.0) {
	  event_type = 1; // PbPb, central
	  event_string = "PbPb central";
	}
	else if (centr > 60.0 && centr < 90.0) {
	  event_type = 2; // PbPb, peripherial
	  event_string = "PbPb peripheral";
	}
	else {
	  event_type = 3; // PbPb, other
	  event_string = "PbPb other";
	}
      } else {
	event_type = 0; // pp
	event_string = "pp";
      }
      std::cout << "ev type: " << event_string << " (" << event_type << ")" << std::endl;
      
      if (event_type == 3)
	vetoEvent;
      
      _histTriggerCounter->fill(event_type, triggerParticles.size());
      
      // Loop over trigger particles
      foreach (const Particle& triggerParticle, triggerParticles) {
	std::cout << "Trigger particle" << std::endl;
	for (int ipt = 0; ipt < PT_BINS; ipt++) {
	  std::cout << "pt bin " << ipt << std::endl;
	  /*if (associatedParticles[ipt].size() > 1 || 
	      (associatedParticles[ipt].size() == 1 && associatedParticles[ipt].at(0) != triggerParticle)) {
	    num_trig[ipt] += 1;
	    _histTriggerCounter[event_type]->fill(ipt, 1);
	    std::cout << "Found some associated particles" << std::endl;*/
	  /*foreach (const Particle& associatedParticle, associatedParticles[ipt]) {
	    if (associatedParticle != triggerParticle) {
	      if (triggerParticle.pt() > associatedParticle.pt()) {
		_histTriggerCounter[event_type]->fill(ipt, 1);
		std::cout << "Adding trigger particle" << std::endl;
		break;
	      }
	    }
	    }*/
	  foreach (const Particle& associatedParticle, associatedParticles[ipt]) {
	    std::cout << "Associated particle" << std::endl;
	    if (associatedParticle != triggerParticle) {
	      if (triggerParticle.pt() > associatedParticle.pt()) {
		double dPhi = triggerParticle.phi() - associatedParticle.phi();
		while (dPhi > 1.5 * M_PI)  { dPhi -= 2 * M_PI; }
		while (dPhi < -0.5 * M_PI) { dPhi += 2 * M_PI; }
		std::cout << "Filling yield histogram for pt = " << pt_limits[ipt] << "-" << pt_limits[ipt + 1] << " GeV/c for delta phi = " << dPhi << std::endl;
		_histYield[event_type][ipt]->fill(dPhi, 1);
	      }
	    }
	    //}
	  }
	  //else {
	  /*if (trigger_added == false) {
	    std::cout << "No associated particles with pt = " << pt_limits[ipt] << "-" << pt_limits[ipt + 1] << " GeV/c found for that trigger particle" << std::endl;
	    }*/
	}
	
      }
      std::cout << std::endl;
      
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
      
      std::cout << "Trigger counter histogram entries: " << 
	_histTriggerCounter->bin(0).sumW() << " " << 
	_histTriggerCounter->bin(1).sumW() << " " << 
	_histTriggerCounter->bin(2).sumW() << std::endl;
      //std::cout << "pp trigger counter histogram entries: " << _histTriggerCounter[0]->numEntries() << std::endl;
      //std::cout << "Central trigger counter histogram entries: " << _histTriggerCounter[1]->numEntries() << std::endl;
      //std::cout << "Peripheral trigger counter histogram entries: " << _histTriggerCounter[2]->numEntries() << std::endl;
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
	  if (_histTriggerCounter->numEntries() == 0 || _histYield[itype][ipt]->numEntries() == 0) {
	    std::cout << _histTriggerCounter->numEntries() << std::endl;
	    std::cout << _histYield[itype][ipt]->numEntries() << std::endl;
	    std::cout << "There are no entries in one of the histograms" << std::endl;
	    continue;
	  }
	  
	  // Scale yield histogram
	  std::cout << "Scaling yield histograms..." << std::endl;
	  if (_histTriggerCounter->bin(itype).sumW() != 0) {
	    std::cout << "Scaling histogram type " << itype << " pt bin " << ipt << "..." << std::endl;
	    std::cout << "Scaling factor: " << _histTriggerCounter->bin(itype).sumW() << std::endl;
	    scale(_histYield[itype][ipt], (1. / _histTriggerCounter->bin(itype).sumW()));
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
	std::cout << "I_AA 0-5%/pp, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << nearSide[1][ipt] / nearSide[0][ipt] << std::endl;
	_histIAA[0]->fillBin(ipt, nearSide[1][ipt] / nearSide[0][ipt]);
	std::cout << "I_AA 60-90%/pp, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << nearSide[2][ipt] / nearSide[0][ipt] << std::endl;
	_histIAA[1]->fillBin(ipt, nearSide[2][ipt] / nearSide[0][ipt]);
	std::cout << "I_CP 0-5%/60-90%, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << nearSide[1][ipt] / nearSide[2][ipt] << std::endl;
	_histIAA[2]->fillBin(ipt, nearSide[1][ipt] / nearSide[2][ipt]);
      }
      std::cout << "Away side" << std::endl;
      for (int ipt = 0; ipt < PT_BINS; ipt++) {	
	std::cout << "I_AA 0-5%/pp, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << awaySide[1][ipt] / awaySide[0][ipt] << std::endl;
	_histIAA[3]->fillBin(ipt, awaySide[1][ipt] / awaySide[0][ipt]);
	std::cout << "I_AA 60-90%/pp, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << awaySide[2][ipt] / awaySide[0][ipt] << std::endl;
	_histIAA[4]->fillBin(ipt, awaySide[2][ipt] / awaySide[0][ipt]);
	std::cout << "I_CP 0-5%/60-90%, pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << awaySide[1][ipt] / awaySide[2][ipt] << std::endl;
	_histIAA[5]->fillBin(ipt, awaySide[1][ipt] / awaySide[2][ipt]);
      }
      
    }
    
  private:
    
    static const int PT_BINS = 4;
    static const int EVENT_TYPES = 3;

    Histo1DPtr _histYield[EVENT_TYPES][PT_BINS];
    Histo1DPtr _histYieldBkgRemoved[EVENT_TYPES][PT_BINS];
    Histo1DPtr _histTriggerCounter;
    Histo1DPtr _histIAA[6];
    int num_trig[PT_BINS];
    int pt_limits[5];
    int event_type;
    
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I930312);
}
