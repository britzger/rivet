// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Centrality.hh"
#include "Rivet/Tools/Cuts.hh"
#include <fstream>
//#include "Centrality.hh"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet {

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
      
      // Centrality projection for this analysis
      Centrality centrality(cfs, Centrality::Method::ImpactParameter, 10000);
      addProjection(centrality, "CENTR");
      
      // Histograms and variables initialization. Do it for each centrality range
      _histNchVsCentr = bookHisto1D(1, 1, 1);
      
    }

    void analyze(const Event& event) {

      const ChargedFinalState& charged = applyProjection<ChargedFinalState>(event, "CFS");
      Particles chargedParticles = charged.particlesByPt();
      
      const Centrality& centralityProjection = applyProjection<Centrality>(event, "CENTR");
      if (!centralityProjection.isValid())
	vetoEvent;
      
      const double c = centralityProjection.getCentrality();
      
      cout << "Centrality: " << c << endl;
      if ((c < 0.) || (c > 80.))
	vetoEvent;
      
      _histNchVsCentr->fill(c, charged.particles().size() * event.weight());
      
    }

    void finalize() {
      // Finishing the centrality calculations
      //const Centrality& centralityProjection = getProjection<Centrality>("CENTR");
      //centralityProjection.finalize();
      // Right scaling of the histograms with their individual weights.
      for (size_t ii = 0; ii < _histNchVsCentr->numBins(); ii++) {
	double scale = _histNchVsCentr->bin(ii).xWidth() / _histNchVsCentr->bin(ii).numEntries();
	_histNchVsCentr->bin(ii).scaleW( scale );
      }
    }


    // @note First time a post() method is implemented and used by Rivet. This is only available on the
    // mcplots-alice-dev.cern.ch and mcplots-alice-dev2.cern.ch machine because the Rivet installation here 
    // is the only one which implements Rivet::Analysis::post() and Rivet::AnalysisHandler::post()
    void post() { }

  private:
    
    Histo1DPtr _histNchVsCentr;
    
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2010_I880049);
}
