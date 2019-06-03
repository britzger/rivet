// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include <cmath>

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief OPAL multiplicities at various energies
  /// @author Peter Richardson
  class OPAL_2002_S5361494 : public Analysis {
  public:

    /// Constructor
    OPAL_2002_S5361494()
      : Analysis("OPAL_2002_S5361494")
    {}

    /// @name Analysis methods
    //@{


    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "CFS");
      declare(InitialQuarks(), "IQF");

      _cLight  = bookCounter("TMP/CLIGHT" );
      _wLight  = bookCounter("TMP/WLIGHT" );
      _cCharm  = bookCounter("TMP/CCHARM" );
      _wCharm  = bookCounter("TMP/WCHARM" );
      _cBottom = bookCounter("TMP/CBOTTOM");
      _wBottom = bookCounter("TMP/WBOTTOM");

      book(_s_hists[0], 1, 1, 1);
      book(_s_hists[1], 1, 1, 2);
      book(_s_hists[2], 1, 1, 3);
      book(_s_hists[3], 1, 1, 4); //< bottom minus light
    }


    void analyze(const Event& event) {
       // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;


      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      }
      else {
        map<int, double> quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }
      const size_t numParticles = cfs.particles().size();
      switch (flavour) {
      case 1: case 2: case 3:
        _wLight ->fill();
        _cLight ->fill(numParticles);
        break;
      case 4:
        _wCharm ->fill();
        _cCharm ->fill(numParticles);
        break;
      case 5:
        _wBottom->fill();
        _cBottom->fill(numParticles);
        break;
      }

    }


    void finalize() {
      // calculate the averages and diffs
      if(_wLight ->numEntries()) scale( _cLight, 1./_wLight->val());
      if(_wCharm ->numEntries()) scale( _cCharm, 1./_wCharm->val());
      if(_wBottom->numEntries()) scale(_cBottom,1./_wBottom->val());
      Counter _cDiff = *_cBottom - *_cLight;

      // fill the histograms
      for (unsigned int ix=1;ix<5;++ix) {
	double val(0.), err(0.0);
	if(ix==1) {
	  val = _cBottom->val();
	  err = _cBottom->err();
	}
	else if(ix==2) {
	  val = _cCharm->val();
	  err = _cCharm->err();
	}
	else if(ix==3) {
	  val = _cLight->val();
	  err = _cLight->err();
	}
	else if(ix==4) {
	  val = _cDiff.val();
	  err = _cDiff.err();
	}

        /// @todo TIDY!
	Scatter2D temphisto(refData(1, 1, ix));
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	    s_hists[ix]->addPoint(x, val, ex, make_pair(err,err));
	  } else {
	    s_hists[ix]->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }

    //@}


  private:

    /// Final plots
    Scatter2DPtr s_hists;

    /// @name Multiplicities
    /// @todo Don't we have a Dbn1D-like type that can do both at once?
    //@{
    CounterPtr _cLight, _wLight;
    CounterPtr _cCharm, _wCharm;
    CounterPtr _cBottom, _wBottom;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2002_S5361494);

}
