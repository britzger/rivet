// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief D0 Run II Z \f$ p_\perp \f$ differential cross-section shape
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  /// @author Frank Siegert
  class D0_2008_S7554427 : public Analysis {

  public:

    /// Default constructor.
    D0_2008_S7554427() : Analysis("D0_2008_S7554427")
    {
      // Run II Z pT
    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      ZFinder zfinder(-MAXRAPIDITY, MAXRAPIDITY, 0.0*GeV, ELECTRON,
                      40.0*GeV, 200.0*GeV, 0.2, true, true);
      addProjection(zfinder, "ZFinder");

      _h_ZpT         = bookHistogram1D(1, 1, 1);
      _h_forward_ZpT = bookHistogram1D(3, 1, 1);
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = e.weight();

      const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() == 1) {
        double yZ = fabs(zfinder.bosons()[0].momentum().rapidity());
        double pTZ = zfinder.bosons()[0].momentum().pT();
        _h_ZpT->fill(pTZ, weight);
        if (yZ > 2.0) {
          _h_forward_ZpT->fill(pTZ, weight);
        }
      }
      else {
        getLog() << Log::DEBUG << "no unique lepton pair found." << endl;
      }

    }



    // Finalize
    void finalize() {
      normalize(_h_ZpT);
      normalize(_h_forward_ZpT);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_ZpT;
    AIDA::IHistogram1D * _h_forward_ZpT;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2008_S7554427);

}
