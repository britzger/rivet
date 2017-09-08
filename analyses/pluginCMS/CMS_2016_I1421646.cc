// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CMS azimuthal decorrelations at 8 TeV
  class CMS_2016_I1421646 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1421646);


    /// Book projections and histograms
    void init() {

      FastJets akt(FinalState(), FastJets::ANTIKT, 0.7);
      addProjection(akt, "antikT");

      {Histo1DPtr tmp; _h_deltaPhi.addHistogram( 200.,  300., book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.addHistogram( 300.,  400., book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.addHistogram( 400.,  500., book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.addHistogram( 500.,  700., book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.addHistogram( 700.,  900., book(tmp, 5, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.addHistogram( 900.,  1100., book(tmp, 6, 1, 1));}
      {Histo1DPtr tmp; _h_deltaPhi.addHistogram( 1100., 4000., book(tmp, 7, 1, 1));}
    }


    /// Per-event analysis
    void analyze(const Event & event) {

      const Jets& jets = apply<JetAlg>(event, "antikT").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 100*GeV);
      if (jets.size() < 2) vetoEvent;
      if (jets[0].pT() < 200*GeV) vetoEvent;

      const double dphi = deltaPhi(jets[0].phi(), jets[1].phi());
      _h_deltaPhi.fill(jets[0].pT(), dphi, 1.0);
    }


    /// Scale histograms
    void finalize() {
      for (Histo1DPtr histo : _h_deltaPhi.histos()) normalize(histo);
    }


  private:

    BinnedHistogram<double> _h_deltaPhi;

  };


  // A hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1421646);

}
