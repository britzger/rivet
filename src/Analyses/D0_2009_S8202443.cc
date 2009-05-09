// -*- C++ -*-
#include "Rivet/Analyses/D0_2009_S8202443.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2009_S8202443::D0_2009_S8202443() :
      _sum_of_weights(0.0), _sum_of_weights_constrained(0.0)
  {
    setBeams(PROTON, ANTIPROTON);
    
    /// @todo Use cross-section from generator
    //setNeedsCrossSection(true);

    // Leptons in constrained tracking acceptance
    std::vector<std::pair<double, double> > etaRanges;
    etaRanges.push_back(make_pair(-2.5, -1.5));
    etaRanges.push_back(make_pair(-1.1, 1.1));
    etaRanges.push_back(make_pair(1.5, 2.5));
    ZFinder zfinder_constrained(etaRanges, 25.0*GeV, ELECTRON,
                                65.0*GeV, 115.0*GeV, 0.2);
    addProjection(zfinder_constrained, "ZFinderConstrained");
    FastJets conefinder_constrained(zfinder_constrained.remainingFinalState(),
                                    FastJets::D0ILCONE, 0.5, 20.0*GeV);
    addProjection(conefinder_constrained, "ConeFinderConstrained");

    // Unconstrained leptons
    ZFinder zfinder(FinalState(), ELECTRON, 65.0*GeV, 115.0*GeV, 0.2);
    addProjection(zfinder, "ZFinder");
    FastJets conefinder(zfinder.remainingFinalState(), FastJets::D0ILCONE, 0.5, 20.0*GeV);
    addProjection(conefinder, "ConeFinder");
  } 


  // Book histograms
  void D0_2009_S8202443::init() {
    _h_jet1_pT_constrained = bookHistogram1D(
        1, 1, 1, "pT of 1st jet (constrained electrons)", "$p_{\\perp}^{\\text{1st jet}}$ [GeV]",
        "$1/\\sigma \\; \\text{d}\\sigma/\\text{d}p_{\\perp}^{\\text{1st jet}}$");
    _h_jet2_pT_constrained = bookHistogram1D(
        3, 1, 1, "pT of 2nd jet (constrained electrons)", "$p_{\\perp}^{\\text{2nd jet}}$ [GeV]",
        "$1/\\sigma \\; \\text{d}\\sigma/\\text{d}p_{\\perp}^{\\text{2nd jet}}$");
    _h_jet3_pT_constrained = bookHistogram1D(
        5, 1, 1, "pT of 3rd jet (constrained electrons)", "$p_{\\perp}^{\\text{3rd jet}}$ [GeV]",
        "$1/\\sigma \\; \\text{d}\\sigma/\\text{d}p_{\\perp}^{\\text{3rd jet}}$");
    _h_jet1_pT = bookHistogram1D(
        2, 1, 1, "pT of 1st jet", "$p_{\\perp}^{\\text{1st jet}}$ [GeV]",
        "$1/\\sigma \\; \\text{d}\\sigma/\\text{d}p_{\\perp}^{\\text{1st jet}}$");
    _h_jet2_pT = bookHistogram1D(
        4, 1, 1, "pT of 2nd jet", "$p_{\\perp}^{\\text{2nd jet}}$ [GeV]",
        "$1/\\sigma \\; \\text{d}\\sigma/\\text{d}p_{\\perp}^{\\text{2nd jet}}$");
    _h_jet3_pT = bookHistogram1D(
        6, 1, 1, "pT of 3rd jet", "$p_{\\perp}^{\\text{3rd jet}}$ [GeV]",
        "$1/\\sigma \\; \\text{d}\\sigma/\\text{d}p_{\\perp}^{\\text{3rd jet}}$");
  }



  // Do the analysis 
  void D0_2009_S8202443::analyze(const Event & e) {
    double weight = e.weight();

    // unconstrained electrons first
    const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
    if (zfinder.particles().size()==1) {
      _sum_of_weights += weight;
      const JetAlg& jetpro = applyProjection<JetAlg>(e, "ConeFinder");
      const Jets& jets = jetpro.jetsByPt(20.0*GeV);
      Jets jets_cut;
      foreach (const Jet& j, jets) {
        if (fabs(j.momentum().pseudorapidity()) < 2.5) {
          jets_cut.push_back(j);
        }
      }
      
      if (jets_cut.size()>0) {
        _h_jet1_pT->fill(jets_cut[0].momentum().pT()/GeV, weight);
      }
      if (jets_cut.size()>1) {
        _h_jet2_pT->fill(jets_cut[1].momentum().pT()/GeV, weight);
      }
      if (jets_cut.size()>2) {
        _h_jet3_pT->fill(jets_cut[2].momentum().pT()/GeV, weight);
      }
    }
    else {
      getLog() << Log::DEBUG << "no unique lepton pair found." << endl;
    }


    // constrained electrons
    const ZFinder& zfinder_constrained = applyProjection<ZFinder>(e, "ZFinderConstrained");
    if (zfinder_constrained.particles().size()==1) {
      _sum_of_weights_constrained += weight;
      const JetAlg& jetpro = applyProjection<JetAlg>(e, "ConeFinderConstrained");
      const Jets& jets = jetpro.jetsByPt(20.0*GeV);
      Jets jets_cut;
      foreach (const Jet& j, jets) {
        if (fabs(j.momentum().pseudorapidity()) < 2.5) {
          jets_cut.push_back(j);
        }
      }
      
      if (jets_cut.size()>0) {
        _h_jet1_pT_constrained->fill(jets_cut[0].momentum().pT()/GeV, weight);
      }
      if (jets_cut.size()>1) {
        _h_jet2_pT_constrained->fill(jets_cut[1].momentum().pT()/GeV, weight);
      }
      if (jets_cut.size()>2) {
        _h_jet3_pT_constrained->fill(jets_cut[2].momentum().pT()/GeV, weight);
      }
    }
    else {
      getLog() << Log::DEBUG << "no unique lepton pair found." << endl;
      vetoEvent(e);
    }
  }



  // Finalize
  void D0_2009_S8202443::finalize() {
    scale(_h_jet1_pT, 1.0/_sum_of_weights);
    scale(_h_jet2_pT, 1.0/_sum_of_weights);
    scale(_h_jet3_pT, 1.0/_sum_of_weights);
    scale(_h_jet1_pT_constrained, 1.0/_sum_of_weights_constrained);
    scale(_h_jet2_pT_constrained, 1.0/_sum_of_weights_constrained);
    scale(_h_jet3_pT_constrained, 1.0/_sum_of_weights_constrained);
  }

}
