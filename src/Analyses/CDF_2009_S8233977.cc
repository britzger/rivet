// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2009_S8233977.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  CDF_2009_S8233977::CDF_2009_S8233977()
    : Analysis("CDF_2009_S8233977")
  { 
    setBeams(PROTON, ANTIPROTON);
    const FinalState fs(-1.0, 1.0, 0.0*GeV);
    const ChargedFinalState cfs(-1.0, 1.0, 0.4*GeV);
    addProjection(fs, "FS");
    addProjection(cfs, "CFS");
    setNeedsCrossSection(true);
  }


  // Book histograms
  void CDF_2009_S8233977::init() {
    _hist_pt_vs_multiplicity = bookProfile1D(1, 1, 1, "Mean track $p_T$ vs multiplicity",
        "$N_\\text{ch}$", "$\\langle p_T \\rangle$ / GeV");
    _hist_pt                 = bookHistogram1D(2, 1, 1, "track $p_T$",
        "$p_T$ / GeV", "$\\text{d}^3 \\sigma / p_T \\text{d}p_T \\text{d}y \\text{d}\\phi$");
    _hist_sumEt              = bookHistogram1D(3, 1, 1, "$\\sum E_T$",
        "$\\sum E_T$ / GeV", "$\\text{d}^3 \\sigma / \\text{d}E_T \\text{d}\\eta \\text{d}\\phi$");
  }


  // Do the analysis
  void CDF_2009_S8233977::analyze(const Event& e) {
    Log log = getLog();

    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(e, "CFS");
    const size_t numParticles = cfs.particles().size();

    // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
    if (numParticles < 1) {
      getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
      vetoEvent;
    }

    // Get the event weight
    const double weight = e.weight();

    foreach (const Particle& p, cfs.particles()) {
      const double pT = p.momentum().pT() / GeV;
      _hist_pt_vs_multiplicity->fill(numParticles, pT, weight);

      // The weight for entries in the pT distribution should be weight/(pT*dPhi*dy).
      //
      // - dPhi = 2*PI
      //
      // - dy depends on the pT: They calculate y assuming the particle has the
      //   pion mass and assuming that eta=1:
      //   dy = 2 * 1/2 * ln [(sqrt(m^2 + (a+1)*pT^2) + a*pT) / (sqrt(m^2 + (a+1)*pT^2) - a*pT)]
      //   with a = sinh(1).
      //
      // sinh(1) = 1.1752012
      // m(charged pion)^2 = (139.57 MeV)^2 = 0.019479785 GeV^2

      //// FIXME: The pT and sum(ET) distributions look slightly different from
      ////        Niccolo's Monte Carlo plots. Still waiting for his answer.
      const double sinh1 = 1.1752012;
      const double apT  = sinh1 * pT;
      const double mPi = 139.57*MeV;
      const double root = sqrt(mPi*mPi + (1+sinh1)*pT*pT);
      const double dy = std::log((root+apT)/(root-apT));
      const double dphi = 2*M_PI;
      _hist_pt->fill(pT, weight/(pT*dphi*dy));
    }
    double sumEt = 0.;
    foreach (const Particle& p, fs.particles()) {
      sumEt += p.momentum().Et();
    }
    _hist_sumEt->fill(sumEt, weight);
  }


  void CDF_2009_S8233977::finalize() {
    // dphi * deta = 2*PI * 2
    //// FIXME: We are normalizing to the data instead of MC cross-section
    //scale(_hist_sumEt, crossSection()/millibarn/(4*M_PI*sumOfWeights()));
    //scale(_hist_pt, crossSection()/barn*1e-24/sumOfWeights());
    normalize(_hist_sumEt, 3.530);
    normalize(_hist_pt, 2.513e-26);
  }


}
