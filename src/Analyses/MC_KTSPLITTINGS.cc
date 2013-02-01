// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetYODA.hh"

namespace Rivet {


  /// @brief MC validation analysis for jet events
  class MC_KTSPLITTINGS : public MC_JetSplittings {
  public:

    MC_KTSPLITTINGS()
      : MC_JetSplittings("MC_KTSPLITTINGS", 4, "Jets")
    {    }


  public:

    void init() {
      FinalState fs;
      FastJets jetpro(fs, FastJets::KT, 0.6);
      addProjection(jetpro, "Jets");

      MC_JetSplittings::init();
    }


    void analyze(const Event& event) {
      MC_JetSplittings::analyze(event);
    }


    void finalize() {
      MC_JetSplittings::finalize();
    }

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_KTSPLITTINGS);

}