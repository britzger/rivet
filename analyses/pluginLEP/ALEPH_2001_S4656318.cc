// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief DELPHI b-fragmentation measurement
  /// @author Hendrik Hoeth
  class ALEPH_2001_S4656318 : public Analysis {
  public:

    /// Constructor
   ALEPH_2001_S4656318()
      : Analysis("ALEPH_2001_S4656318")
    {
    }


    /// @name Analysis methods
    //@{

    /// Book projections and histograms
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      _histXbweak     = bookHisto1D(1, 1, 1);
      _histXbprim     = bookHisto1D(1, 1, 2);
      _histMeanXbweak = bookProfile1D(7, 1, 1);
      _histMeanXbprim = bookProfile1D(7, 1, 2);
    }


    void analyze(const Event& e) {
      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed ncharged cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");

      // Get event weight for histo filling
      const double weight = e.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);


      foreach (const GenParticle* p, particles(e.genEvent())) {
        const GenVertex* pv = p->production_vertex();
        const GenVertex* dv = p->end_vertex();
        if (PID::isBottomHadron(p->pdg_id())) {
          const double xp = p->momentum().e()/meanBeamMom;

          // If the B-hadron has a parton as parent, call it primary B-hadron:
          if (pv) {
            bool is_primary = false;
            for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin(); pp != pv->particles_in_const_end() ; ++pp) {
              if (PID::isParton((*pp)->pdg_id())) is_primary = true;
            }
            if (is_primary) {
              _histXbprim->fill(xp, weight);
              _histMeanXbprim->fill(_histMeanXbprim->bin(0).xMid(), xp, weight);
            }
          }

          // If the B-hadron has no B-hadron as a child, it decayed weakly:
          if (dv) {
            bool is_weak = true;
            for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin() ;
                 pp != dv->particles_out_const_end() ; ++pp) {
              if (PID::isBottomHadron((*pp)->pdg_id())) {
                is_weak = false;
              }
            }
            if (is_weak) {
              _histXbweak->fill(xp, weight);
              _histMeanXbweak->fill(_histMeanXbweak->bin(0).xMid(), xp, weight);
            }
          }

        }
      }
    }


    // Finalize
    void finalize() {
      normalize(_histXbprim);
      normalize(_histXbweak);
    }


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.

    Histo1DPtr _histXbprim;
    Histo1DPtr _histXbweak;

    Profile1DPtr _histMeanXbprim;
    Profile1DPtr _histMeanXbweak;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALEPH_2001_S4656318);

}