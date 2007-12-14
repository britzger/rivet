// -*- C++ -*-
#ifndef RIVET_VetoedFinalState_HH
#define RIVET_VetoedFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {

  /// Project all final state particles except for those listed by PDG code.
  class VetoedFinalState : public FinalState {

  public:
    
    /// Typedef for a vetoing entry.
    typedef map<long, BinaryCut> VetoDetails;

    /// The default constructor. Must specify a FinalState projection 
    /// object which is assumed to live through the run.
    VetoedFinalState(FinalState& fsp)
      : _fsproj(fsp)
    {
      addProjection(_fsproj);
    }

    /// You can add a map of ID plus a pair containing \f$ p_{Tmin} \f$ and
    /// \f$ p_{Tmax} \f$ - these define the range of particles to be vetoed.
    VetoedFinalState(FinalState& fsp, const VetoDetails& vetocodes)
      : _fsproj(fsp), _vetoCodes(vetocodes)
    {
      addProjection(_fsproj);
    }
    

  public:
    /// Return the name of the projection
    string getName() const {
      return "VetoedFinalState";
    }

  protected:
    
    /// Apply the projection on the supplied event.
    void project(const Event& e);
    
    /// Compare projections.
    int compare(const Projection& p) const;
    
  public:
    
    /// Get the list of particle IDs and \f$ p_T \f$ ranges to veto.
    const VetoDetails& getVetoDetails() const {
      return _vetoCodes;
    }
  
    /// Add a particle ID and \f$ p_T \f$ range to veto. Particles with \f$ p_T \f$ 
    /// IN the given range will be rejected.
    VetoedFinalState& addVetoDetail(const long id, const double ptmin, const double ptmax) {
      BinaryCut ptrange(ptmin, ptmax);
      _vetoCodes.insert(make_pair(id, ptrange));
      return *this;
    }

    /// Add a particle/antiparticle pair to veto in a given \f$ p_T \f$ range. Given a single ID, both
    /// the particle and its conjugate antiparticle will be rejected if their \f$ p_T \f$ is IN the given range.
    VetoedFinalState& addVetoPairDetail(const long id, const double ptmin, const double ptmax) {
      addVetoDetail(id,  ptmin, ptmax);
      addVetoDetail(-id, ptmin, ptmax);
      return *this;
    }

    /// Add a particle/antiparticle pair to veto. Given a single ID, both the particle and its corresponding 
    /// antiparticle (for all \f$ p_T \f$ values) will be vetoed.
    VetoedFinalState& addVetoPairId(const long id) {
      addVetoId(id);
      addVetoId(-id);
      return *this;
    }

    /// Add a particle ID to veto (all \f$ p_T \f$ range will be vetoed).
    VetoedFinalState& addVetoId(const long id) {
      BinaryCut ptrange(0.0, numeric_limits<double>::max());
      _vetoCodes.insert(make_pair(id, ptrange));
      return *this;
    }

    /// Set the list of particle IDs and \f$ p_T \f$ ranges to veto.
    VetoedFinalState& setVetoDetails(const VetoDetails& ids) {
      _vetoCodes = ids;
      return *this;
    }

    /// Clear the list of particle IDs and ranges to veto.
    VetoedFinalState& clearVetoDetails() {
      _vetoCodes.clear();
      return *this;
    }
    
  private:
    
    /// The projector for the full final state.
    FinalState _fsproj;
    
    /// The final-state particles.
    VetoDetails _vetoCodes;
  
  };

  
}


#endif
