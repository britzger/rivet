// -*- C++ -*-
#ifndef RIVET_ZFinder_HH
#define RIVET_ZFinder_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Convenience finder of leptonically decaying Zs
  ///
  /// Chain together different projections as convenience for finding Z's
  /// from two leptons in the final state, including photon clustering.
  ///
  /// @todo Inherit directly from ParticleFinder, not FinalState
  /// @todo Alias then rename as Dileptons
  class ZFinder : public FinalState {
  public:

    enum ClusterPhotons { NOCLUSTER=0, CLUSTERNODECAY=1, CLUSTERALL };
    enum PhotonTracking { NOTRACK=0, TRACK=1 };

    /// @name Constructors
    //@{

    /// Constructor taking cuts object
    /// @param inputfs Input final state
    /// @param cuts lepton cuts
    /// @param pid type of the leptons
    /// @param minmass,maxmass mass window
    /// @param dRmax maximum dR of photons around leptons to take into account
    ///  for Z reconstruction (only relevant if one of the following are true)
    /// @param clusterPhotons whether such photons are supposed to be
    ///  clustered to the lepton objects and thus Z mom
    /// @param trackPhotons whether such photons should be added to _theParticles
    ///  (cf. _trackPhotons)
    ZFinder(const FinalState& inputfs,
            const Cut & cuts,
            PdgId pid,
            double minmass, double maxmass,
            double dRmax=0.1,
            ClusterPhotons clusterPhotons=CLUSTERNODECAY,
            PhotonTracking trackPhotons=NOTRACK,
            double masstarget=91.2*GeV);

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(ZFinder);

    //@}


    /// Access to the found bosons
    ///
    /// @note Currently either 0 or 1 boson can be found.
    const Particles& bosons() const { return _bosons; }

    /// Access to the found boson (assuming it exists).
    const Particle boson() const { return _bosons[0]; }

    /// Access to the Z constituent clustered leptons
    ///
    /// For example, to make more fine-grained cuts on the clustered leptons.
    /// The positive charge constituent is first in the list (if not empty), and
    /// the negative one second.
    const Particles& constituents() const { return _constituents; }

    /// Access to the Z constituent clustered leptons, sorted by a comparison functor
    ///
    /// Unlike the no-arg version, this returns by value (i.e. is less efficient)
    template <typename CMP>
    Particles constituents(const CMP& cmp) const {
      Particles rtn = constituents();
      std::sort(rtn.begin(), rtn.end(), cmp);
      return rtn;
    }

    /// Access to all DressedLeptons in the fiducial region.
    ///
    /// This includes those DressedLeptons that could not
    /// be paired up with any other DressedLepton to form a Z candidate.
    const vector<DressedLepton>& allLeptons() const { return _allLeptons; }

    /// Access to the particles other than the Z leptons and clustered photons
    ///
    /// Useful for e.g. input to a jet finder
    const VetoedFinalState& remainingFinalState() const;


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  public:

    /// Clear the projection
    void clear() {
      _theParticles.clear();
      _bosons.clear();
      _constituents.clear();
      _allLeptons.clear();
    }


  private:
    /// Mass cuts to apply to clustered leptons (cf. InvMassFinalState)
    double _minmass, _maxmass, _masstarget;

    /// Switch for tracking of photons (whether to add them to _theParticles)
    /// This is relevant when the ZFinder::_theParticles are to be excluded
    /// from e.g. the input to a jet finder, to specify whether the clustered
    /// photons are to be excluded as well.
    /// (Yes, some experiments make a difference between clusterPhotons and
    /// trackPhotons!)
    PhotonTracking _trackPhotons;

    /// Lepton flavour
    PdgId _pid;

    /// list of found bosons (currently either 0 or 1)
    Particles _bosons;

    /// Clustered leptons
    vector<Particle> _constituents;

    ///Dressed leptons
    vector<DressedLepton> _allLeptons;

  };


}



#endif
