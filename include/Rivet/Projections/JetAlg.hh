// -*- C++ -*-
#ifndef RIVET_JetAlg_HH
#define RIVET_JetAlg_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

namespace Rivet {


  /// Abstract base class for projections which can return a set of {@link Jet}s.
  class JetAlg : public Projection {
  public:

    /// Constructor
    JetAlg(const FinalState& fs);

    JetAlg() {};

    /// Clone on the heap.
    virtual const Projection* clone() const = 0;

    /// Destructor
    virtual ~JetAlg() { }

    /// @brief Include invisible particles in jet construction.
    /// The default behaviour is that jets are only constructed from visible
    /// (i.e. charged under an SM gauge group) particles. Some jet studies,
    /// including those from ATLAS, use a definition in which neutrinos from hadron
    /// decays are included (via MC correction) in the experimental jet definition.
    /// Setting this flag to true avoids the automatic restriction to a VisibleFinalState.
    void useInvisibles(bool useinvis=true) {
      _useInvisibles = useinvis;
    }


    /// Get jets in no guaranteed order, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts
    /// @todo Can't this be a const Cut& arg?
    virtual Jets jets(const Cut & c = Cuts::open()) const {
      const Jets rawjets = _jets(0.0); // arg means no pT cut
      // Just return a copy of rawjets if the cut is open
      if (c == Cuts::open()) return rawjets;
      // If there is a non-trivial cut...
      Jets rtn;
      rtn.reserve(size());
      foreach (const Jet& j, rawjets)
        if (c->accept(j)) rtn.push_back(j);
      return rtn;
    }

    /// Get the jets, ordered by supplied sorting function object, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    template <typename F>
    Jets jets(F sorter, const Cut & c = Cuts::open()) const {
      /// @todo Will the vector be efficiently std::move'd by value through this function chain?
      return sortBy(jets(c), sorter);
    }

    /// Get the jets, ordered by supplied sorting function object, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    template <typename F>
    Jets jets(const Cut & c ,  F sorter) const {
      /// @todo Will the vector be efficiently std::move'd by value through this function chain?
      return sortBy(jets(c), sorter);
    }


    /// Get the jets, ordered by \f$ p_T \f$, with optional cuts.
    ///
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    ///
    /// This is a very common use-case, so is available as syntatic sugar for jets(c, cmpMomByPt).
    /// @todo The other sorted accessors should be removed in a cleanup.
    Jets jetsByPt(const Cut & c = Cuts::open()) const {
      return jets(c, cmpMomByPt);
    }


    /// @name Old sorted jet accessors
    /// @deprecated Use the versions with sorter function arguments
    //@{

    /// Get the jets, ordered by \f$ |p| \f$, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @deprecated Use the version with a sorter function argument.
    DEPRECATED("Use the version with a sorter function argument.")
    Jets jetsByP(const Cut & c = Cuts::open()) const {
      return jets(c, cmpMomByP);
    }

    /// Get the jets, ordered by \f$ E \f$, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @deprecated Use the version with a sorter function argument.
    DEPRECATED("Use the version with a sorter function argument.")
    Jets jetsByE(const Cut & c = Cuts::open()) const {
      return jets(c, cmpMomByE);
    }

    /// Get the jets, ordered by \f$ E_T \f$, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @deprecated Use the version with a sorter function argument.
    DEPRECATED("Use the version with a sorter function argument.")
    Jets jetsByEt(const Cut & c = Cuts::open()) const {
      return jets(c, cmpMomByEt);
    }

    //@}


    /// @name Old jet accessors
    /// @deprecated Use the versions with Cut arguments
    //@{

    /// Get jets in no guaranteed order, with optional cuts on \f$ p_\perp \f$ and rapidity.
    ///
    /// @deprecated Use the version with a Cut argument
    /// @note Returns a copy rather than a reference, due to cuts
    DEPRECATED("Use the version with a Cut argument.")
    Jets jets(double ptmin, double ptmax=MAXDOUBLE,
              double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
              RapScheme rapscheme=PSEUDORAPIDITY) const {
      if (rapscheme == PSEUDORAPIDITY) {
        return jets((Cuts::pT >= ptmin) & (Cuts::pT < ptmax) & (Cuts::rapIn(rapmin, rapmax)));
      } else if (rapscheme == RAPIDITY) {
        return jets((Cuts::pT >= ptmin) & (Cuts::pT < ptmax) & (Cuts::etaIn(rapmin, rapmax)));
      }
      throw LogicError("Unknown rapidity scheme. This shouldn't be possible!");
    }

    /// Get the jets, ordered by \f$ p_T \f$, with a cut on \f$ p_\perp \f$.
    ///
    /// @deprecated Use the version with a Cut argument
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    ///
    /// This is a very common use-case, so is available as syntatic sugar for jets(Cuts::pT >= ptmin, cmpMomByPt).
    /// @todo The other sorted accessors should be removed in a cleanup.
    Jets jetsByPt(double ptmin) const {
      return jets(Cuts::pT >= ptmin, cmpMomByPt);
    }

    //@}


  protected:

    /// @brief Internal pure virtual method for getting jets in no guaranteed order.
    /// An optional cut on min \f$ p_\perp \f$ is applied in this function, since that is
    /// directly supported by FastJet and it seems a shame to not make use of that. But
    /// all other jet cuts are applied at the @c ::jets() function level.
    /// @todo Remove the ptmin cut
    virtual Jets _jets(double ptmin=0) const = 0;


  public:

    /// Number of jets.
    virtual size_t size() const = 0;
    /// Determine if the jet collection is empty.
    bool empty() const { return size() != 0; }

    /// Clear the projection.
    virtual void reset() = 0;

    typedef Jet entity_type;
    typedef Jets collection_type;

    /// Template-usable interface common to FinalState.
    collection_type entities() const { return jets(); }

    /// Do the calculation locally (no caching).
    virtual void calc(const Particles& constituents, const Particles& tagparticles=Particles()) = 0;


  protected:

    /// Perform the projection on the Event.
    virtual void project(const Event& e) = 0;

    /// Compare projections.
    virtual int compare(const Projection& p) const = 0;


  protected:

    /// Flag to determine whether or not the VFS wrapper is to be used.
    bool _useInvisibles;

  };


}

#endif
