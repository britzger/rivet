// -*- C++ -*-
#include "Rivet/Projections/DISLepton.hh"

namespace Rivet {


  int DISLepton::compare(const Projection& p) const {
    const DISLepton& other = pcast<DISLepton>(p);
    return mkNamedPCmp(other, "Beam") || mkNamedPCmp(other, "LFS") ||
      mkNamedPCmp(other, "IFS");
  }


  void DISLepton::project(const Event& e) {

    // Find incoming lepton beam
    const ParticlePair& inc = applyProjection<Beam>(e, "Beam").beams();
    bool firstIsLepton = PID::isLepton(inc.first.pid());
    bool secondIsLepton = PID::isLepton(inc.second.pid());
    if (firstIsLepton && !secondIsLepton) {
      _incoming = inc.first;
    } else if (!firstIsLepton && secondIsLepton) {
      _incoming = inc.second;
    } else {
      fail();
      return;
    }

    // If no graph-connected scattered lepton, use the hardest
    // (preferably same-flavour) prompt FS lepton in the event.

    const Particles fsleptons =
      applyProjection<FinalState>(e, "LFS").particles(isLepton, cmpMomByE);
    Particles sfleptons =
      filter_select(fsleptons, Cuts::pid == _incoming.pid());
    MSG_DEBUG("SF leptons = " << sfleptons.size() << ", all leptons = "
              << fsleptons.size());
    if ( sfleptons.empty() ) sfleptons = fsleptons;

    if ( _isolDR > 0.0 ) {
      const Particles & other =
        applyProjection<FinalState>(e, "IFS").particles();
      while (!sfleptons.empty()) {
        bool skip = false;
        Particle testlepton = sfleptons.front();
        for ( auto p: other ) {
          if ( skip ) break;
          if ( deltaR(p, testlepton) < _isolDR ) skip = true;
          for ( auto c : testlepton.constituents() ) {
            if ( c.genParticle() == p.genParticle() ) {
              skip = false;
              break;
            }
          }
        }
        if ( !skip ) break;
        sfleptons.erase(sfleptons.begin());
      }
    }

    if ( !sfleptons.empty() ) {
      _outgoing = sfleptons.front();
    } else {
      fail();
    }

  }


}
