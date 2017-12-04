// -*- C++ -*-
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  // Separate-FS version
  DressedLeptons::DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
                                 double dRmax, const Cut& cut, bool useDecayPhotons)
    : FinalState(cut),
      _dRmax(dRmax), _fromDecay(useDecayPhotons)
  {
    setName("DressedLeptons");

    IdentifiedFinalState photonfs(photons, PID::PHOTON);
    addProjection(photonfs, "Photons");

    IdentifiedFinalState leptonfs(bareleptons);
    leptonfs.acceptIdPairs({PID::ELECTRON, PID::MUON, PID::TAU});
    addProjection(leptonfs, "Leptons");
  }


  // Single-FS version
  DressedLeptons::DressedLeptons(const FinalState& barefs,
                                 double dRmax, const Cut& cut, bool useDecayPhotons)
    : DressedLeptons(barefs, barefs, dRmax, cut, useDecayPhotons)
  {     }




  int DressedLeptons::compare(const Projection& p) const {
    // Compare the two as final states (for pT and eta cuts)
    const DressedLeptons& other = dynamic_cast<const DressedLeptons&>(p);
    int fscmp = FinalState::compare(other);
    if (fscmp != EQUIVALENT) return fscmp;

    const PCmp phcmp = mkNamedPCmp(p, "Photons");
    if (phcmp != EQUIVALENT) return phcmp;

    const PCmp sigcmp = mkNamedPCmp(p, "Leptons");
    if (sigcmp != EQUIVALENT) return sigcmp;

    return (cmp(_dRmax, other._dRmax) ||
            cmp(_fromDecay, other._fromDecay));
  }


  void DressedLeptons::project(const Event& e) {
    _theParticles.clear();

    // Get bare leptons
    const FinalState& signal = applyProjection<FinalState>(e, "Leptons");
    Particles bareleptons = signal.particles();
    if (bareleptons.empty()) return;

    // Initialise DL collection with bare leptons
    vector<DressedLepton> allClusteredLeptons;
    for (const Particle& bl : bareleptons)
      allClusteredLeptons += DressedLepton(bl);

    // If the radius is 0 or negative, don't even attempt to cluster
    if (_dRmax > 0) {
      // Match each photon to its closest charged lepton within the dR cone
      const FinalState& photons = applyProjection<FinalState>(e, "Photons");
      for (const Particle& photon : photons.particles()) {
        // Ignore photon if it's from a hadron/tau decay and we're avoiding those
        if (!_fromDecay && photon.fromDecay()) continue;
        const FourMomentum& p_P = photon.momentum();
        double dRmin = _dRmax;
        int idx = -1;
        for (size_t i = 0; i < bareleptons.size(); ++i) {
          // Only cluster photons around *charged* signal particles
          if (bareleptons[i].charge3() == 0) continue;
          // Find the closest lepton
          const FourMomentum& p_l = bareleptons[i].momentum();
          double dR = deltaR(p_l, p_P);
          if (dR < dRmin) {
            dRmin = dR;
            idx = i;
          }
        }
        if (idx > -1) {
          allClusteredLeptons[idx].addPhoton(photon);
        }
      }
    }

    // Fill the canonical particles collection with the composite DL Particles
    for (const DressedLepton& lepton : allClusteredLeptons) {
      if (accept(lepton)) _theParticles.push_back(lepton);
    }

  }


}
