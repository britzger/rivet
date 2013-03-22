// -*- C++ -*-
#include "Rivet/Projections/LeptonClusters.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  LeptonClusters::LeptonClusters(const FinalState& photons, const FinalState& signal,
                 double dRmax, bool cluster,
                 const std::vector<std::pair<double, double> >& etaRanges,
                 double pTmin) :
    FinalState(etaRanges, pTmin),
    _dRmax(dRmax), _cluster(cluster)
  {
    setName("LeptonClusters");

    IdentifiedFinalState photonfs(photons);
    photonfs.acceptId(PID::PHOTON);
    addProjection(photonfs, "Photons");
    addProjection(signal, "Signal");
  }


  int LeptonClusters::compare(const Projection& p) const {
    // Compare the two as final states (for pT and eta cuts)
    const LeptonClusters& other = dynamic_cast<const LeptonClusters&>(p);
    int fscmp = FinalState::compare(other);
    if (fscmp != EQUIVALENT) return fscmp;

    const PCmp phcmp = mkNamedPCmp(p, "Photons");
    if (phcmp != EQUIVALENT) return phcmp;

    const PCmp sigcmp = mkNamedPCmp(p, "Signal");
    if (sigcmp != EQUIVALENT) return sigcmp;

    return (cmp(_dRmax, other._dRmax) || cmp(_cluster, other._cluster));
  }


  void LeptonClusters::project(const Event& e) {
    _theParticles.clear();
    _clusteredLeptons.clear();

    const FinalState& signal = applyProjection<FinalState>(e, "Signal");
    Particles bareleptons = signal.particles();
    if (bareleptons.empty()) return;

    vector<ClusteredLepton> allClusteredLeptons;
    for (size_t i = 0; i < bareleptons.size(); ++i) {
      allClusteredLeptons.push_back(ClusteredLepton(bareleptons[i]));
    }

    const FinalState& photons = applyProjection<FinalState>(e, "Photons");
    foreach (const Particle& photon, photons.particles()) {
      const FourMomentum p_P = photon.momentum();
      double dRmin=_dRmax;
      int idx = -1;
      for (size_t i = 0; i < bareleptons.size(); ++i) {
        FourMomentum p_l = bareleptons[i].momentum();
        // Only cluster photons around *charged* signal particles
        if (PID::threeCharge(bareleptons[i].pdgId()) == 0) continue;
        // Geometrically match momentum vectors
        double dR = deltaR(p_l, p_P);
        if (dR < dRmin) {
          dRmin = dR;
          idx = i;
        }
      }
      if (idx > -1) {
        if (_cluster) allClusteredLeptons[idx].addPhoton(photon, _cluster);
      }
    }

    foreach (const ClusteredLepton& lepton, allClusteredLeptons) {
      if (accept(lepton)) {
        _clusteredLeptons.push_back(lepton);
        _theParticles.push_back(lepton.constituentLepton());
        _theParticles.insert(_theParticles.end(),
                             lepton.constituentPhotons().begin(),
                             lepton.constituentPhotons().end());
      }
    }
  }
}
