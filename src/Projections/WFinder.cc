// -*- C++ -*-
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/LeptonClusters.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  WFinder::WFinder(double etaMin, double etaMax,
                   double pTmin,
                   PdgId pid,
                   double minmass, double maxmass,
                   double missingET,
                   double dRmax, bool clusterPhotons, bool trackPhotons,
                   double masstarget,
                   bool useTransverseMass) {
    vector<pair<double, double> > etaRanges;
    etaRanges += std::make_pair(etaMin, etaMax);
    _init(etaRanges, pTmin, pid, minmass, maxmass, missingET,
          dRmax, clusterPhotons, trackPhotons, masstarget, useTransverseMass);
  }


  WFinder::WFinder(const std::vector<std::pair<double, double> >& etaRanges,
                   double pTmin,
                   PdgId pid,
                   double minmass, double maxmass,
                   double missingET,
                   double dRmax, bool clusterPhotons, bool trackPhotons,
                   double masstarget,
                   bool useTransverseMass) {
    _init(etaRanges, pTmin, pid, minmass, maxmass, missingET,
          dRmax, clusterPhotons, trackPhotons, masstarget, useTransverseMass);
  }


  void WFinder::_init(const std::vector<std::pair<double, double> >& etaRanges,
                      double pTmin,
                      PdgId pid,
                      double minmass, double maxmass,
                      double missingET,
                      double dRmax, bool clusterPhotons, bool trackPhotons,
                      double masstarget,
                      bool useTransverseMass)
  {
    setName("WFinder");

    _minmass = minmass;
    _maxmass = maxmass;
    _masstarget = masstarget;
    _pid = pid;
    _trackPhotons = trackPhotons;
    _useTransverseMass = useTransverseMass;

    // Check that the arguments are legal
    assert(abs(_pid) == ELECTRON || abs(_pid) == MUON);
    _nu_pid = abs(_pid) + 1;
    assert(abs(_nu_pid) == NU_E || abs(_nu_pid) == NU_MU);

    // Don't make pT or eta cuts on the neutrino
    IdentifiedFinalState neutrinos;
    neutrinos.acceptNeutrinos();
    addProjection(neutrinos, "Neutrinos");

    // Lepton clusters
    FinalState fs;
    IdentifiedFinalState bareleptons(fs);
    bareleptons.acceptIdPair(pid);
    LeptonClusters leptons(fs, bareleptons, dRmax,
                           clusterPhotons,
                           etaRanges, pTmin);
    addProjection(leptons, "LeptonClusters");

    // Add MissingMomentum proj to calc MET
    MissingMomentum vismom(fs);
    addProjection(vismom, "MissingET");
    // Set ETmiss
    _etMiss = missingET;

    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(*this);
    addProjection(remainingFS, "RFS");
  }


  /////////////////////////////////////////////////////


  const FinalState& WFinder::remainingFinalState() const {
    return getProjection<FinalState>("RFS");
  }

  int WFinder::compare(const Projection& p) const {
    PCmp LCcmp = mkNamedPCmp(p, "LeptonClusters");
    if (LCcmp != EQUIVALENT) return LCcmp;

    const WFinder& other = dynamic_cast<const WFinder&>(p);
    return (cmp(_minmass, other._minmass) || cmp(_maxmass, other._maxmass) ||
            cmp(_useTransverseMass, other._useTransverseMass) ||
            cmp(_etMiss, other._etMiss) ||
            cmp(_pid, other._pid) || cmp(_trackPhotons, other._trackPhotons));

  }


  void WFinder::project(const Event& e) {
    clear();

    const LeptonClusters& leptons = applyProjection<LeptonClusters>(e, "LeptonClusters");
    const FinalState& neutrinos = applyProjection<FinalState>(e, "Neutrinos");

    // Make and register an invariant mass final state for the W decay leptons
    vector<pair<PdgId, PdgId> > l_nu_ids;
    l_nu_ids += make_pair(abs(_pid), -abs(_nu_pid));
    l_nu_ids += make_pair(-abs(_pid), abs(_nu_pid));
    InvMassFinalState imfs(FinalState(), l_nu_ids, _minmass, _maxmass, _masstarget);
    imfs.useTransverseMass(_useTransverseMass);
    ParticleVector tmp;
    tmp.insert(tmp.end(), leptons.clusteredLeptons().begin(), leptons.clusteredLeptons().end());
    tmp.insert(tmp.end(), neutrinos.particles().begin(), neutrinos.particles().end());
    imfs.calc(tmp);

    if (imfs.particlePairs().size() < 1) return;

    ParticlePair Wconstituents(imfs.particlePairs()[0]);
    Particle p1(Wconstituents.first), p2(Wconstituents.second);

    if (PID::threeCharge(p1)==0) {
      _constituentLeptons += p2;
      _constituentNeutrinos += p1;
    }
    else {
      _constituentLeptons += p1;
      _constituentNeutrinos += p2;
    }

    FourMomentum pW = p1.momentum() + p2.momentum();
    const int w3charge = PID::threeCharge(p1) + PID::threeCharge(p2);
    assert(abs(w3charge) == 3);
    const int wcharge = w3charge/3;

    stringstream msg;
    string wsign = (wcharge == 1) ? "+" : "-";
    string wstr = "W" + wsign;
    msg << wstr << " reconstructed from: " << endl
        << "   " << p1.momentum() << " " << p1.pdgId() << endl
        << " + " << p2.momentum() << " " << p2.pdgId() << endl;

    // Check missing ET
    const MissingMomentum& vismom = applyProjection<MissingMomentum>(e, "MissingET");
    /// @todo Restrict missing momentum eta range? Use vectorET()?
    if (vismom.scalarET() < _etMiss) {
      getLog() << Log::DEBUG << "Not enough missing ET: " << vismom.scalarET()/GeV
               << " GeV vs. " << _etMiss/GeV << " GeV" << endl;
      return;
    }

    // Make W Particle and insert into particles list
    const PdgId wpid = (wcharge == 1) ? WPLUSBOSON : WMINUSBOSON;
    _bosons.push_back(Particle(wpid, pW));

    // find the LeptonClusters and neutrinos which survived the IMFS cut such that we can
    // extract their original particles
    foreach (const Particle& p, _constituentNeutrinos) {
      _theParticles.push_back(p);
    }
    foreach (const Particle& p, _constituentLeptons) {
      foreach (const ClusteredLepton& l, leptons.clusteredLeptons()) {
        if (p.pdgId()==l.pdgId() && p.momentum()==l.momentum()) {
          _theParticles.push_back(l.constituentLepton());
          if (_trackPhotons) {
            _theParticles.insert(_theParticles.end(),
                                 l.constituentPhotons().begin(),
                                 l.constituentPhotons().end());
          }
        }
      }
    }

//    /// Apply the RFS here so that it can be acquired later via the remainingFinalState method
//    applyProjection<FinalState>(e, "RFS");
  }


}
