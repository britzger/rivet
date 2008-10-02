#ifndef RIVET_PARTICLENAME_HH 
#define RIVET_PARTICLENAME_HH

#include "Rivet/Rivet.hh"

namespace Rivet {

  /// Enumeration of available beam particles (using PDG IDs where available)
  enum ParticleName { 
    ELECTRON = 11, 
    POSITRON = -11, 
    PROTON = 2212, 
    ANTIPROTON = -2212, 
    PHOTON = 22, 
    NEUTRON = 2112, 
    ANTINEUTRON = 2112, 
    MUON = 13, 
    ANTIMUON = -13,
    NU_E = 12, 
    NU_EBAR = -12,
    NU_MU = 14, 
    NU_MUBAR = -14,
    NU_TAU = 16, 
    NU_TAUBAR = -16, 
    PIPLUS = 211, 
    PIMINUS = -211,
    TAU = 15, 
    ANTITAU = -15,
    EMINUS = 11, 
    EPLUS = -11, 
    P = 2212, 
    PBAR = -2212,
    GLUON = 21,
    GAMMA = 22,
    WPLUSBOSON = 24,
    WMINUSBOSON = -24,
    ZBOSON = 23,
    HIGGS = 25,
    ANY = 10000,
    PHOTOELECTRON,
    PHOTOPOSITRON,
    PHOTOMUON,     
    PHOTOANTIMUON,
    PHOTOTAU,      
    PHOTOANTITAU
  };


  /// Typedef for a map of beam particle name enums to strings.
  typedef std::map<ParticleName, std::string> ParticleNameMap;


  /// Typedef for a map of beam particle name strings to enums.
  typedef std::map<std::string, ParticleName> ParticleNameMapR;


  /// Function which returns a map from beam particle enums to the corresponding name strings.
  inline ParticleNameMap getParticleNamesMap() {
    ParticleNameMap bpmap;
    bpmap[ELECTRON] = "ELECTRON";
    bpmap[POSITRON] = "POSITRON";
    bpmap[PROTON] = "PROTON";
    bpmap[ANTIPROTON] = "ANTIPROTON";
    bpmap[PHOTON] = "PHOTON";
    bpmap[NEUTRON] = "NEUTRON";
    bpmap[ANTINEUTRON] = "ANTINEUTRON";
    bpmap[MUON] = "MUON";
    bpmap[ANTIMUON] = "ANTIMUON";
    bpmap[NU_E] = "NU_E";
    bpmap[NU_EBAR] = "NU_EBAR";
    bpmap[NU_MU] = "NU_MU";
    bpmap[NU_MUBAR] = "NU_MUBAR";
    bpmap[NU_TAU] = "NU_TAU";
    bpmap[NU_TAUBAR] = "NU_TAUBAR";
    bpmap[PIPLUS] = "PIPLUS"; 
    bpmap[PIMINUS] = "PIMINUS";
    bpmap[TAU] = "TAU"; 
    //bpmap[WPLUSBOSON] = "WPLUSBOSON";
    //bpmap[WMINUSBOSON] = "WMINUSBOSON";
    //bpmap[ZBOSON] = "ZBOSON";
    //bpmap[HIGGS] = "HIGGS";
    bpmap[ANTITAU] = "ANTITAU";
    bpmap[PHOTOELECTRON] = "PHOTOELECTRON";
    bpmap[PHOTOPOSITRON] = "PHOTOPOSITRON";
    bpmap[PHOTOMUON] = "PHOTOMUON";
    bpmap[PHOTOANTIMUON] = "PHOTOANTIMUON";
    bpmap[PHOTOTAU] = "PHOTOTAU"; 
    bpmap[PHOTOANTITAU] = "PHOTOANTITAU";
    bpmap[ANY] = "*";
    return bpmap;
  }

  /// Function which returns a map from beam particle name strings to the corresponding enums.
  inline ParticleNameMapR getParticleNamesRMap() {
    ParticleNameMap bpmap = getParticleNamesMap();
    ParticleNameMapR bpmapr;
    for (ParticleNameMap::const_iterator bp = bpmap.begin(); bp != bpmap.end(); ++bp) {
      bpmapr[bp->second] = bp->first;
    }
    return bpmapr;
  }


  /// Typedef for a collection of beam particle name enums.
  typedef std::vector<ParticleName> ParticleNameList;


  /// Function which returns a vector of all the beam particle values in 
  /// the ParticleName enum.
  inline ParticleNameList getParticleNameEnums() {
    ParticleNameList names;
    ParticleNameMap bpmap = getParticleNamesMap();
    for (ParticleNameMap::const_iterator bp = bpmap.begin(); bp != bpmap.end(); ++bp) {
      names.push_back(bp->first);
    }
    return names;
  }


  /// Function which returns a vector of all the beam particle values in 
  /// the ParticleName enum.
  inline ParticleName getParticleNameEnum(const std::string& pname) {
    return Rivet::getParticleNamesRMap()[pname];
  }



  /// Function which returns a vector of all the beam particle name strings.
  inline std::vector<std::string> getParticleNames() {
    vector<string> names;
    ParticleNameMap bpmap = getParticleNamesMap();
    for (ParticleNameMap::const_iterator bp = bpmap.begin(); bp != bpmap.end(); ++bp) {
      names.push_back(bp->second);
    }
    return names;
  }


  /// Print a ParticleName as a string.
  inline std::string toString(const ParticleName& p) {
    return getParticleNamesMap()[p];
  }

  /// Allow ParticleName to be passed to an iostream.
  inline std::ostream& operator<<(std::ostream& os, const ParticleName& p) {
    os << toString(p);
    return os;
  }

  /////////////////////////////////////////////////

  

  /// Typedef for a pair of beam particle names.
  typedef std::pair<ParticleName, ParticleName> BeamPair;


  /// Print a BeamPair as a string.
  inline std::string toString(const BeamPair& pair) {
    string out = "[" + toString(pair.first) + ", " + toString(pair.second) + "]";
    return out;
  }

  /// Allow BeamPair to be passed to an iostream.
  inline std::ostream& operator<<(std::ostream& os, const BeamPair& bp) {
    os << toString(bp);
    return os;
  }

  
}


#endif
