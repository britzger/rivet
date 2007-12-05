// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/ZEUS_2001_S4815815.hh"
using namespace Rivet;

#include "Rivet/RivetAIDA.hh"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////


// Book histograms
void ZEUS_2001_S4815815::init() {
  /// @todo This doesn't seem to correspond to the plots in the paper (SPIRES 4730372)
  _histJetEt1 = bookHistogram1D("JetET1", "Jet transverse energy", 11, 14.0, 75.0);
}


// Do the analysis
void ZEUS_2001_S4815815::analyze(const Event& event) {
  Log& log = getLog();
  log << Log::DEBUG << "Starting analyzing" << endl;

  // Analyse and print some info
  const FastJets& jets = event.applyProjection(_jetsproj);
  const size_t nj = jets.getNJets();
  log << Log::INFO << "Jet multiplicity = " << nj << endl;

  // Fill histograms
  typedef vector<fastjet::PseudoJet> Jets;
  Jets jetList = jets.getJets(); // was getJetsEt()
  for (Jets::const_iterator j = jetList.begin(); j != jetList.end(); ++j) {
    _histJetEt1->fill(j->perp(), event.weight() );
  }
  
  // Finished
  log << Log::DEBUG << "Finished analyzing" << endl;
}


// Finalize
void ZEUS_2001_S4815815::finalize() { }
