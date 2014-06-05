// -*- C++ -*-
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  /// @brief BaBar pion, kaon and proton production in the continuum
  /// @author Peter Richardson
  class BABAR_2013_I1238276 : public Analysis {
  public:

    BABAR_2013_I1238276()
      : Analysis("BABAR_2013_I1238276")
    { }


    void analyze(const Event& e) {
      const double weight = e.weight();

      // Loop through charged FS particles and look for charmed mesons/baryons
      const ChargedFinalState& fs = applyProjection<ChargedFinalState>(e, "FS");

      const Beam beamproj = applyProjection<Beam>(e, "Beams");
      const ParticlePair& beams = beamproj.beams();
      FourMomentum mom_tot = beams.first.momentum() + beams.second.momentum();
      LorentzTransform cms_boost(-mom_tot.boostVector());
      MSG_DEBUG("CMS Energy sqrt s = " << beamproj.sqrtS());

      foreach (const Particle& p, fs.particles()) {
	// check if prompt or not
	const GenParticle* pmother = p.genParticle();
	GenVertex* ivertex = pmother->production_vertex();
	bool prompt = true;
	while (ivertex) {
	  int n_inparts = ivertex->particles_in_size();
	  if(n_inparts<1) break;
	  const HepMC::GenVertex::particles_in_const_iterator iPart_invtx = ivertex->particles_in_const_begin();
	  pmother = (*iPart_invtx);                   // first mother particle
	  int mother_pid = abs(pmother->pdg_id());
	  if(mother_pid==PID::K0S || mother_pid==PID::LAMBDA) {
	    prompt = false;
	    break;
	  }
	  else if (mother_pid<6) {
	    break;
	  }
	  ivertex = pmother->production_vertex();
	}


        // momentum in CMS frame
        const double mom = cms_boost.transform(p.momentum()).vector3().mod();
        const int PdgId = p.abspid();
        MSG_DEBUG("pdgID = " << PdgId << " Momentum = " << mom);
        switch (PdgId) {
	case PID::PIPLUS:
	  if(prompt) _histPion_no_dec->fill(mom,weight);
	  _histPion_dec   ->fill(mom,weight);
	  break;
	case PID::KPLUS:
	  if(prompt) _histKaon_no_dec->fill(mom,weight);
	  _histKaon_dec   ->fill(mom,weight);
	  break;
	case PID::PROTON:
	  if(prompt) _histProton_no_dec->fill(mom,weight);
	  _histProton_dec   ->fill(mom,weight);
	default :
	  break;
        }
      }
    } // analyze


    void finalize() {

      scale(_histPion_no_dec  ,1./sumOfWeights());
      scale(_histKaon_no_dec  ,1./sumOfWeights());
      scale(_histProton_no_dec,1./sumOfWeights());
      scale(_histPion_dec     ,1./sumOfWeights());
      scale(_histKaon_dec     ,1./sumOfWeights());
      scale(_histProton_dec   ,1./sumOfWeights());
    } // finalize


    void init() {
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "FS");

      _histPion_no_dec   = bookHisto1D(1,1,1);
      _histKaon_no_dec   = bookHisto1D(1,1,2);
      _histProton_no_dec = bookHisto1D(1,1,3);
      _histPion_dec      = bookHisto1D(2,1,1);
      _histKaon_dec      = bookHisto1D(2,1,2);
      _histProton_dec    = bookHisto1D(2,1,3);

    } // init

  private:

    //@{
    // Histograms for continuum data (sqrt(s) = 10.52 GeV)
    // no K_S and Lambda decays
    Histo1DPtr _histPion_no_dec;
    Histo1DPtr _histKaon_no_dec;
    Histo1DPtr _histProton_no_dec;
    // including decays
    Histo1DPtr _histPion_dec;
    Histo1DPtr _histKaon_dec;
    Histo1DPtr _histProton_dec;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2013_I1238276);

}
