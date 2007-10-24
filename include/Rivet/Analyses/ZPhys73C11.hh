// -*- C++ -*-
#ifndef RIVET_ZPHYS73C11_HH
#define RIVET_ZPHYS73C11_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {


  /// Implementation of DELPHI event shape paper
  class ZPhys73C11 : public Analysis {

  public:

    /// Default constructor.
    inline ZPhys73C11()
      : _cnfsproj(), _cfsproj(_cnfsproj),
        _cspherproj(_cfsproj), _cnspherproj(_cnfsproj), 
        _cthrustproj(_cfsproj), _cnthrustproj(_cnfsproj), 
        _cparisiproj(_cfsproj), _cnparisiproj(_cnfsproj),
        _chemiproj(_cfsproj, _cthrustproj), 
        _cnhemiproj(_cfsproj, _cnthrustproj)
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(_cfsproj);
      addProjection(_cnfsproj);
      addProjection(_beamsproj);
      addProjection(_cspherproj);
      addProjection(_cnspherproj);
      addProjection(_cthrustproj);
      addProjection(_cnthrustproj);
      addProjection(_cparisiproj);
      addProjection(_cnparisiproj);
      addProjection(_chemiproj);
      addProjection(_cnhemiproj);
    }

    /// Factory method.
    static Analysis* create() { 
      return new ZPhys73C11(); 
    }

    /// Get the name of this analysis.
    inline string getName() const {
      return "ZPhys73C11";
    }

    virtual void init();

    virtual void analyze(const Event& event);

    virtual void finalize();


  private:

    /// Hide the assignment operator
    ZPhys73C11& operator=(const ZPhys73C11&);


  private:

    /// The final state projector.
    FinalState _cnfsproj;

    /// Charged final state projector.
    ChargedFinalState _cfsproj;

    /// Projection to get the beams.
    Beam _beamsproj;

    /// Sphericity projections.
    Sphericity _cspherproj, _cnspherproj;

    /// Thrust projections.
    Thrust _cthrustproj, _cnthrustproj;

    /// Parisi tensor (C and D params) projections.
    ParisiTensor _cparisiproj, _cnparisiproj;

    /// Projections to calculate event hemisphere masses and broadenings.
    Hemispheres _chemiproj, _cnhemiproj;


    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_histPtTInC,  *_histPtTInCN;
    AIDA::IHistogram1D *_histPtTOutC, *_histPtTOutCN;
    AIDA::IHistogram1D *_histPtSInC,  *_histPtSInCN;
    AIDA::IHistogram1D *_histPtSOutC, *_histPtSOutCN;

    AIDA::IHistogram1D *_histRapidityTC, *_histRapidityTCN;
    AIDA::IHistogram1D *_histRapiditySC, *_histRapiditySCN;

    AIDA::IHistogram1D *_histScaledMom, *_histLogScaledMom;

    AIDA::IHistogram1D *_histPtTOutVsXp, *_histPtTInVsXp;

    AIDA::IHistogram1D *_hist1MinusTC, *_hist1MinusTCN; 
    AIDA::IHistogram1D *_histTMajorC, *_histTMajorCN; 
    AIDA::IHistogram1D *_histTMinorC, *_histTMinorCN; 
    AIDA::IHistogram1D *_histOblatenessC, *_histOblatenessCN; 

    AIDA::IHistogram1D *_histSphericityC, *_histSphericityCN;
    AIDA::IHistogram1D *_histAplanarityC, *_histAplanarityCN;
    AIDA::IHistogram1D *_histPlanarityC, *_histPlanarityCN;

    AIDA::IHistogram1D *_histCParamC, *_histCParamCN;
    AIDA::IHistogram1D *_histDParamC, *_histDParamCN;

    AIDA::IHistogram1D *_histHemiMassDC, *_histHemiMassDCN;
    AIDA::IHistogram1D *_histHemiMassHC, *_histHemiMassHCN;
    AIDA::IHistogram1D *_histHemiMassLC, *_histHemiMassLCN;
               
    AIDA::IHistogram1D *_histHemiBroadWC, *_histHemiBroadWCN;
    AIDA::IHistogram1D *_histHemiBroadNC, *_histHemiBroadNCN;
    AIDA::IHistogram1D *_histHemiBroadTC, *_histHemiBroadTCN;
    AIDA::IHistogram1D *_histHemiBroadDC, *_histHemiBroadDCN;

    AIDA::IHistogram1D *_histDiffRate2DurhamC, *_histDiffRate2DurhamCN;
    AIDA::IHistogram1D *_histDiffRate2JadeC,   *_histDiffRate2JadeCN; 
    AIDA::IHistogram1D *_histDiffRate3DurhamC, *_histDiffRate3DurhamCN;
    AIDA::IHistogram1D *_histDiffRate3JadeC,   *_histDiffRate3JadeCN;
    AIDA::IHistogram1D *_histDiffRate4DurhamC, *_histDiffRate4DurhamCN;
    AIDA::IHistogram1D *_histDiffRate4JadeC,   *_histDiffRate4JadeCN;

    AIDA::IHistogram1D *_histEEC, *_histAEEC;
    //@}
  
  };

}

#endif
