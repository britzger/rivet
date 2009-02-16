// -*- C++ -*-
#ifndef RIVET_H1_1995_S3167097_HH
#define RIVET_H1_1995_S3167097_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalStateHCM.hh"
#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {


  /// @brief Measures energy flow in DIS? To be checked!
  /// @todo Check this analysis!
  /// @author Leif Lonnblad
  class H1_1995_S3167097 : public Analysis {
  public:

    /// Constructor.
    H1_1995_S3167097() { 
      setBeams(ELECTRON, PROTON);
      const DISLepton& lepton = addProjection(DISLepton(ELECTRON, POSITRON), "Lepton");
      const DISKinematics& diskin = addProjection(DISKinematics(lepton, PROTON), "Kinematics");
      const FinalStateHCM& fshcm = addProjection(FinalStateHCM(diskin), "FS");
      addProjection(CentralEtHCM(fshcm), "Y1HCM");
      //addCut("x", MORE_EQ, _xmin);
      //addCut("x", LESS_EQ, _xmax);
    }

  public:

    /// Factory method
    static Analysis* create() { return new H1_1995_S3167097(); }

    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "3167097";
    }
    /// A short description of the analysis.
    string summary() const {
      return "TODO";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "H1";
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "1995";
    }
    //@}
    

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

    
  private:

    /// Calculate the bin number from the DISKinematics projection.
    int _getbin(const DISKinematics& dk);


    /// Some integer constants used.
    static const size_t _nb = 24, _nbin = 9;
    
    /// Some double constants used.
    static const double _xmin, _xmax;

    /// Histograms for the \f$ E_T \f$ flows
    vector<AIDA::IHistogram1D*> _hEtFlow, _hEtFlowStat;

    /// Histograms for averages in different kinematical bins.
    AIDA::IHistogram1D *_hAvEt, *_hAvX, *_hAvQ2, *_hN;

    /// Helper vector;
    vector<double> _nev;

  };

}

#endif
