// -*- C++ -*-
#ifndef RIVET_Sphericity_HH
#define RIVET_Sphericity_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/AxesDefinition.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Event.hh"
#include "Rivet/RivetCLHEP.hh"

namespace Rivet {

  /**
     @brief Calculate the sphericity event shape.
     
     The sphericity tensor (or quadratic momentum tensor) is defined as 
     \f[ 
     S^{\alpha \beta} = \frac{\sum_i p_i^\alpha p_i^\beta}{\sum_i |\mathbf{p}_i|^2} 
     \f],
     where the Greek indices are spatial components and the Latin indices are used
     for sums over particles. From this, the sphericity, aplanarity and planarity can be
     calculated by combinations of eigenvalues.
     
     Defining the three eigenvalues
     \f$ \lambda_1 \ge \lambda_2 \ge \lambda_3 \f$, with \f$ \lambda_1 + \lambda_2 + \lambda_3 = 1 \f$, 
     the sphericity is
     \f[ 
     S = \frac{3}{2} (\lambda_2 + \lambda_3) 
     \f]
     
     The aplanarity is \f$ A = \frac{3}{2}\lambda_3 \f$ and the planarity
     is \f$ P = \frac{2}{3}(S-2A) = \lambda_2 - \lambda_3 \f$. The eigenvectors define a 
     set of spatial axes comparable with the thrust axes, but more sensitive to 
     high momentum particles due to the quadratic sensitivity of the tensor to
     the particle momenta.
     
     Since the sphericity is quadratic in the particle momenta, it is not an
     infrared safe observable in perturbative QCD. This can be fixed by adding
     a regularizing power of \f$r\f$ to the definition:
     \f[ 
     S^{\alpha \beta} = 
     \frac{\sum_i |\mathbf{p}_i|^{r-2} p_i^\alpha p_i^\beta}
     {\sum_i |\mathbf{p}_i|^r} 
     \f]
     
     \f$r\f$ is available as a constructor argument on this class and will be
     taken into account by the Cmp<Projection> operation, so a single analysis
     can use several sphericity projections with different \f$r\f$ values without
     fear of a clash.
  */
  class Sphericity : public AxesDefinition {

  public:

    /// Constructor. Supplied FinalState projection must live throughout the run.
    Sphericity(FinalState& fsp, double rparam=2.0)
      : _sphericity(-1), _planarity(-1), _aplanarity(-1), _regparam(rparam), 
        _fsproj(fsp)
    { 
      addProjection(_fsproj);

      /// @todo Do this properly! Currently no calculation performed, so this result is meaningless.
      _sphAxes.push_back( Vector3(0,0,1) );
      _sphAxes.push_back( Vector3(0,0,1) );
      _sphAxes.push_back( Vector3(0,0,1) );
    }

    /// Return the name of the projection
    string getName() const {
      return "Sphericity";
    }

  protected:

    /// Perform the projection on the Event
    void project(const Event& e);

    /// Compare with other projections
    int compare(const Projection& p) const;

  public:

    /// @name Access the event shapes by name
    /// @{
    /// Sphericity
    const double sphericity() const { return _sphericity; }
    /// Planarity
    const double planarity() const { return _planarity; }
    /// Aplanarity
    const double aplanarity() const { return _aplanarity; }
    /// @}

    /// @name Access the sphericity basis vectors
    /// @{
    /// Sphericity axis
    const Vector3& sphericityAxis() const { return _sphAxes[0]; }
    /// Sphericity major axis
    const Vector3& sphericityMajorAxis() const { return _sphAxes[1]; }
    /// Sphericity minor axis
    const Vector3& sphericityMinorAxis() const { return _sphAxes[2]; }
    /// @}

    ///@{ AxesDefinition axis accessors.
    const Vector3& axis1() const { return sphericityAxis(); }
    const Vector3& axis2() const { return sphericityMajorAxis(); }
    const Vector3& axis3() const { return sphericityMinorAxis(); }
    ///@}


    /// @name Access the momentum tensor eigenvalues
    /// @{
    inline const double lambda1() const { return _lambdas[0]; }
    inline const double lambda2() const { return _lambdas[1]; }
    inline const double lambda3() const { return _lambdas[2]; }
    /// @}


  private:
    /// Eigenvalues.
    double _lambdas[3];

    /// The event shape scalars.
    double _sphericity, _planarity, _aplanarity;

    /// Sphericity axes.
    vector<Vector3> _sphAxes;

    /// Regularizing parameter, used to force infra-red safety.
    const double _regparam;

    /// The FinalState projection used by this projection.
    FinalState _fsproj;

  };

}


#endif
