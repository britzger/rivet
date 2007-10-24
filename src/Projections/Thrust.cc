// -*- C++ -*-

#include "Rivet/Projections/Thrust.hh"
#include "Rivet/RivetCLHEP.hh"


namespace Rivet {


  void Thrust::calcT(const vector<Vector3>& p, double& t, Vector3& taxis) const {
    t = 0.0;
    Vector3 tv, ptot;
    vector<Vector3> cpm(4);
    const size_t psize = p.size();
    for (size_t k = 1; k < psize; ++k) {
      for (size_t j = 0; j < k; ++j) {
        tv = p[j].cross(p[k]);
        ptot = Vector3();
        for (size_t l = 0; l < psize; ++l) {
          if (l != j && l != k) {
            if (p[l].dot(tv) > 0.0) { 
              ptot += p[l];
            } else {
              ptot -= p[l];
            }
          }
        }
        cpm[0] = ptot - p[j] - p[k];
        cpm[1] = ptot - p[j] + p[k];
        cpm[2] = ptot + p[j] - p[k];
        cpm[3] = ptot + p[j] + p[k];
        for (size_t i = 0; i < 4; ++i) {
          double tval = cpm[i].mag2();
          if (tval > t) {
            t = tval;
            taxis = cpm[i];
          }
        }
      } // j loop
    } // k loop
  }



  void Thrust::calcM(const vector<Vector3>& p, double& m, Vector3& maxis) const {
    double mval;
    m = 0.0;
    Vector3 tv, ptot;
    vector<Vector3> cpm;
    for (unsigned int j = 0; j < p.size(); ++j) {
      tv = p[j];
      ptot = Vector3();
      for (unsigned int l = 0; l < p.size(); ++l) {
        if (l != j) {
          if (p[l].dot(tv) > 0.0) { 
            ptot += p[l];
          } else {
            ptot -= p[l];
          }
        }
      }
      cpm.clear();
      cpm.push_back(ptot - p[j]);
      cpm.push_back(ptot + p[j]);
      for (vector<Vector3>::iterator it = cpm.begin(); it != cpm.end(); ++it) {
        mval = it->mag2();
        if (mval > m) {
          m = mval;
          maxis = *it;
        }
      }
    } // j loop
  }



  void Thrust::calcThrust(const FinalState& fs) {
    // If we've already been here, stop
    //if (_calculatedThrust) return;

    // Make a vector of the three-momenta in the final state
    vector<Vector3> threeMomenta;
    double momentumSum(0.0);
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      Vector3 p3 = p->getMomentum().vect();
      threeMomenta.push_back(p3);
      momentumSum += threeMomenta.back().mag();
    }

    // Clear the caches
    _thrusts.clear();
    _thrustAxes.clear(); 
    //_calculatedThrust = true; // Pre-emptive strike to avoid repeats in each special case


    // If there are fewer than 2 visible particles, we can't do much
    if (threeMomenta.size() < 2) {
      for (int i = 0; i < 3; ++i) {
        _thrusts.push_back(-1);
        _thrustAxes.push_back(Vector3());
      }
      return;
    }


    // Handle special case of thrust = 1 if there are only 2 particles
    if (threeMomenta.size() == 2) {
      Vector3 axis;
      _thrusts.push_back(1.0);
      _thrusts.push_back(0.0);
      _thrusts.push_back(0.0);
      axis = threeMomenta[0].unit();
      if (axis.z() < 0) axis = -axis;
      _thrustAxes.push_back(axis);
      _thrustAxes.push_back(axis.orthogonal());
      _thrustAxes.push_back( _thrustAxes[0].cross(_thrustAxes[1]) );
      return;
    }


    // Handle special case where thrust = max{x1,x2,x3} if there are 3 particles
    if (threeMomenta.size() == 3) {
      Vector3 axis;
      // Order by magnitude
      /// @todo rsort with mag2() functor is much neater
      if (threeMomenta[0].mag2() < threeMomenta[1].mag2()) std::swap(threeMomenta[0], threeMomenta[1]);
      if (threeMomenta[0].mag2() < threeMomenta[2].mag2()) std::swap(threeMomenta[0], threeMomenta[2]);
      if (threeMomenta[1].mag2() < threeMomenta[2].mag2()) std::swap(threeMomenta[1], threeMomenta[2]);
      // Thrust
      axis = threeMomenta[0].unit();
      if (axis.z() < 0) axis = -axis;
      _thrusts.push_back(2.0 * threeMomenta[0].mag() / momentumSum);
      _thrustAxes.push_back(axis);
      // Thrust major (independent part of next largest momentum)
      axis = ( threeMomenta[1] - (axis.dot(threeMomenta[1])) * axis ).unit();
      if (axis.x() < 0) axis = -axis;
      _thrusts.push_back( (fabs(threeMomenta[1].dot(axis)) + 
                           fabs(threeMomenta[2].dot(axis)))  / momentumSum);
      _thrustAxes.push_back(axis);
      // Thrust minor
      _thrusts.push_back(0.0);
      _thrustAxes.push_back( _thrustAxes[0].cross(_thrustAxes[1]) );
      return;
    }


    // If the special cases don't apply, we have to use a general method. Here
    // we explicitly calculate in units of MeV using an algorithm based on 
    // Brandt/Dahmen Z Phys C1 (1978) and 'tasso' code from HERWIG. Re-coded
    // from Herwig++ implementation by Stefan Gieseke by Andy Buckley.
    // NB. special case with >= 4 coplanar particles will still fail. 
    // probably not too important... 
    /// @todo Thrust assumes all momenta are in the CoM system: no explicit boost is performed.

    // Temporary variables for calcs
    Vector3 axis; double val;

    // Get thrust
    calcT(threeMomenta, val, axis);
    _thrusts.push_back(sqrt(val) / momentumSum);
    if (axis.z() < 0) axis = -axis;
    _thrustAxes.push_back(axis.unit()); 

    // Get thrust major 
    threeMomenta.clear();
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      // Get the part of each 3-momentum which is perpendicular to the thrust axis
      const Vector3 v = p->getMomentum().vect();
      const Vector3 vpar = v.dot(axis.unit()) * axis.unit();
      threeMomenta.push_back(v - vpar);
    }
    calcM(threeMomenta, val, axis);
    _thrusts.push_back(sqrt(val) / momentumSum);
    if (axis.x() < 0) axis = -axis;
    _thrustAxes.push_back(axis.unit()); 

    // Get thrust minor
    if (_thrustAxes[0].dot(_thrustAxes[1]) < 1e-10) {
      axis = _thrustAxes[0].cross(_thrustAxes[1]);
      _thrustAxes.push_back(axis);
      val = 0.0;
      for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
        val += fabs(axis * p->getMomentum().vect());
      }
      _thrusts.push_back(val / momentumSum);
    } else {
      _thrusts.push_back(-1.0);
      _thrustAxes.push_back(Vector3());
    }

  }


  void Thrust::project(const Event& e) {
    const FinalState& fs = e.applyProjection(_fsproj);
    calcThrust(fs);
  }


}
