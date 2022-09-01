/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_ATMVERTINTERP_OBSATMVERTINTERPTLAD_H_
#define UFO_OPERATORS_ATMVERTINTERP_OBSATMVERTINTERPTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/atmvertinterp/ObsAtmVertInterpParameters.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// AtmVertInterp observation operator
class ObsAtmVertInterpTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsAtmVertInterpTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsAtmVertInterpTLAD";}

  typedef ObsAtmVertInterpParameters Parameters_;

  ObsAtmVertInterpTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAtmVertInterpTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  oops::Variables simulatedVars() const override {return operatorVars_;}

  int & toFortran() {return keyOperAtmVertInterp_;}
  const int & toFortran() const {return keyOperAtmVertInterp_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperAtmVertInterp_;
  oops::Variables varin_;
  oops::Variables operatorVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_ATMVERTINTERP_OBSATMVERTINTERPTLAD_H_
