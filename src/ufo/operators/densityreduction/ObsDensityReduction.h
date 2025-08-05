/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_DENSITYREDUCTION_OBSDENSITYREDUCTION_H_
#define UFO_OPERATORS_DENSITYREDUCTION_OBSDENSITYREDUCTION_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/base/Locations.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/densityreduction/DensityReductionBase.h"
#include "ufo/operators/densityreduction/ObsDensityReductionParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Reduce the density of GeoVaL locations prior to calling GetValues.
///
/// \details This operator can be used to reduce the number of locations passed
/// to the GetValues class, leading to a faster execution time at the cost of
/// potentially reduced accuracy in the subsequent H(x) calculations.
/// The `locations` and `computeReducedVars` functions are overridden in order to
/// achieve this.
///
/// The ObsDensityReduction operator wraps an underlying observation operator,
/// selected with the `operator` parameter.
/// The ObsDensityReduction operator passes to the underlying operator a set of GeoVaLs
/// of length equal to the number of ObsSpace locations, but requests from OOPS a smaller
/// number of GeoVaLs.
///
/// The GeoVaLs passed to the underlying operator are referred to as the 'reduced' GeoVaLs,
/// and those retrieved from OOPS are referred to as the 'sampled' GeoVaLs.
/// There are, in general, more reduced than sampled GeoVaLs, and the
/// `computeReducedVars` method actually increases the number of GeoVaLs present.
/// This is contrary to the behaviour in other operators for which there are more
/// sampled than reduced GeoVaLs, but we retain the nomenclature here for consistency.
///
/// The algorithm used is chosen with the `algorithm` parameter and is a subclass
/// of the `DensityReductionBase` class.
class ObsDensityReduction : public ObsOperatorBase,
                            private util::ObjectCounter<ObsDensityReduction> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsDensityReductionParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;
  static const std::string classname() {return "ufo::ObsDensityReduction";}

  ObsDensityReduction(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsDensityReduction() override;

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override;

  void computeReducedVars(const oops::Variables & reducedVars, GeoVaLs & geovals) const override;

  Locations_ locations() const override;

 private:
  void print(std::ostream &) const override;

 private:
  const ioda::ObsSpace& odb_;
  std::unique_ptr<ObsOperatorBase> operator_;
  std::unique_ptr<DensityReductionBase> algorithm_;
  oops::Variables requiredVars_;

  // Flag indicating that the locations() method has been called.
  // This is necessary to ensure that ufo ctests of this operator work correctly.
  mutable bool locationsCalled_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_DENSITYREDUCTION_OBSDENSITYREDUCTION_H_
