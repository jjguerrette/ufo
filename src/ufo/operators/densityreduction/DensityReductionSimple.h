/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_DENSITYREDUCTION_DENSITYREDUCTIONSIMPLE_H_
#define UFO_OPERATORS_DENSITYREDUCTION_DENSITYREDUCTIONSIMPLE_H_

#include <vector>

#include "oops/util/parameters/NumericConstraints.h"
#include "ufo/operators/densityreduction/DensityReductionBase.h"

namespace ufo {

/// Parameters associated with the DensityReductionSimple class.
class DensitySimpleParameters : public DensityReductionParametersBase {
  OOPS_CONCRETE_PARAMETERS(DensitySimpleParameters, DensityReductionParametersBase)
 public:
  oops::Parameter<int> reductionFactor{"reduction factor",
      "length of blocks of GeoVaLs for which values at one location "
      "are copied to the other locations",
      1,
      this,
      {oops::minConstraint(1)}};
};

/// \brief Simple GeoVaL density reduction algorithm.
///
/// \details Given a `reduction factor` k, divide the observation locations into consecutive blocks
/// of length k. Retrieve the GeoVaL associated with the first location and copy it to the
/// subsequent k - 1 locations.
/// If k is equal to 1, there is no change relative to an observation operator that does
/// not use a modified `locations` method.
/// There is no guarantee that the observation locations are spatially close, so this
/// algorithm achieves a speed-up in GetValues at the potential cost to accuracy in H(x).
class DensityReductionSimple : public DensityReductionBase {
 public:
  typedef DensitySimpleParameters Parameters_;

  explicit DensityReductionSimple(const Parameters_ &,
                                  const ioda::ObsSpace &);
  ~DensityReductionSimple() {}

  void fillModifiedLocations(std::vector<float> &,
                             std::vector<float> &,
                             std::vector<util::DateTime> &,
                             std::vector<util::Range<size_t>> &) const override;

 private:
  const Parameters_ params_;
  mutable int nsamples_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_DENSITYREDUCTION_DENSITYREDUCTIONSIMPLE_H_
