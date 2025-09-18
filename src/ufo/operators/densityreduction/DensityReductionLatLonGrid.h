/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_DENSITYREDUCTION_DENSITYREDUCTIONLATLONGRID_H_
#define UFO_OPERATORS_DENSITYREDUCTION_DENSITYREDUCTIONLATLONGRID_H_

#include <vector>

#include "oops/util/parameters/NumericConstraints.h"
#include "ufo/operators/densityreduction/DensityReductionBase.h"

namespace ufo {

/// Parameters associated with the DensityReductionLatLonGrid class.
class DensityLatLonGridParameters : public DensityReductionParametersBase {
  OOPS_CONCRETE_PARAMETERS(DensityLatLonGridParameters, DensityReductionParametersBase)
 public:
  oops::RequiredParameter<float> latMin{"latitude lower bound", this};
  oops::RequiredParameter<float> latMax{"latitude upper bound", this};
  oops::RequiredParameter<float> lonMin{"longitude lower bound", this};
  oops::RequiredParameter<float> dLat{"latitude grid length", this,
                                      {oops::exclusiveMinConstraint(0.0f)}};
  oops::RequiredParameter<float> dLon{"longitude grid length", this,
                                      {oops::exclusiveMinConstraint(0.0f)}};
};

/// \brief Lat-lon grid GeoVaL density reduction algorithm.
///
/// \details Divide the model domain into a regular lat-lon grid whose
/// dimensions and resolution are specified by the user.
/// One observation location is chosen in each grid box and passed to GetValues
/// in order to draw the GeoVaLs.
/// Any additional locations in that grid box are assigned copies of the GeoVaLs
/// that were produced for the selected location.
class DensityReductionLatLonGrid : public DensityReductionBase {
 public:
  typedef DensityLatLonGridParameters Parameters_;

  explicit DensityReductionLatLonGrid(const Parameters_ &,
                                      const ioda::ObsSpace &);
  ~DensityReductionLatLonGrid() {}

  void fillModifiedLocations(std::vector<float> &,
                             std::vector<float> &,
                             std::vector<util::DateTime> &,
                             std::vector<util::Range<size_t>> &) const override;

 private:
  const Parameters_ params_;

  // Number of modified lat/lon/time samples.
  mutable int numSamples_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_DENSITYREDUCTION_DENSITYREDUCTIONLATLONGRID_H_
