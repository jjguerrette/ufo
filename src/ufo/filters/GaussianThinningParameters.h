/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_GAUSSIANTHINNINGPARAMETERS_H_
#define UFO_FILTERS_GAUSSIANTHINNINGPARAMETERS_H_

#include <string>
#include <utility>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/FilterParametersBase.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

enum class DistanceNorm {
  GEODESIC, MAXIMUM
};

struct DistanceNormParameterTraitsHelper {
  typedef DistanceNorm EnumType;
  static constexpr char enumTypeName[] = "DistanceNorm";
  static constexpr util::NamedEnumerator<DistanceNorm> namedValues[] = {
    { DistanceNorm::GEODESIC, "geodesic" },
    { DistanceNorm::MAXIMUM, "maximum" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::DistanceNorm> :
    public EnumParameterTraits<ufo::DistanceNormParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// \brief Options controlling the operation of the Gaussian_Thinning filter.
class GaussianThinningParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(GaussianThinningParameters, FilterParametersBase)

 public:
  /// Reimplemented to detect incompatible options.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  // Horizontal grid

  /// Cell size (in km) along the meridians. Thinning in the horizontal direction is disabled if
  /// this parameter is negative.
  // FIXME(wsmigaj): for consistency with vertical thinning, I think it would be better to
  // interpret absence of the horizontal_mesh setting as "don't do horizontal thinning"
  // rather than using an arbitrary default grid size. Leaving things as they are for backward
  // compatibility, for now.
  oops::Parameter<float> horizontalMesh{"horizontal_mesh", defaultHorizontalMesh(), this};
  /// True to use a reduced grid, with fewer cells at high latitudes.
  /// False to use a regular grid, with the same number of cells at all latitudes.
  oops::Parameter<bool> useReducedHorizontalGrid{"use_reduced_horizontal_grid", true, this};
  /// True to set the number of zonal bands so that the band width is as close as possible to
  /// \c horizontal_mesh, and the number of bins in each zonal band so that the bin width in the
  /// zonal direction is as close as possible to that in the meridional direction.
  /// False to set the number of zonal bands so that the band width is as small as possible, but
  /// no smaller than \c horizontal_mesh, and the bin width in the zonal direction is as small as
  /// possible, but no smaller than in the meridional direction.
  ///
  /// Defaults to \c false unless the \c ops_compatibility_mode option is enabled, in which case
  /// it's set to \c true.
  oops::OptionalParameter<bool> roundHorizontalBinCountToNearest{
    "round_horizontal_bin_count_to_nearest", this};
  /// Set this option to \c true to calculate partioning of longitude bins explicitly using
  /// horizontal mesh distance. By default, with this option set to \c false, calculating the number
  /// of longitude bins per latitude bin index involves the integer number of latitude
  /// bins. Setting this option to \c true adopts the Met Office OPS method whereby the
  /// integer number of latitude bins is replaced, in the calculation of longitude bins, by the
  /// Earth half-circumference divided by the horizontal mesh distance.
  oops::Parameter<bool> partitionLongitudeBinsUsingMesh{
    "partition_longitude_bins_using_mesh", false, this};

  // Vertical grid

  /// Cell size in the vertical direction. Thinning in the vertical direction is disabled
  /// if this parameter is not specified or negative.
  oops::Parameter<float> verticalMesh{"vertical_mesh", -1.0f, this};
  /// Lower bound of the vertical coordinate interval split into cells of size \c vertical_mesh.
  oops::Parameter<float> verticalMin{"vertical_min", 100.0f, this};
  /// Upper bound of the vertical coordinate interval split into cells of size \c vertical_mesh.
  /// This parameter is rounded upwards to the nearest multiple of \c vertical_mesh starting from
  /// \c vertical_min.
  oops::Parameter<float> verticalMax{"vertical_max", 110000.0f, this};
  /// Observation vertical coordinate.
  oops::Parameter<std::string> verticalCoord{"vertical_coordinate", "air_pressure", this};

  // Temporal grid

  /// Cell size in the temporal direction. Temporal thinning is disabled if this this parameter is
  /// not specified or set to 0.
  oops::OptionalParameter<util::Duration> timeMesh{"time_mesh", this};
  /// Lower bound of the time interval split into cells of size \c time_mesh. Temporal thinning is
  /// disabled if this parameter is not specified.
  oops::OptionalParameter<util::DateTime> timeMin{"time_min", this};
  /// Upper bound of the time interval split into cells of size \c time_mesh.
  /// This parameter is rounded upwards to the nearest multiple of \c time_mesh starting from
  /// \c time_min. Temporal thinning is disabled if this parameter is not specified.
  oops::OptionalParameter<util::DateTime> timeMax{"time_max", this};

  // Observation categories

  /// A string-valued or integer-valued variable. Observations with different values of that
  /// variable are thinned separately.
  ///
  /// Note: The filter will automatically detect if the chosen variable was also used to group
  /// observations into records when the ObsSpace was constructed, and if so, avoid exchanging
  /// data with other MPI processes, since in these circumstances each process can thin its
  /// observations independently from others.
  ///
  /// The variable used to group observations into records can be set with the
  /// `obs space.obsdatain.obsgrouping.group variable` YAML option.
  oops::OptionalParameter<Variable> categoryVariable{"category_variable", this};

  // Selection of observations to retain

  /// Variable storing observation priorities. Among all observations in a cell, only those with
  /// the highest priority are considered as candidates for retaining. If not specified, all
  /// observations are assumed to have equal priority.
  oops::OptionalParameter<Variable> priorityVariable{"priority_variable", this};

  /// Determines which of the highest-priority observations lying in a cell is retained.
  ///
  /// Allowed values:
  ///
  /// - \c geodesic: retain the observation closest to the cell centre in the horizontal direction
  ///   (the vertical coordinate and time are ignored)
  /// - \c maximum: retain the observation lying furthest from the cell's bounding box in the
  ///   system of coordinates in which the cell is a unit cube (all dimensions along which thinning
  ///   is enabled are taken into account).
  ///
  /// Defaults to \c geodesic unless the \c ops_compatibility_mode option is enabled, in which case
  /// it's set to \c maximum.
  oops::OptionalParameter<DistanceNorm> distanceNorm{"distance_norm", this};

  /// Set this option to \c true to make the filter produce identical results as the Ops_Thinning
  /// subroutine from the Met Office OPS system when both are run serially (on a single process).
  ///
  /// The filter behavior is modified in several ways:
  ///
  /// - The \c round_horizontal_bin_count_to_nearest option is set to \c true.
  ///
  /// - The \c distance_norm option is set to \c maximum.
  ///
  /// - Bin indices are calculated by rounding values away from rather towards zero. This can alter
  ///   the bin indices assigned to observations lying at bin boundaries.
  ///
  /// - The bin lattice is assumed to cover the whole real axis (for times and pressures) or the
  ///   [-360, 720] degrees interval (for longitudes) rather than just the intervals [\c time_min,
  ///   \c time_max], [\c pressure_min, \c pressure_max] and [0, 360] degrees, respectively. This
  ///   may cause observations lying at the boundaries of the latter intervals to be put in bins of
  ///   their own, which is normally undesirable.
  ///
  /// - A different (non-stable) sorting algorithm is used to order observations before inspection.
  ///   This can alter the set of retained observations if some bins contain multiple equally good
  ///   observations (with the same priority and distance to the cell center measured with the
  ///   selected norm). If this happens for a significant fraction of bins, it may be a sign the
  ///   criteria used to rank observations (the priority and the distance norm) are not specific
  ///   enough.
  oops::Parameter<bool> opsCompatibilityMode{"ops_compatibility_mode", false, this};

  /// Option to choose how to treat observations where there are multiple filter variables. If true,
  /// treats an observation location as valid if any filter variables have not been rejected.
  /// If false, observations are treated as valid only if all filter variables have passed QC.
  /// This is an optional parameter, if omitted the default value is true.
  oops::Parameter<bool>
    thinIfAnyFilterVariablesAreValid{"thin_if_any_filter_variables_are_valid", true, this};

 private:
  static float defaultHorizontalMesh() {
    return static_cast<float>(2 * M_PI * Constants::mean_earth_rad / 360.0);
  }
};

}  // namespace ufo

#endif  // UFO_FILTERS_GAUSSIANTHINNINGPARAMETERS_H_
