/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/densityreduction/DensityReductionSimple.h"

namespace ufo {

static DensityReductionMaker<DensityReductionSimple> makerDensityReductionSimple_("simple");

DensityReductionSimple::DensityReductionSimple(const Parameters_ & params,
                                               const ioda::ObsSpace & odb)
  : DensityReductionBase(params, odb),
    params_(params)
{}

void DensityReductionSimple::fillModifiedLocations
(std::vector<float> & latsModified,
 std::vector<float> & lonsModified,
 std::vector<util::DateTime> & timesModified,
 std::vector<util::Range<size_t>> & pathsGroupedByLocation) const {
  const int nlocs = odb_.nlocs();

  // Get lat/lon/time from ObsSpace.
  std::vector<float> lons(nlocs);
  std::vector<float> lats(nlocs);
  std::vector<util::DateTime> times(nlocs);
  odb_.get_db("MetaData", "longitude", lons);
  odb_.get_db("MetaData", "latitude", lats);
  odb_.get_db("MetaData", "dateTime", times);

  // Fill reduced lat/lon/time.
  const int reductionFactor = params_.reductionFactor;
  // Count the number of sampled locations.
  int numSamples = 0;
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (jloc % reductionFactor == 0) {
      numSamples++;
    }
  }

  // Fill latsModified, lonsModified and timesModified after resizing them
  // to the correct number of locations.
  latsModified.resize(numSamples);
  lonsModified.resize(numSamples);
  timesModified.resize(numSamples);
  size_t jlocModified = 0;
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (jloc % reductionFactor == 0) {
      lonsModified[jlocModified] = lons[jloc];
      latsModified[jlocModified] = lats[jloc];
      timesModified[jlocModified] = times[jloc];
      jlocModified++;
    }
  }

  // Set up pathsGroupedByLocation.
  pathsGroupedByLocation.resize(nlocs);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    const size_t jlocModified = jloc / reductionFactor;
    pathsGroupedByLocation[jloc].begin = jlocModified;
    pathsGroupedByLocation[jloc].end = jlocModified + 1;
  }
  updateNumSamples(numSamples);
}

}  // namespace ufo
