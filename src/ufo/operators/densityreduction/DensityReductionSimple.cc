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
  nsamples_ = 0;
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (jloc % reductionFactor == 0) {
      nsamples_++;
    }
  }

  // Fill latsModified, lonsModified and timesModified after resizing them
  // to the correct number of locations.
  latsModified.resize(nsamples_);
  lonsModified.resize(nsamples_);
  timesModified.resize(nsamples_);
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
}

void DensityReductionSimple::fillReducedGeoVaLs(GeoVaLs & geovals) const {
  const oops::Variables & gvars = geovals.getVars();

  // Lengths of each GeoVaL.
  std::vector<size_t> nLevelsGeoVaLs(gvars.size());
  for (size_t i = 0; i < gvars.size(); ++i) {
    nLevelsGeoVaLs[i] = geovals.nlevs(gvars[i], GeoVaLFormat::SAMPLED);
  }

  // Retrieve pathsGroupedByLocation.
  std::vector<util::Range<size_t>> pathsGroupedByLocation;
  geovals.getProfileIndicesGroupedByLocation(gvars[0], pathsGroupedByLocation);
  const size_t nlocs = pathsGroupedByLocation.size();

  // Fill in reduced vectors of GeoVaLs.
  std::vector<double> sampledValues(nsamples_);
  std::vector<double> reducedValues(nlocs);
  const int reductionFactor = params_.reductionFactor;
  for (size_t i = 0; i < gvars.size(); ++i) {
    const auto var = gvars[i];
    for (size_t jlev = 0; jlev < nLevelsGeoVaLs[i]; ++jlev) {
      geovals.getAtLevel(sampledValues, var, jlev, GeoVaLFormat::SAMPLED);
      for (size_t jloc = 0; jloc < nlocs; ++jloc) {
        reducedValues[jloc] = sampledValues[pathsGroupedByLocation[jloc].begin];
      }
      geovals.putAtLevel(reducedValues, var, jlev, GeoVaLFormat::REDUCED);
    }
  }
}

}  // namespace ufo
