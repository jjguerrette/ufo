/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/MaxColumnGeoVaLs.h"

#include <cmath>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<MaxColumnGeoVaLs> makerMaxGeoVaLs_("MaxColumnGeoVaLs");

// -----------------------------------------------------------------------------

MaxColumnGeoVaLs::MaxColumnGeoVaLs(const eckit::LocalConfiguration config)
  : invars_() {
  options_.deserialize(config);

  // Include required variables
  invars_ += Variable("GeoVaLs/"+options_.geovar.value());
}

// -----------------------------------------------------------------------------

MaxColumnGeoVaLs::~MaxColumnGeoVaLs() {}

// -----------------------------------------------------------------------------

void MaxColumnGeoVaLs::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // oops::Log::trace() << "MaxColumnGeoVaLs::compute started" << std::endl;

  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nlevs = in.nlevs(invars_[0]);

  // simulated value
  std::vector<std::vector<float>> gval(nlocs, std::vector<float>(nlevs));
  std::vector<float> gvalATLevel(nlocs);
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    // establish bottom-to-top orientation
    // const size_t level = nlevs - ilev - 1;
    // in.get(invars_[0], level, gvalATLevel);
    in.get(invars_[0], ilev, gvalATLevel);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      gval[iloc][ilev] = gvalATLevel[iloc];
    }
  }

  // oops::Log::trace() << "MaxColumnGeoVaLs::compute, assigning max to out" << std::endl;

  float gvalMax;
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    // determine maximum value in column
    gvalMax = *std::max_element(gval[iloc].begin(), gval[iloc].end());
    // same for all filter variables
    for (size_t ivar = 0; ivar < out.nvars(); ++ivar) {
      out[ivar][iloc] = gvalMax;
    }
  }
  // oops::Log::trace() << "MaxColumnGeoVaLs::compute, out: " << out << std::endl;
  // oops::Log::trace() << "MaxColumnGeoVaLs::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & MaxColumnGeoVaLs::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
