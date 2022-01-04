/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/MaxColumnCloudFraction.h"

#include <cmath>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<MaxColumnCloudFraction> makerMAXCldFrac_("MaxColumnCloudFraction");

// -----------------------------------------------------------------------------

MaxColumnCloudFraction::MaxColumnCloudFraction(const eckit::LocalConfiguration config)
  : invars_() {
  // Include required variables
  invars_ += Variable("cloud_area_fraction_in_atmosphere_layer@GeoVaLs");
}

// -----------------------------------------------------------------------------

MaxColumnCloudFraction::~MaxColumnCloudFraction() {}

// -----------------------------------------------------------------------------

void MaxColumnCloudFraction::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {

  // Get dimensions
  size_t varsize = out.nvars();
  size_t nlocs = in.nlocs();
  size_t nlevs = in.nlevs(Variable("cloud_area_fraction_in_atmosphere_layer@GeoVaLs"));

  // simulated cloud fraction
  std::vector<std::vector<float>> CFx(nlocs, std::vector<float>(nlevs));
  std::vector<float> CFatLevel(nlocs);
  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
    const size_t level = nlevs - ilev - 1;
    in.get(Variable("cloud_area_fraction_in_atmosphere_layer@GeoVaLs"), level, CFatLevel);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      CFx[iloc][ilev] = CFatLevel[iloc];
    }
  }

  double CFxmax;
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    // determine maximum cloud fraction in column
    CFxmax = *std::max_element(CFx[iloc].begin(), CFx[iloc].end());
    // same for all variables
    for (size_t ivar = 0; ivar < varsize; ++ivar) {
      out[ivar][iloc] = CFxmax;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & MaxColumnCloudFraction::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
