/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SatWindsLNVDCheck.h"

#include <math.h>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<SatWindsLNVDCheck> makerObsFuncSatWindsLNVDCheck_("SatWindsLNVDCheck");

// -----------------------------------------------------------------------------

SatWindsLNVDCheck::SatWindsLNVDCheck(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::debug() << "SatWindsLNVDCheck: config = " << conf << std::endl;
  // Initialize options
  options_.deserialize(conf);

  // We need to retrieve the observed wind components.
  invars_ += Variable("eastward_wind@ObsValue");
  invars_ += Variable("northward_wind@ObsValue");

  // Typical use would be HofX group, but during testing, we include option for GsiHofX
  std::string test_hofx = options_.test_hofx.value();
  invars_ += Variable("eastward_wind@" + test_hofx);
  invars_ += Variable("northward_wind@" + test_hofx);

  // TODO(gthompsn): Need to include a check that whatever HofX group name used actually exists.
}

// -----------------------------------------------------------------------------

SatWindsLNVDCheck::~SatWindsLNVDCheck() {}

// -----------------------------------------------------------------------------

void SatWindsLNVDCheck::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue(missing);

  // Ensure that only one output variable is expected.
  ASSERT(out.nvars() == 1);

  // Retrieve SatWinds observations of wind components
  std::vector<float> u, v;
  in.get(Variable("eastward_wind@ObsValue"), u);
  in.get(Variable("northward_wind@ObsValue"), v);
  // Retrieve Model HofX wind components
  std::string test_hofx = options_.test_hofx.value();
  std::vector<float> um, vm;
  in.get(Variable("eastward_wind@" + test_hofx), um);
  in.get(Variable("northward_wind@" + test_hofx), vm);

  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (u[jj] != missing && v[jj] != missing) {
      if ((u[jj]*u[jj] + v[jj]*v[jj]) > 1.01f) {
        out[0][jj] = sqrt((u[jj]-um[jj])*(u[jj]-um[jj]) + (v[jj]-vm[jj])*(v[jj]-vm[jj]))
                    / log(sqrt(u[jj]*u[jj] + v[jj]*v[jj]));
      } else {
        out[0][jj] = sqrt((u[jj]-um[jj])*(u[jj]-um[jj]) + (v[jj]-vm[jj])*(v[jj]-vm[jj]));
      }
    } else {
      out[0][jj] = missing;
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SatWindsLNVDCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
