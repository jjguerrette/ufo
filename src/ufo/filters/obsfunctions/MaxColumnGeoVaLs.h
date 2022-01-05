/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MAXCOLUMNGEOVALS_H_
#define UFO_FILTERS_OBSFUNCTIONS_MAXCOLUMNGEOVALS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

/// \brief parameters controlling MaxColumnGeoVaLs

class MaxColumnGeoVaLsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MaxColumnGeoVaLsParameters, Parameters)

 public:
  /// GeoVaLs variable name
  oops::RequiredParameter<std::string> geovar{"geophysical variable", this};
};

// -----------------------------------------------------------------------------

/// \brief maximum cloud fraction in profile

class MaxColumnGeoVaLs : public ObsFunctionBase<float> {
 public:
  explicit MaxColumnGeoVaLs(const eckit::LocalConfiguration);
  ~MaxColumnGeoVaLs();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  MaxColumnGeoVaLsParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MAXCOLUMNGEOVALS_H_
