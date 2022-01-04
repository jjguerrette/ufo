/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MAXCOLUMNCLOUDFRACTION_H_
#define UFO_FILTERS_OBSFUNCTIONS_MAXCOLUMNCLOUDFRACTION_H_

#include <string>
#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

/// \brief maximum cloud fraction in profile

class MaxColumnCloudFraction : public ObsFunctionBase<float> {
 public:
  explicit MaxColumnCloudFraction(const eckit::LocalConfiguration);
  ~MaxColumnCloudFraction();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MAXCOLUMNCLOUDFRACTION_H_
