/*
 * (C) Copyright 2025 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_DENSITYREDUCTION_OBSDENSITYREDUCTIONPARAMETERS_H_
#define UFO_OPERATORS_DENSITYREDUCTION_OBSDENSITYREDUCTIONPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

class DensityReductionParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DensityReductionParametersWrapper, oops::Parameters)
 public:
  /// Name of the density reduction algorithm.
  /// Valid names are specified using a `DensityReductionMaker` in subclasses
  /// of DensityReductionBase.

  oops::RequiredPolymorphicParameter<DensityReductionParametersBase, DensityReductionFactory>
    densityReductionName{"name", this};
};

/// Configuration options recognized by the density reduction operator.
class ObsDensityReductionParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsDensityReductionParameters, ObsOperatorParametersBase)

 public:
  /// Configuration options for the underlying observation operator.
  oops::RequiredParameter<eckit::LocalConfiguration> oper{"operator", this};

  /// Density reduction algorithm.
  oops::RequiredParameter<DensityReductionParametersWrapper> algorithm{"algorithm", this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_DENSITYREDUCTION_OBSDENSITYREDUCTIONPARAMETERS_H_
