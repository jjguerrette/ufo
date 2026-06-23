/*
 * (C) Copyright 2023 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <string>
#include <vector>

#include "ufo/predictors/ObsMetaDataPredictor.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/missingValues.h"

#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<ObsMetaDataPredictor> makerFuncObsMetaDataPredictor_(\
"obsMetadataPredictor");

// -----------------------------------------------------------------------------

ObsMetaDataPredictor::ObsMetaDataPredictor(const Parameters_ & parameters,
const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    order_(parameters.order.value().value_or(1.0f)),
    variable_(parameters.varName),
    functional_form_(parameters.functional_form),
    multiplier_(parameters.multiplier) {
  // predictor name is a variable name
  name() = variable_;
  if (parameters.order.value() != boost::none) {
    // override the predictor name to distinguish between predictors of different orders
    name() = name() +
    (functional_form_ == FunctionalForm::POLYNOMIAL ? "" :
     functional_form_ == FunctionalForm::COS ? "_cos" : "_sin") +
    "_order_" + std::to_string(order_);
  }
}

// -----------------------------------------------------------------------------

void ObsMetaDataPredictor::compute(const ioda::ObsSpace & odb,
                                   const GeoVaLs &,
                                   const ObsDiagnostics &,
                                   const ObsBias &,
                                   ioda::ObsVector & out) const {
  const size_t nlocs = out.nlocs();
  const size_t nvars = out.nvars();

  std::vector<float> obsMetaDataPred(nlocs, 0.0);
  const int imiss = util::missingValue<int>();
  const float fmiss = util::missingValue<float>();
  const double dmiss = util::missingValue<double>();

  // retrieve the predictor

  if (odb.dtype("MetaData", variable_) == ioda::ObsDtype::Integer) {
    std::vector<int> obsMetaDataPred2(nlocs, 0);
    odb.get_db("MetaData", variable_, obsMetaDataPred2);
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      if (obsMetaDataPred2[jloc] == imiss) {
        obsMetaDataPred[jloc] = fmiss;
      } else {
        obsMetaDataPred[jloc] = static_cast<float>(obsMetaDataPred2[jloc])*1.0f;
      }
    }
  } else {
    odb.get_db("MetaData", variable_, obsMetaDataPred);
  }

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      if (obsMetaDataPred[jloc] == fmiss) {
        // missing values do not contribute to the predictor coefficient
        out[jloc*nvars+jvar] = dmiss;
      } else {
        switch (functional_form_) {
          case FunctionalForm::POLYNOMIAL:
            out[jloc*nvars+jvar] = std::pow(obsMetaDataPred[jloc]*multiplier_, order_);
            break;
          case FunctionalForm::COS:
            out[jloc*nvars+jvar] = std::cos(M_PI * obsMetaDataPred[jloc]*multiplier_ * order_);
            break;
          case FunctionalForm::SIN:
            out[jloc*nvars+jvar] = std::sin(M_PI * obsMetaDataPred[jloc]*multiplier_ * order_);
            break;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
