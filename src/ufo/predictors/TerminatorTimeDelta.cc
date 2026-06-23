/*
 * (C) Copyright 2025 Tomorrow.io
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/predictors/TerminatorTimeDelta.h"

namespace ufo {

static PredictorMaker<TerminatorTimeDelta>
       makerFuncTerminatorTimeDelta_("terminator_time_delta");

// -----------------------------------------------------------------------------

TerminatorTimeDelta::TerminatorTimeDelta(const Parameters_ & parameters,
  const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    order_(parameters.order),
    day_night_(parameters.day_night),
    functional_form_(parameters.functional_form) {
  // override the predictor name to distinguish between TimeSinceTerminator predictors of
  // different orders, functional form, and three day or night possibilities.
  name() = name() +
      (day_night_ == DayNightType::NIGHT ? "_night" :
       day_night_ == DayNightType::DAY ? "_day" : "_night_and_day") +
      (functional_form_ == FunctionalForm::POLYNOMIAL ? "_polynomial" :
       functional_form_ == FunctionalForm::COS ? "_cos" : "_sin") +
      "_order_" + std::to_string(order_);
}

// -----------------------------------------------------------------------------
void TerminatorTimeDelta::compute(const ioda::ObsSpace & odb,
                                             const GeoVaLs &,
                                             const ObsDiagnostics &,
                                             const ObsBias &,
                                             ioda::ObsVector & out) const {
  const std::size_t nlocs = out.nlocs();
  const std::size_t nvars = out.nvars();
  const float fmiss = util::missingValue<float>();
  const double dmiss = util::missingValue<double>();

  std::vector<int> day_or_night_qualifier(nlocs, 0);
  odb.get_db("MetaData", "dayOrNightQualifier", day_or_night_qualifier);
  // dayOrNightQualifier:description = "True (1) if earth is between the sun and spacecraft
  // (night), False (0) if earth is between the spacecraft and sun (day)"

  std::vector<float> delta_t_since_last(nlocs, 0.0);
  odb.get_db("MetaData", "timeSinceLastTerminatorCrossing", delta_t_since_last);
  // timeSinceLastTerminatorCrossing:description = "Time since last terminator crossing (seconds)"

  std::vector<float> delta_t_to_next(nlocs, 0.0);
  odb.get_db("MetaData", "timeToNextTerminatorCrossing", delta_t_to_next);
  // timeToNextTerminatorCrossing:description = "Time to next terminator crossing (seconds)"

  // KLUDGE: due to limitations in the converter, delta_t_to_next and delta_t_since_last are
  // sometimes negative. Either case is unphysical.  When this happens, and the other delta_t is
  // equal to its opposite, we get a total_delta_t of 0.0, and an incorrect sign for
  // delta_t_predictor. The solution is to set both to zero in this case, then handle the case of
  // total_delta_t == 0.0 below.
  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (delta_t_since_last[jloc] != fmiss && delta_t_to_next[jloc] != fmiss) {
      if (delta_t_since_last[jloc] == -delta_t_to_next[jloc]) {
        delta_t_since_last[jloc] = 0.0;
        delta_t_to_next[jloc] = 0.0;
      }
      // check that delta_t's are non-negative besides the case above
      if (delta_t_since_last[jloc] < 0.0) {
        throw eckit::Exception("delta_t_since_last is negative: " +
          std::to_string(delta_t_since_last[jloc]));
      }
      if (delta_t_to_next[jloc] < 0.0) {
        throw eckit::Exception("delta_t_to_next is negative: " +
          std::to_string(delta_t_to_next[jloc]));
      }
    }
  }

  float delta_t_predictor;
  float total_delta_t;
  double predictor_value;

  // only normalize delta_t if the functional form is cosine or sine,
  // which places the predictor in the range of -1 to 1
  bool normalize = (
    functional_form_ == FunctionalForm::COS ||
    functional_form_ == FunctionalForm::SIN);

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    // the delta_t_predictor is:
    // * < 0 if the satellite is in the night hemisphere
    // * > 0 if the satellite is in the day hemisphere
    // * == 0 if the satellite is at the night-to-day terminator
    // when normalized, the delta_t_predictor:
    // * varies between -1 and 1
    // * == ±1 if the satellite is at the day-to-night terminator
    // when not normalized, the delta_t_predictor is converted to hours

    delta_t_predictor = fmiss;
    if (delta_t_since_last[jloc] != fmiss && delta_t_to_next[jloc] != fmiss) {
      if (day_or_night_qualifier[jloc] == 1) {
        // night (negative)
        delta_t_predictor = - delta_t_to_next[jloc];
      } else if (day_or_night_qualifier[jloc] == 0) {
        // day (positive)
        delta_t_predictor = delta_t_since_last[jloc];
      }
      if (normalize) {
        total_delta_t = delta_t_since_last[jloc] + delta_t_to_next[jloc];
        if (total_delta_t > 0.0) {
          delta_t_predictor /= total_delta_t;
        }
      } else {
        delta_t_predictor /= 3600.0;
      }
    }

    // Calculate the predictor value
    predictor_value = fmiss;
    if (delta_t_predictor != fmiss) {
      // Check that delta_t_predictor is within the range of -1 to 1
      if (normalize && (delta_t_predictor < -1.0 || delta_t_predictor > 1.0)) {
        throw eckit::Exception("delta_t_predictor is out of range: " +
          std::to_string(delta_t_predictor));
      }
      switch (functional_form_) {
        case FunctionalForm::POLYNOMIAL:
          predictor_value = std::pow(delta_t_predictor, order_);
        break;
      case FunctionalForm::COS:
        // scale x by pi to get a range of -pi to pi
        predictor_value = std::cos(M_PI * delta_t_predictor * order_);
        break;
      case FunctionalForm::SIN:
        // scale x by pi to get a range of -pi to pi
          predictor_value = std::sin(M_PI * delta_t_predictor * order_);
        break;
      }
    }

    // Check that predictor_value is within the range of -1 to 1
    if (normalize && (predictor_value < -1.0 || predictor_value > 1.0)) {
      throw eckit::Exception("predictor_value is out of range: " +
        std::to_string(predictor_value));
    }

    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      if (predictor_value == fmiss) {
        out[jloc*nvars+jvar] = dmiss;
      } else if (
          (day_or_night_qualifier[jloc] == 1 && day_night_ == DayNightType::NIGHT) ||
          (day_or_night_qualifier[jloc] == 0 && day_night_ == DayNightType::DAY) ||
          day_night_ == DayNightType::BOTH) {
        out[jloc*nvars+jvar] = predictor_value;
      } else {
        out[jloc*nvars+jvar] = 0.0;
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
