/*
 * (C) Copyright 2025 Tomorrow.io
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_TERMINATORTIMEDELTA_H_
#define UFO_PREDICTORS_TERMINATORTIMEDELTA_H_

#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/Parameter.h"

#include "ufo/predictors/PredictorBase.h"

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

enum class DayNightType {
  DAY, NIGHT, BOTH
};

struct DayNightTypeParameterTraitsHelper {
  typedef DayNightType EnumType;
  static constexpr char enumTypeName[] = "DayNightType";
  static constexpr util::NamedEnumerator<DayNightType> namedValues[] = {
    { DayNightType::DAY, "day" },
    { DayNightType::NIGHT, "night" },
    { DayNightType::BOTH, "both" }
  };
};

enum class FunctionalForm {
  POLYNOMIAL, COS, SIN
};

struct FunctionalFormParameterTraitsHelper {
  typedef FunctionalForm EnumType;
  static constexpr char enumTypeName[] = "FunctionalForm";
  static constexpr util::NamedEnumerator<FunctionalForm> namedValues[] = {
    { FunctionalForm::POLYNOMIAL, "polynomial" },
    { FunctionalForm::COS, "cos" },
    { FunctionalForm::SIN, "sin" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::DayNightType> :
    public EnumParameterTraits<ufo::DayNightTypeParameterTraitsHelper>
{};

template <>
struct ParameterTraits<ufo::FunctionalForm> :
    public EnumParameterTraits<ufo::FunctionalFormParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

// -----------------------------------------------------------------------------

/// Configuration parameters of the TerminatorTimeDelta predictor.
class TerminatorTimeDeltaParameters : public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(TerminatorTimeDeltaParameters, PredictorParametersBase);

 public:
  /// Order of the term.  For polynomial, this is the exponent.  For cos and sin, this is the
  /// order of the term in the Fourier series (coefficient of the operand, e.g. cos(x), cos(2x,
  /// etc...).
  oops::RequiredParameter<float> order{"order", this};
  /// Controls whether this predictor applies to satellite positions that are or are not in Earth's
  /// shadow ("night" or "day" or "both"). Default is "both".
  oops::Parameter<DayNightType> day_night{"day or night", DayNightType::BOTH, this};
  /// Functional form of the predictor, either polynomial, cos, or sin.
  oops::RequiredParameter<FunctionalForm> functional_form{"functional form", this};
};

// -----------------------------------------------------------------------------

/**
 * This predictor is used to fit residual errors as a function of time deltas since the
 * last, and to the next, terminator crossing. The data must contain these time deltas
 * in the variables "MetaData/timeSinceLastTerminatorCrossing" and
 * "MetaData/timeToNextTerminatorCrossing", respectively. Three member variables are
 * used to store the order of the term in the series being calculated (order_), the
 * day/night selector, and the functional form, respectively. These are read from the
 * yaml configuration file.
 *
 * This predictor is particularly useful when the fractions of time a sattelite spends
 * in day and night are variable, as in a non-sun-synchronous orbit with a transient solar
 * beta angle. It may also be useful when those fractions are nearly constant, but unequal,
 * as in a sun-synchronous orbit. The functional form is left open to the user in case additional
 * flexibility is needed for future satellite missions.
 */

class TerminatorTimeDelta : public PredictorBase {
 public:
  /// The type of parameters accepted by the constructor of this predictor.
  /// This typedef is used by the PredictorFactory.
  typedef TerminatorTimeDeltaParameters Parameters_;

  TerminatorTimeDelta(const Parameters_ &, const oops::Variables &);

  void compute(const ioda::ObsSpace &,
               const GeoVaLs &,
               const ObsDiagnostics &,
               const ObsBias &,
               ioda::ObsVector &) const override;

 private:
  double order_;
  DayNightType day_night_;
  FunctionalForm functional_form_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_TERMINATORTIMEDELTA_H_
