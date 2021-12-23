/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_POLYNOMIAL2D_H_
#define UFO_FILTERS_OBSFUNCTIONS_POLYNOMIAL2D_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------

/// \brief Options controlling Polynomial2D function
class Polynomial2DParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(Polynomial2DParameters, Parameters)

 public:

  /// independent Parameters for each term of the polynomial
  oops::RequiredParameter<std::vector<Polynomial2DTermParameters>>
    fittingTerms{"fitting terms", this};

  /// independent coefficients for each filter variable
  oops::RequiredParameter<std::vector<LinearFit2DCoefficientsParameters>>
    fittingCoefficients{"fitting coefficients", this};

  /// x and y variables in the 2D polynomial
  oops::RequiredParameter<Variable> xvar{"xvar", this};
  oops::RequiredParameter<Variable> yvar{"yvar", this};

  /// filter variables, must be same as the parent filter variables,
  /// including missing entirely if that is the case. Only the filter
  /// variables will be processed, such that non-matching fitting coefficients
  /// vector members will be ignored.  See Polynomial2D for examples.
  oops::OptionalParameter<std::vector<Variable>> filterVariables{
    "filter variables", this};
};

class Polynomial2DTermParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(Polynomial2DTermParameters, Parameters)

 public:
  /// exponents of x and y for a single term in the polynomial
  /// should be a vector of size-2
  oops::RequiredParameter<std::vector<float>> exponents{"exponents", this};
}:

class LinearFit2DCoefficientsParameters: public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearFit2DCoefficientsParameters, Parameters)

 public:
  /// filter variable to which these polynomial coefficients apply
  oops::RequiredParameter<std::string>> variable{"name", this};

  /// single channel, if applicable
  oops::OptionalParameter<int> channel{"channel", this};

  /// coefficients of the polynomial, ordered identically to "fitting terms"
  /// and of the same length
  oops::RequiredParameter<std::vector<float>> values{"values", this};
}:

// -----------------------------------------------------------------------------

/// \brief 2D polynomial function
///
/// 1 + c1 x + c2 y + c3 x^2 + c4 xy + c5 y^2 + ...
/// Notes:

// -----------------------------------------------------------------------------

/// ### example configuration ###
///   Notes:
///   () The pseudo-yaml example below is for a single contiguous configuration file with anchors
///   () Although the examples below are for the Perform Action filter with "assign error" as the
///      action, this ObsFunction can be used in any application that requires a quantity
///      parameterized as a 2D polynomial.
///   () Make sure the "filter variables" are prescribed identically to the parent filter (e.g.,
///      using yaml anchors like "*FilterVariables" below).
///   () The example below is for a three-term degree-1 polynomial.  More terms are added by
///      extending options.(fitting terms) and options.(fitting coefficients)[:].values.
///   () Alternatively, the output can be parameterized on a single variable, x or y, by
///      setting the exponents of the other term, y or x, repsectively, to zero in all terms.
///
///      _polynomial 2d options anchor: &Polynomial2DAnchor
///        fitting terms:
///        - exponents: [0, 0]
///        - exponents: [1, 0]
///        - exponents: [0, 1]
///        fitting coefficients:
///        - name: brightness_temperature
///          channel: 1
///          values: [0., 1., 2.5]
///        - name: brightness_temperature
///          channel: 2
///          values: [0., 1., 2.5]
///        - name: brightness_temperature # never selected in below examples
///          channel: 5
///          values: [0., 1., 2.5]
///        - name: air_temperature
///          values: [1., 1., 1.5]
///        - name: specific_humidity
///          values: [1., 0., 1.5]
///
///     # (1) channels 1 and 2 for brightness_temperature are selected
///     - filter: Perform Action
///       filter variables: &radianceFilterVariables
///       - name: brightness_temperature
///         channels: 1,2
///       action:
///         name: assign error
///         error function:
///           name: Polynomial2D@ObsFunction
///           options: &Polynomial2DAnchor
///             filter variables: *radianceFilterVariables
///             xvar: XRadGroupName/XRadVariableName
///             yvar: YRadGroupName/YRadVariableName
///
///     # (2) air_temperature and specific_humidity are selected
///
///     - filter: Perform Action
///       filter variables: &conventionalFilterVariables
///       - name: air_temperature
///       - name: specific_humidity
///       action:
///         name: assign error
///         error function:
///           name: Polynomial2D@ObsFunction
///           options: *Polynomial2DAnchor
///             filter variables: *conventionalFilterVariables
///             xvar: XConvGroupName/XConvVariableName
///             yvar: YConvGroupName/YConvVariableName
///
///     # (3) defaults to simulated variables in obs operator
///
///     - filter: Perform Action
///       action:
///         name: assign error
///         error function:
///           name: Polynomial2D@ObsFunction
///           options: *Polynomial2DAnchor
///             xvar: XSimulatedGroupName/XSimulatedVariableName
///             yvar: YSimulatedGroupName/YSimulatedVariableName

class Polynomial2D : public ObsFunctionBase<float> {
 public:
  explicit Polynomial2D(const eckit::LocalConfiguration);
  ~Polynomial2D();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  Polynomial2DParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_POLYNOMIAL2D_H_
