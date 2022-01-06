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

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

// -----------------------------------------------------------------------------

/// \brief parameters controlling a single term in a Polynomial2D equation

class Polynomial2DTermParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(Polynomial2DTermParameters, Parameters)

 public:
  /// exponents of x and y for a single term in the polynomial
  /// should be a vector of size-2
  oops::RequiredParameter<std::vector<int>> exponents{"exponents", this};
};

// -----------------------------------------------------------------------------

/// \brief coefficients for terms in a linear fit

class LinearFitCoefficientsParameters: public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearFitCoefficientsParameters, Parameters)

 public:
  /// filter variable to which these polynomial coefficients apply
  oops::RequiredParameter<std::string> name{"name", this};

  /// single channel, if applicable.  omit for non-channel filter variables.
  oops::OptionalParameter<int> channel{"channel", this};

  /// coefficients of the polynomial, ordered identically to "polynomial terms"
  /// and of the same length or smaller. only the first values.size() "polynomial terms"
  /// will be used to form the output for this filter variable.
  oops::RequiredParameter<std::vector<double>> values{"values", this};
};

// -----------------------------------------------------------------------------

/// \brief parameters that determine a full Polynomial2D equation
class Polynomial2DParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(Polynomial2DParameters, Parameters)

 public:
  /// filter variables, must be same as the parent filter variables,
  /// including missing entirely if that is the case. Only the filter
  /// variables will be processed, such that non-matching 'polynomial coefficients'
  /// members will be ignored.  See below for examples.
  oops::OptionalParameter<std::vector<Variable>> filterVariables{
    "filter variables", this};

  /// x and y variables in the 2D polynomial
  oops::RequiredParameter<Variable> xvar{"xvar", this};
  oops::RequiredParameter<Variable> yvar{"yvar", this};

  /// independent Parameters for each term of the polynomial
  oops::RequiredParameter<std::vector<Polynomial2DTermParameters>>
    polynomialTerms{"polynomial terms", this};

  /// independent coefficients for each filter variable
  oops::RequiredParameter<std::vector<LinearFitCoefficientsParameters>>
    polynomialCoefficients{"polynomial coefficients", this};
};

// -----------------------------------------------------------------------------

/// \brief 2D polynomial equation of any degree
///
/// 1 + c1 x + c2 y + c3 x^2 + c4 xy + c5 y^2 + ...
///
/// ### important notes ###
/// () The pseudo-yaml example below is for a single contiguous configuration file with anchors
/// () While the example below only has one Polynomial2DAnchor, multiple such anchors are likely
///    needed, one fore each observation type.  That requirement is driven by overlap in
///    filter variable names and/or channel numbers between observation types.
/// () Make sure the "filter variables" are prescribed identically to the parent filter (e.g.,
///    using yaml anchors like "*FilterVariables" below).
/// () The example below includes a 6-term degree-2 polynomial. The first 2 brightness_temperature
///    channels use all 6 of the terms. Channel 5, air_temperature, and specific_humidity use
///    only the leading 3 of those 6 terms, equivalent to a linear fit in 2-dimensional space.
///    More terms are added by extending options.(polynomial terms) and
///    options.(polynomial coefficients)[:].values.
/// () Alternatively, the output can be parameterized on a single variable, x or y, by
///    setting the exponents of the other term, y or x, repsectively, to zero in all terms.
/// () Although the examples below are for the Perform Action filter with "assign error" as the
///    action, this ObsFunction can be used in any application that requires a quantity
///    parameterized as a 2D polynomial.
///
/// ### example configuration ###
///    _polynomial 2d options anchor: &Polynomial2DAnchor
///      polynomial terms:
///      - exponents: [0, 0]
///      - exponents: [1, 0]
///      - exponents: [0, 1]
///      - exponents: [2, 0]
///      - exponents: [1, 1]
///      - exponents: [0, 2]
///      polynomial coefficients:
///      - name: brightness_temperature
///        channel: 1
///        values: [0., 1., 2.5, 0.25, 0.1, 0.3]
///      - name: brightness_temperature
///        channel: 2
///        values: [0., 1., 2.5, 0.5, 0.1, 0.5]
///      - name: brightness_temperature # never selected in below examples
///        channel: 5
///        values: [0., 1., 2.5]
///      - name: air_temperature
///        values: [1., 1., 1.5]
///      - name: specific_humidity
///        values: [1., 0., 1.5]
///
///   observations:
///   # (1) channels 1 and 2 for brightness_temperature are selected
///   - obs filters:
///     - filter: Perform Action
///       filter variables: &radianceFilterVariables
///       - name: brightness_temperature
///         channels: 1,2
///       action:
///         name: assign error
///         error function:
///           name: Polynomial2D@ObsFunction
///           channels: 1,2
///           options: &Polynomial2DAnchor
///             filter variables: *radianceFilterVariables
///             xvar: MetaData/latitude
///             yvar: GeoVaLs/surface_wind_speed
///
///   # (2) air_temperature and specific_humidity are selected individually
///   - obs filters:
///     - filter: Perform Action
///       filter variables: &airtempFilterVariables
///       - name: air_temperature
///       action:
///         name: assign error
///         error function:
///           name: Polynomial2D@ObsFunction
///           options: *Polynomial2DAnchor
///             filter variables: *airtempFilterVariables
///             xvar: XConvGroupName/XConvVariableName
///             yvar: YConvGroupName/YConvVariableName
///     - filter: Perform Action
///       filter variables: &spechumFilterVariables
///       - name: specific_humidity
///       action:
///         name: assign error
///         error function:
///           name: Polynomial2D@ObsFunction
///           options: *Polynomial2DAnchor
///             filter variables: *spechumFilterVariables
///             xvar: XConvGroupName/XConvVariableName
///             yvar: YConvGroupName/YConvVariableName
///
///   # (3) defaults to simulated variables in obs operator, can only handle a single variable
///   - obs filters:
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
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_POLYNOMIAL2D_H_
