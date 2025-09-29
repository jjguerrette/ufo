/*
 * (C) Crown copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORRECONDITIONER_H_
#define UFO_ERRORS_OBSERRORRECONDITIONER_H_

#include <Eigen/Core>

#include "ioda/ObsVector.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace ufo {

enum class ObsErrorReconditionerMethod {
  MINIMUMEIGENVALUE, RIDGEREGRESSION, NORECONDITIONING
};

struct ObsErrorReconditionerMethodParameterTraitsHelper {
  typedef ObsErrorReconditionerMethod EnumType;
  static constexpr char enumTypeName[] = "ReconditionMethod";
  static constexpr util::NamedEnumerator<ObsErrorReconditionerMethod> namedValues[] = {
    { ObsErrorReconditionerMethod::MINIMUMEIGENVALUE, "Minimum Eigenvalue" },
    { ObsErrorReconditionerMethod::RIDGEREGRESSION, "Ridge Regression" },
    { ObsErrorReconditionerMethod::NORECONDITIONING, "No reconditioning" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::ObsErrorReconditionerMethod> :
    public EnumParameterTraits<ufo::ObsErrorReconditionerMethodParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// \brief Parameters for Recondition method
class ObsErrorReconditionerParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(ObsErrorReconditionerParameters, Parameters)
 public:
  /// Method of reconditioning, either Ridge regression or minimum eigenvalue
  oops::Parameter<ObsErrorReconditionerMethod> ReconMethod{"recondition method",
                                                 "recondition method (options:"
                                                 " Minimum Eigenvalue,"
                                                 " Ridge Regression,"
                                                 " No reconditioning)",
                                                 ObsErrorReconditionerMethod::NORECONDITIONING,
                                                 this};
  /// Target fraction of condition number
  oops::OptionalParameter<float> kFrac{"fraction", this};
  /// Threshold for Minimum EigenValue
  oops::OptionalParameter<float> Threshold{"threshold", this};
  /// Shift for Ridge Regression
  oops::OptionalParameter<float> Shift{"shift", this};
};

class ObsErrorReconditioner {
 public:
  /// The type of parameters for this class.
  typedef ObsErrorReconditionerParameters Parameters_;

  explicit ObsErrorReconditioner(const Parameters_ &);

  // method
  void recondition(Eigen::MatrixXd &R) const;

 private:
  /// Configuration as a Parameters_
  Parameters_ params_;
};

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORRECONDITIONER_H_
