/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_DENSITYREDUCTION_DENSITYREDUCTIONBASE_H_
#define UFO_OPERATORS_DENSITYREDUCTION_DENSITYREDUCTIONBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/base/Locations.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/HasParameters_.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variables.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace ufo {

/// \brief DensityReduction parameters base class.
class DensityReductionParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(DensityReductionParametersBase, oops::Parameters)
 public:
  /// Name of the density reduction algorithm.
  /// Valid names are specified using a `DensityReductionMaker` in subclasses
  /// of DensityReductionBase.
  oops::RequiredParameter<std::string> reductionName{"name", this};
};

/// \brief Concrete class containing the options specified by the DensityReductionParametersBase.
class GenericDensityReductionParameters : public DensityReductionParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericDensityReductionParameters, DensityReductionParametersBase)
};

/// \brief DensityReduction base class.
///
/// Subclasses of this class must implement the `fillModifiedLocations` and
/// `fillReducedGeoVaLs` methods.
class DensityReductionBase {
 public:
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  explicit DensityReductionBase(const DensityReductionParametersBase &,
                                const ioda::ObsSpace &);
  virtual ~DensityReductionBase() {}

  virtual void fillModifiedLocations(std::vector<float> &,
                                     std::vector<float> &,
                                     std::vector<util::DateTime> &,
                                     std::vector<util::Range<size_t>> &) const = 0;

  virtual void fillReducedGeoVaLs(GeoVaLs &) const;

 protected:
  const ioda::ObsSpace & odb_;

  mutable int numSamples_;

  void updateNumSamples(const int n) const {
    numSamples_ = n;
  }
};

/// DensityReduction factory.
class DensityReductionFactory {
 public:
  static std::unique_ptr<DensityReductionBase> create(const DensityReductionParametersBase &,
                                                      const ioda::ObsSpace &);

  static std::unique_ptr<DensityReductionParametersBase> createParameters(const std::string &name);

  /// \brief Return the names of all subclasses that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~DensityReductionFactory() = default;

 protected:
  explicit DensityReductionFactory(const std::string &);

 private:
  virtual std::unique_ptr<DensityReductionBase> make(const DensityReductionParametersBase &,
                                                     const ioda::ObsSpace &) = 0;

  virtual std::unique_ptr<DensityReductionParametersBase> makeParameters() const = 0;

  static std::map <std::string, DensityReductionFactory*> & getMakers() {
    static std::map <std::string, DensityReductionFactory*> makers_;
    return makers_;
  }
};

template<class T>
class DensityReductionMaker : public DensityReductionFactory {
  typedef oops::TParameters_IfAvailableElseFallbackType_t<T, GenericDensityReductionParameters>
    Parameters_;

  std::unique_ptr<DensityReductionBase> make(const DensityReductionParametersBase & params,
                                             const ioda::ObsSpace & odb) override {
    const auto & stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::unique_ptr<DensityReductionBase>
      (new T(stronglyTypedParams, odb));
  }

  std::unique_ptr<DensityReductionParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit DensityReductionMaker(const std::string & name)
    : DensityReductionFactory(name) {}
};

}  // namespace ufo

#endif  // UFO_OPERATORS_DENSITYREDUCTION_DENSITYREDUCTIONBASE_H_
