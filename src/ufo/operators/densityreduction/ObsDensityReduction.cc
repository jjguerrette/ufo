/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/densityreduction/ObsDensityReduction.h"

#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"
#include "ufo/operators/densityreduction/ObsDensityReductionParameters.h"
#include "ufo/SampledLocations.h"
#include "ufo/ScopedDefaultGeoVaLFormatChange.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsDensityReduction> obsDensityReductionMaker_("Density Reduction");
// -----------------------------------------------------------------------------

ObsDensityReduction::ObsDensityReduction(const ioda::ObsSpace & odb,
                                         const Parameters_ & parameters)
  : ObsOperatorBase(odb), odb_(odb), locationsCalled_(false)
{
  oops::Log::trace() << "ObsDensityReduction constructor start" << std::endl;

  const eckit::LocalConfiguration & operatorConfig = parameters.oper.value();
  ObsOperatorParametersWrapper operatorParams;
  operatorParams.deserialize(operatorConfig);
  operator_ = std::unique_ptr<ObsOperatorBase>
    (ObsOperatorFactory::create(odb, operatorParams.operatorParameters));
  requiredVars_ += operator_->requiredVars();

  const DensityReductionParametersWrapper & algorithmConfig = parameters.algorithm.value();
  algorithm_ = std::unique_ptr<DensityReductionBase>
    (DensityReductionFactory::create(algorithmConfig.densityReductionName, odb_));

  oops::Log::trace() << "ObsDensityReduction constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsDensityReduction::~ObsDensityReduction() {
  oops::Log::trace() << "ObsDensityReduction destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDensityReduction::simulateObs(const GeoVaLs & gv,
                                      ioda::ObsVector & ovec,
                                      ObsDiagnostics & ydiags,
                                      const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsDensityReduction:::simulateObs start" << std::endl;

  ScopedDefaultGeoVaLFormatChange change(gv, GeoVaLFormat::REDUCED);
  operator_->simulateObs(gv, ovec, ydiags, qc_flags);
  ScopedDefaultGeoVaLFormatChange changeback(gv, GeoVaLFormat::SAMPLED);

  oops::Log::trace() << "ObsDensityReduction::simulateObs done" <<  std::endl;
}

// -----------------------------------------------------------------------------

oops::ObsVariables ObsDensityReduction::simulatedVars() const {
  return operator_->simulatedVars();
}

// -----------------------------------------------------------------------------

ObsDensityReduction::Locations_ ObsDensityReduction::locations() const {
  oops::Log::trace() << "ObsDensityReduction::locations start" << std::endl;

  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  std::vector<float> lonsModified;
  std::vector<float> latsModified;
  std::vector<util::DateTime> timesModified;
  std::vector<util::Range<size_t>> pathsGroupedByLocation;

  algorithm_->fillModifiedLocations(latsModified,
                                    lonsModified,
                                    timesModified,
                                    pathsGroupedByLocation);

  // Set the internal flag indicating that this method has been called.
  locationsCalled_ = true;

  oops::Log::trace() << "ObsDensityReduction::locations done" <<  std::endl;

  return SampledLocations_
    (std::make_unique<SampledLocations>
     (lonsModified, latsModified, timesModified, odb_.distribution(),
      std::move(pathsGroupedByLocation)));
}

void ObsDensityReduction::computeReducedVars(const oops::Variables & vars,
                                             GeoVaLs & geovals) const {
  oops::Log::trace() << "ObsDensityReduction::computeReducedVars start" << std::endl;

  oops::Variables vars_to_reduce = requiredVars_;
  const oops::Variables & gvars = geovals.getVars();

  // Ensure all variables to reduce are present in the GeoVaLs.
  for (const auto & v : vars.variables()) {
    ASSERT(gvars.has(v));
  }

  // Determine variable lengths.
  std::vector<size_t> nLevelsGeoVaLs(gvars.size());
  for (size_t i = 0; i < gvars.size(); ++i) {
    nLevelsGeoVaLs[i] = geovals.nlevs(gvars[i], GeoVaLFormat::SAMPLED);
  }

  // Add reduced GeoVaLs.
  if (geovals.getReducedVars().size() == 0) {
    geovals.addReducedVars(gvars, nLevelsGeoVaLs);
  }

  // Call the custom locations() routine if it has not already been called.
  // This is necessary to ensure that ufo ctests of this operator always enter this routine.
  if (!locationsCalled_) {
    this->locations();
  }

  // Fill reduced GeoVaLs.
  algorithm_->fillReducedGeoVaLs(geovals);

  // In this operator, GeoVaLs are set to the 'reduced' format by
  // default. That format has one GeoVaL per observation location.
  geovals.setDefaultFormat(GeoVaLFormat::REDUCED);

  oops::Log::trace() << "ObsDensityReduction::computeReducedVars done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDensityReduction::print(std::ostream & os) const {
  os << "ObsDensityReduction with operator:" << std::endl;
  os << *operator_ << std::endl;
}

}  // namespace ufo
