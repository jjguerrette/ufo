/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/densityreduction/DensityReductionBase.h"

#include <map>
#include <string>

#include "oops/util/Logger.h"

namespace ufo {

DensityReductionBase::DensityReductionBase(const DensityReductionParametersBase & params,
                                           const ioda::ObsSpace & odb)
  : odb_(odb)
{}

void DensityReductionBase::fillReducedGeoVaLs(GeoVaLs & geovals) const {
  const oops::Variables & gvars = geovals.getVars();

  // Lengths of each GeoVaL.
  std::vector<size_t> nLevelsGeoVaLs(gvars.size());
  for (size_t i = 0; i < gvars.size(); ++i) {
    nLevelsGeoVaLs[i] = geovals.nlevs(gvars[i], GeoVaLFormat::SAMPLED);
  }

  // Retrieve pathsGroupedByLocation.
  std::vector<util::Range<size_t>> pathsGroupedByLocation;
  geovals.getProfileIndicesGroupedByLocation(gvars[0], pathsGroupedByLocation);
  const size_t nlocs = pathsGroupedByLocation.size();

  // Fill in reduced vectors of GeoVaLs.
  std::vector<double> sampledValues(numSamples_);
  std::vector<double> reducedValues(nlocs);
  for (size_t i = 0; i < gvars.size(); ++i) {
    const auto var = gvars[i];
    for (size_t jlev = 0; jlev < nLevelsGeoVaLs[i]; ++jlev) {
      geovals.getAtLevel(sampledValues, var, jlev, GeoVaLFormat::SAMPLED);
      for (size_t jloc = 0; jloc < nlocs; ++jloc) {
        reducedValues[jloc] = sampledValues[pathsGroupedByLocation[jloc].begin];
      }
      geovals.putAtLevel(reducedValues, var, jlev, GeoVaLFormat::REDUCED);
    }
  }
}

DensityReductionFactory::DensityReductionFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end())
    throw eckit::BadParameter
      (name + " already registered in ufo::DensityReductionFactory.", Here());
  getMakers()[name] = this;
}

std::unique_ptr<DensityReductionBase>
DensityReductionFactory::create(const DensityReductionParametersBase & params,
                                const ioda::ObsSpace & odb) {
  oops::Log::trace() << "DensityReductionBase::create starting" << std::endl;
  const std::string & name = params.reductionName;
  typename std::map<std::string, DensityReductionFactory*>::iterator jloc = getMakers().find(name);
  if (jloc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in ufo::DensityReductionFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<DensityReductionBase> ptr =
    jloc->second->make(params, odb);
  oops::Log::trace() << "DensityReductionBase::create done" << std::endl;
  return ptr;
}

std::unique_ptr<DensityReductionParametersBase>
DensityReductionFactory::createParameters(const std::string & name) {
  oops::Log::trace() << "DensityReductionBase::createParameters starting" << std::endl;
  typename std::map<std::string, DensityReductionFactory*>::iterator jloc = getMakers().find(name);
  if (jloc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in ufo::DensityReductionFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<DensityReductionParametersBase> ptr = jloc->second->makeParameters();
  oops::Log::trace() << "DensityReductionBase::createParameters done" << std::endl;
  return ptr;
}

}  // namespace ufo
