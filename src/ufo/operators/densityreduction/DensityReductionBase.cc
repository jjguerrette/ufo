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
