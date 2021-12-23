/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SETSURFACETYPE_H_
#define UFO_FILTERS_OBSFUNCTIONS_SETSURFACETYPE_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/utils/SurfaceReportConstants.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options applying to the setting of the observation operator surface type
///
class SetSurfaceTypeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SetSurfaceTypeParameters, Parameters)

 public:
  /// Minimum ice fraction required for assignment of sea-ice surface type (default 0.2)
  /// Example: For AIRS, override default and set to 0.01
  ///          MinIceFrac: 0.01
  oops::Parameter<float> MinIceFrac{"MinIceFrac", 0.2f, this};

  /// Minimum water fraction for assignment of sea surface type (default 0.99)
  /// Example: Override default and set to 0.8
  ///          MinWaterFrac: 0.8
  oops::Parameter<float> MinWaterFrac{"MinWaterFrac", 0.99f, this};

  /// Use reported Surface Types (default false)
  /// Example: To use the reported surface types set
  ///          UseReportSurface: true
  oops::Parameter<bool> UseReportSurface{"UseReportSurface", false, this};

  /// Change name of variable storing the reported surface type
  /// Default is 'land_sea@MetaData'
  /// Example: to set to SSMIS report type (which is surface_flag@MetaData) set
  ///          SurfaceReport Name: surface_flag@MetaData
  oops::Parameter<std::string> SurfaceMetaDataName{"SurfaceReport Name", "land_sea@MetaData", this};

  /// Use reported Surface Elevation (default false)
  /// Example: To use the AAPPreported elevation set
  ///          UseReportElevation: true
  oops::Parameter<bool> UseReportElevation{"UseReportElevation", false, this};

  /// Use AAPP Surface Types (default false)
  /// Example: To use the AAPP Surface Type set
  ///          UseAAPPSurfaceType: true
  oops::Parameter<bool> UseAAPPSurfaceClass{"UseAAPPSurfaceClass", false, this};

  /// Use reported Surface Water Fraction (default false)
  /// Example: To use reported Surface Water Fraction set
  ///          UseSurfaceWaterFraction: true
  oops::Parameter<bool> UseSurfaceWaterFraction{"UseSurfaceWaterFraction", false, this};

  /// Assumed limit of seaice regardless of ice fraction
  oops::Parameter<float> IceLimitHard{"IceLimitHard", 72.0, this};

  /// Assumed limit of seaice where any sea ice is observed
  oops::Parameter<float> IceLimitSoft{"IceLimitSoft", 55.0, this};

  /// Surface height above which surface is deemed to be Highland for purposes of QC
  oops::Parameter<float> HighlandHeight{"HighlandHeight", 1000.0, this};

  /// Default Obs Operator surface types to map to (RTTOV is the default here)
  /// Example: To use the RTTOV sea surface type (1) as the ObsOperator sea surface type set
  ///          SurfaceTypeSea: 1
  oops::Parameter<int> SurfaceTypeLand{"SurfaceTypeLand", 0, this};
  oops::Parameter<int> SurfaceTypeSea{"SurfaceTypeSea", 1, this};
  oops::Parameter<int> SurfaceTypeSeaIce{"SurfaceTypeSeaIce", 2, this};
};

///
/// \brief Set observation operator surface type based on model and observation data
///
class SetSurfaceType : public ObsFunctionBase<float> {
 public:
  explicit SetSurfaceType(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~SetSurfaceType() override;

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const override;

  const ufo::Variables & requiredVariables() const override;

 private:
  ufo::Variables invars_;
  SetSurfaceTypeParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SETSURFACETYPE_H_
