# (C) Crown copyright 2025 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( densityreduction_files
      DensityReductionBase.cc
      DensityReductionBase.h
      DensityReductionSimple.cc
      DensityReductionSimple.h
      ObsDensityReduction.cc
      ObsDensityReduction.h
      ObsDensityReductionParameters.h
      ObsDensityReductionTLAD.h
      ObsDensityReductionTLAD.cc
)
PREPEND( _p_densityreduction_files     "operators/densityreduction"     ${densityreduction_files} )

set ( densityreduction_src_files
      ${_p_densityreduction_files}
      PARENT_SCOPE
)
