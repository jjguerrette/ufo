/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/Polynomial2D.h"

#include <cmath>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

static ObsFunctionMaker<Polynomial2D> makerSCIIR_("Polynomial2D");

// -----------------------------------------------------------------------------

Polynomial2D::Polynomial2D(const eckit::LocalConfiguration config)
  : invars_() {
  // Initialize options_
  options_.deserialize(config);

  // Include required variables
  invars_ += Variable(options_.xVariable.value());
  invars_ += Variable(options_.yVariable.value());
}

// -----------------------------------------------------------------------------

Polynomial2D::~Polynomial2D() {}

// -----------------------------------------------------------------------------

void Polynomial2D::compute(const ObsFilterData & data,
                           ioda::ObsDataVector<float> & out) const {

  // Get nlocs dimension
  size_t nlocs = data.nlocs();

  // Get x and y variable data
  const Variable &xvar = options_.xvar.value();
  ioda::ObsDataVector<float> xvals(data.obsspace(), xvar.toOopsVariables());
  data.get(xvar, xvals);

  const Variable &yvar = options_.yvar.value();
  ioda::ObsDataVector<float> yvals(data.obsspace(), yvar.toOopsVariables());
  data.get(yvar, yvals);

  // Identify filter variables
  Variables filtervars;
  if (options_.filterVariables.value() != boost::none) {
  // read filter variables
    for (const Variable &var : *options_.filterVariables.value())
      filtervars += var;
  } else {
  // if no filter variables explicitly specified, filter out all variables
    filtervars += Variables(data.obsspace().obsvariables());
  }

  // Allocate re-usable placeholder variables
  double c, xi, yi, px, py;
  size_t cvar, ovar;
  std::vector<double> exponents(2);
  LinearFit2DCoefficientsParameters &coeffOpt;
  Polynomial2DTermParameters &term;
  std::vector<double> &coeff;

  // fill in output values
  out.zero()

  for (size_t ivar = 0; ivar < filtervars.size(); ++ivar) {
    // get this filter variable
    Variable fVar = filtervars[ivar];
    ovar = ivar;

    // handle scalar and 1-D Variable cases
    for (size_t jvar = 0; jvar < fVar.size(); ++jvar) {
      if (jvar > 0) {
        ovar = jvar;
      }

      // look for the correct coefficients to match this filter variable
      cvar = -1
      for (size_t cvar_ = 0; cvar_ < options_.polynomialCoefficients.size(); ++cvar_) {
        coeffOpt = options_.polynomialCoefficients.value()[cvar_];
        if (coeffOpt.name.value() == fVar.variable()) {
          if (coeffOpt.channel.value() != boost::none) {
            if (coeffOpt.channel.value() == fVar.channels()[jvar]) {
              cvar = cvar_;
              break;
            }
          } else {
            cvar = cvar_;
            break;
          }
        }
      }
      if (cvar < 0) {
         continue;
      }
      // coeffOpt = options_.polynomialCoefficients.value()[cvar];
      // const std::vector<double> &coeff = coeffOpt.values.value();
      coeff = options_.polynomialCoefficients.value()[cvar].values.value();

      // xvar and yvar indices
      size_t ixvar = std::min(ovar, xvar.size() - 1);
      size_t iyvar = std::min(ovar, yvar.size() - 1);
      // note: use xvals[ixvar][iloc] and yvals[iyvar][iloc] below

      // only loop over available coefficients
      for (size_t iterm = 0; jterm < coeff.size(); ++iterm) {
        c = coeff[iterm]

        term = options.polynomialTerms.value()[iterm];
        px = term.exponents.value()[0];
        py = term.exponents.value()[1];

        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          xi = std::pow(xvals[ixvar][iloc], px);
          yi = std::pow(yvals[iyvar][iloc], py);
          xi *= yi;
          xi *= ci;

          out[ovar][iloc] += (float) xi;
        }
      }
    }
  }
}

// code for producing max cloud_area_fraction_in_atmosphere_layer@GeoVaLs
//  // Get dimensions
//  size_t outVarSize = out.nvars();
//  size_t nlocs = data.nlocs();
//  size_t nlevs = data.nlevs(Variable("cloud_area_fraction_in_atmosphere_layer@GeoVaLs"));
//
//  std::vector<std::vector<float>> CFx(nlocs, std::vector<float>(nlevs));
//  std::vector<float> CFatLevel(nlocs);
//  for (size_t ilev = 0; ilev < nlevs; ++ilev) {
//    const size_t level = nlevs - ilev - 1;
//    data.get(Variable("cloud_area_fraction_in_atmosphere_layer@GeoVaLs"), level, CFatLevel);
//    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
//      CFx[iloc][ilev] = CFatLevel[iloc];
//    }
//  }
//
//  double CFxmax;
//  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
//    // use maximum cloud fraction in column as a proxy for modeled cloud fraction
//    // could alternatively integrate vertically with some weighting function
//    CFxmax = *std::max_element(CFx[iloc].begin(), CFx[iloc].end());
//
//    // same for all variables
//    for (size_t ivar = 0; ivar < outVarSize; ++ivar) {
//      out[ivar][iloc] = CFxmax;
//    }
//  }

// -----------------------------------------------------------------------------

const ufo::Variables & Polynomial2D::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
