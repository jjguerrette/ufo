/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/Polynomial2D.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

static ObsFunctionMaker<Polynomial2D> makerPoly2D_("Polynomial2D");

// -----------------------------------------------------------------------------

Polynomial2D::Polynomial2D(const eckit::LocalConfiguration config)
  : invars_() {
  // Initialize options_
  options_.deserialize(config);

  // Include required variables
  invars_ += options_.xvar.value();
  invars_ += options_.yvar.value();
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

  // Identify filter variables (same behavior as FilterBase::FilterBase)
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
  const std::vector<Polynomial2DTermParameters>
    &allTerms = options_.polynomialTerms.value();
  const std::vector<LinearFitCoefficientsParameters>
    &allCoeff = options_.polynomialCoefficients.value();

  double c, xi, yi, zi, px, py;
  size_t cvar, ovar;

  // initialize output to zero
  // should this be missing instead?
  out.zero();

  // fill in output values
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
      cvar = -1;
      for (size_t cvar_ = 0; cvar_ < allCoeff.size(); ++cvar_) {
        if (allCoeff[cvar].name.value() == fVar.variable()) {
          if (allCoeff[cvar].channel.value() != boost::none) {
            if (allCoeff[cvar].channel.value() == fVar.channels()[jvar]) {
              cvar = cvar_;
              break;
            }
          } else {
            cvar = cvar_;
            break;
          }
        }
      }
      // skip filter variables for which there are no coefficients
      if (cvar < 0) {
         continue;
      }
      const std::vector<double> &coeff = allCoeff[cvar].values.value();

      // xvar and yvar indices
      size_t ixvar = std::min(ovar, xvar.size() - 1);
      size_t iyvar = std::min(ovar, yvar.size() - 1);

      // allocate intermediate output variables
      std::vector<double> z(coeff.size());
      std::vector<double> zabs(coeff.size());

      // evaluate polynomial
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        // only use as many terms as there are available coefficients
        for (size_t iterm = 0; iterm < coeff.size(); ++iterm) {
          px = allTerms[iterm].exponents.value()[0];
          xi = std::pow(xvals[ixvar][iloc], px);

          py = allTerms[iterm].exponents.value()[1];
          yi = std::pow(yvals[iyvar][iloc], py);

          z[iterm] = coeff[iterm];
          z[iterm] *= xi;
          z[iterm] *= yi;
          zabs[iterm] = std::abs(z[iterm]);
        }

        // re-order z from smallest to largest magnitude before sum
        std::sort(z.begin(), z.end(), [&zabs](size_t i1, size_t i2) {return zabs[i1] < zabs[i2];});
        out[ovar][iloc] = std::accumulate(z.begin(), z.end(), 0.0f);
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & Polynomial2D::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
