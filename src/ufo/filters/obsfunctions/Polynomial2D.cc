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
#include "oops/util/missingValues.h"
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
  oops::Log::trace() << "Polynomial2D::compute started" << std::endl;
  const float missing = util::missingValue(missing);

  // Get nlocs dimension
  size_t nlocs = data.nlocs();

  // Get x and y variable data
  const Variable &xvar = options_.xvar.value();
  oops::Log::trace() << "Polynomial2D::compute, retrieving xvar: " << xvar << std::endl;
  // TODO(JJG): would be nice to make this work for filter-variable resolved variables (like ObsValue)
  // ioda::ObsDataVector<float> xvals(data.obsspace(), xvar.toOopsVariables());
  // ioda::ObsDataVector<float> xvals(data.obsspace(), xvar.toOopsVariables(), xvar.group());
  // std::vector<float> xvals(nlocs);
  std::vector<float> xvals;
  data.get(xvar, xvals);

  const Variable &yvar = options_.yvar.value();
  oops::Log::trace() << "Polynomial2D::compute, retrieving yvar: " << yvar << std::endl;
  // TODO(JJG): would be nice to make this work for filter-variable resolved variables (like ObsValue)
  // ioda::ObsDataVector<float> yvals(data.obsspace(), yvar.toOopsVariables());
  // ioda::ObsDataVector<float> yvals(data.obsspace(), yvar.toOopsVariables(), yvar.group());
  // std::vector<float> yvals(nlocs);
  std::vector<float> yvals;
  data.get(yvar, yvals);

  oops::Log::trace() << "Polynomial2D::compute, initialize outVars and filterVars" << std::endl;

  // Identify filter variables (same behavior as FilterBase::FilterBase)
  Variables outVars(out.varnames());
  Variables filterVars;
  if (options_.filterVariables.value() != boost::none) {
    // read filter variables
    for (const Variable &var : *options_.filterVariables.value())
      filterVars += var;
  } else {
    // if no filter variables explicitly specified, filter out all variables
    filterVars += Variables(data.obsspace().obsvariables());
  }

  oops::Log::trace() << "Polynomial2D::compute, outVars: " << outVars << std::endl;
  oops::Log::trace() << "Polynomial2D::compute, filterVars: " << filterVars << std::endl;

  // Allocate re-usable placeholder variables
  const std::vector<Polynomial2DTermParameters>
    &allTerms = options_.polynomialTerms.value();
  const std::vector<LinearFitCoefficientsParameters>
    &allVarCoeffs = options_.polynomialCoefficients.value();

  size_t cvar, ovar;
  int px, py;
  double xi, yi;

  // fill in output values
  // outVars.size() and filterVars.size() are always 1
  // for (size_t ivar = 0; ivar < outVars.size(); ++ivar) {
  for (size_t ivar = 0; ivar < filterVars.size(); ++ivar) {
    // get output variable and filter variable
    // Variable outVar = outVars[ivar];
    Variable filterVar = filterVars[ivar];

    oops::Log::debug() << "Polynomial2D::compute, filterVar.variable: " << filterVar.variable() << std::endl;

    // for (size_t jvar = 0; jvar < outVar.size(); ++jvar) {
    for (size_t jvar = 0; jvar < filterVar.size(); ++jvar) {
      // handle scalar and 1-D Variable cases (i.e., channels)
      ovar = jvar;

      // initialize to missing
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        out[ovar][iloc] = missing;
      }
      // oops::Log::debug() << "Polynomial2D::compute, ovar: " << ovar << std::endl;

      // look for the correct coefficients to match this filter variable
      cvar = allVarCoeffs.size();
      for (size_t csrch = 0; csrch < allVarCoeffs.size(); ++csrch) {
        if (allVarCoeffs[csrch].name.value() == filterVar.variable()) {
          if (allVarCoeffs[csrch].channel.value() != boost::none) {
            // if (allVarCoeffs[csrch].channel.value() == outVar.channels()[jvar]) {
            if (allVarCoeffs[csrch].channel.value() == filterVar.channels()[jvar]) {
              // oops::Log::debug() << "Polynomial2D::compute, filterVar.channels()[jvar]: " << filterVar.channels()[jvar] << std::endl;
              cvar = csrch;
              break;
            }
          } else {
            cvar = csrch;
            break;
          }
        }
      }

      oops::Log::debug() << "Polynomial2D::compute, cvar: " << cvar << std::endl;

      // skip filter variables for which there are no coefficients
      if (cvar >= allVarCoeffs.size()) {
        continue;
      }

      // initialize coefficients
      // oops::Log::debug() << "Polynomial2D::compute, allVarCoeffs.size(): " << allVarCoeffs.size() << std::endl;

      const std::vector<double> &coeff = allVarCoeffs[cvar].values.value();
      size_t nCoeff = coeff.size();

      // oops::Log::debug() << "Polynomial2D::compute, nCoeff:" << nCoeff << std::endl;

      // xvar and yvar variable indices
      // size_t ixvar = std::min(ovar, xvar.size() - 1);
      // size_t iyvar = std::min(ovar, yvar.size() - 1);

      // allocate intermediate output variables
      std::vector<std::vector<double>> z(nlocs, std::vector<double>(nCoeff));

      // evaluate polynomial
      // only use as many terms as there are available coefficients
      for (size_t icof = 0; icof < nCoeff; ++icof) {
        // oops::Log::debug() << "Polynomial2D::compute, icof: " << icof << std::endl;
        px = allTerms[icof].exponents.value()[0];
        py = allTerms[icof].exponents.value()[1];
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          // xi = std::pow(xvals[ixvar][iloc], px);
          xi = std::pow(xvals[iloc], px);

          // yi = std::pow(yvals[iyvar][iloc], py);
          yi = std::pow(yvals[iloc], py);

          z[iloc][icof] = coeff[icof];
          z[iloc][icof] *= xi;
          z[iloc][icof] *= yi;
        }
      }

      // oops::Log::debug() << "Polynomial2D::compute, accumulating z" << std::endl;
      // oops::Log::debug() << "Polynomial2D::compute, nlocs: " << nlocs << std::endl;

      std::vector<double> zabs(nCoeff);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        // oops::Log::debug() << "Polynomial2D::compute, iloc: " << iloc << std::endl;
        for (size_t iterm = 0; iterm < nCoeff; ++iterm) {
          // oops::Log::debug() << "Polynomial2D::compute, iterm: " << iterm << std::endl;
          zabs[iterm] = std::abs(z[iloc][iterm]);
        }

        // re-order z from smallest to largest magnitude before sum
        // oops::Log::debug() << "Polynomial2D::compute, sort" << std::endl;

        std::sort(z[iloc].begin(), 
                  z[iloc].end(),
                  [&zabs](const size_t &i1, const size_t &i2) {return zabs[i1] < zabs[i2];});

        // oops::Log::debug() << "Polynomial2D::compute, accumulate" << std::endl;

        out[ovar][iloc] = static_cast<float>(
          std::accumulate(z[iloc].begin(), z[iloc].end(), 0.0));
      }
      oops::Log::debug() << "Polynomial2D::compute, completed ovar: " << ovar << std::endl;

    }
  }
  oops::Log::trace() << "Polynomial2D::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & Polynomial2D::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
