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

#include "eckit/exception/Exceptions.h"

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
  size_t nlocs = out.nlocs();
  size_t nvars = out.nvars();

  // Get x and y variable data
  const Variable &xvar = options_.xvar.value();
  // oops::Log::debug() << "Polynomial2D::compute, retrieving xvar: " << xvar << std::endl;
  // TODO(JJG): would be nice to make this work for filter-variable resolved variables (like ObsValue)
  // ioda::ObsDataVector<float> xvals(data.obsspace(), xvar.toOopsVariables());
  // ioda::ObsDataVector<float> xvals(data.obsspace(), xvar.toOopsVariables(), xvar.group());
  // std::vector<float> xvals(nlocs);
  std::vector<float> xvalsf;
  data.get(xvar, xvalsf);

  const Variable &yvar = options_.yvar.value();
  // oops::Log::debug() << "Polynomial2D::compute, retrieving yvar: " << yvar << std::endl;
  // TODO(JJG): would be nice to make this work for filter-variable resolved variables (like ObsValue)
  // ioda::ObsDataVector<float> yvals(data.obsspace(), yvar.toOopsVariables());
  // ioda::ObsDataVector<float> yvals(data.obsspace(), yvar.toOopsVariables(), yvar.group());
  // std::vector<float> yvals(nlocs);
  std::vector<float> yvalsf;
  data.get(yvar, yvalsf);

  std::vector<double> xvals(nlocs), yvals(nlocs);
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    xvals[iloc] = static_cast<double>(xvalsf[iloc]);
    yvals[iloc] = static_cast<double>(yvalsf[iloc]);
  }

  // oops::Log::debug() << "Polynomial2D::compute, initialize outVars and filterVars" << std::endl;

  // Identify filter variables (same behavior as FilterBase::FilterBase)
  Variables filterVars;
  if (options_.filterVariables.value() != boost::none) {
    // read filter variables
    for (const Variable &var : *options_.filterVariables.value())
      filterVars += var;
  } else {
    // if no filter variables explicitly specified, filter out all variables
    filterVars += Variables(data.obsspace().obsvariables());
  }

  // Establish output variables
  // requires changes to channels handling in Variable::toOopsVariables and Variables::toOopsVariables
  Variables outVars(out.varnames());


  // oops::Log::debug() << "Polynomial2D::compute, outVars: " << outVars << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, filterVars: " << filterVars << std::endl;

  // Retrieve exponents for all available polynomial terms
  const std::vector<Polynomial2DTermParameters>
    &allTerms = options_.polynomialTerms.value();
  size_t nTerms = allTerms.size();
  std::vector<double> px(nTerms), py(nTerms);
  for (size_t iterm = 0; iterm < nTerms; ++iterm) {
    px[iterm] = static_cast<double>(allTerms[iterm].exponents.value()[0]);
    py[iterm] = static_cast<double>(allTerms[iterm].exponents.value()[1]);
  }

  // Retrieve coefficient parameters
  const std::vector<LinearFitCoefficientsParameters>
    &allVarCoeffs = options_.polynomialCoefficients.value();
  size_t nVarCoeffs = allVarCoeffs.size();

  // fill in output values
  // outVars.size() and filterVars.size() are assumed to be 1, because Perform Action filters
  // require their 'error function' to be a scalar Variable, which may optionally apply to
  // multiple channels, but not multiple non-channel filter variables
  ASSERT(outVars.size() == 1);
  ASSERT(filterVars.size() == 1);
  size_t fvar = 0;

  // get filter and output variables
  Variable filterVar = filterVars[fvar];
  // oops::Log::debug() << "Polynomial2D::compute, filterVar: " << filterVar << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, filterVar.size(): " << filterVar.size() << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, filterVar.variable(): " << filterVar.variable() << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, filterVar.channels(): " << filterVar.channels() << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, filterVar.channels().size(): " << filterVar.channels().size() << std::endl;
  Variable outVar = outVars[fvar];
  // oops::Log::debug() << "Polynomial2D::compute, outVar: " << outVar << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, outVar.size(): " << outVar.size() << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, outVar.variable(): " << outVar.variable() << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, outVar.channels(): " << outVar.channels() << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, outVar.channels().size(): " << outVar.channels().size() << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, out.nvars(): " << out.nvars() << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, out.nlocs(): " << out.nlocs() << std::endl;
  // oops::Log::debug() << "Polynomial2D::compute, data.nlocs(): " << data.nlocs() << std::endl;

  // initialize all output values to missing
  for (size_t ovar = 0; ovar < nvars; ++ovar) {
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[ovar][iloc] = missing;
    }
  }

  // Allocate re-usable placeholder variables
  // handle scalar and 1-D Variable cases (i.e., channels)
  size_t cvar;
  for (size_t ovar = 0; ovar < nvars; ++ovar) {
    // look for the correct coefficients to match this filter variable and channel (if applicable)
    cvar = nVarCoeffs;
    for (size_t csrch = 0; csrch < nVarCoeffs; ++csrch) {
      if (allVarCoeffs[csrch].name.value() == filterVar.variable()) {
        if (allVarCoeffs[csrch].channel.value() != boost::none) {
          ASSERT(outVar.channels().size() > ovar);
          if (allVarCoeffs[csrch].channel.value().get() == outVar.channels()[ovar]) {
            cvar = csrch;
            break;
          }
        } else {
          cvar = csrch;
          break;
        }
      }
    }

    // skip filter variables for which there are no coefficients
    if (cvar < 0 || cvar >= nVarCoeffs) {
      continue;
    }

    // initialize coefficients
    const std::vector<double> &coeff = allVarCoeffs[cvar].values.value();
    size_t nCoeff = coeff.size();

    // xvar and yvar variable indices
    // size_t ixvar = std::min(ovar, xvar.size() - 1);
    // size_t iyvar = std::min(ovar, yvar.size() - 1);

    // allocate intermediate output variables
    std::vector<double> z(nCoeff); // zabs(nCoeff);
    // std::vector<std::pair<double, double>> zp(nCoeff);
    // double z_;

    // evaluate polynomial
    // only use as many terms as there are available coefficients
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      for (size_t iterm = 0; iterm < nCoeff; ++iterm) {
        z[iterm] = coeff[iterm];
        // z[iterm] *= std::pow(xvals[ixvar][iloc], px[iterm]);
        // z[iterm] *= std::pow(yvals[iyvar][iloc], py[iterm]);
        z[iterm] *= std::pow(xvals[iloc], px[iterm]);
        z[iterm] *= std::pow(yvals[iloc], py[iterm]);
        // zabs[iterm] = std::abs(z[iterm]);

        // z_ = coeff[iterm];
        // // z_ *= std::pow(xvals[ixvar][iloc], px[iterm]);
        // // z_ *= std::pow(yvals[iyvar][iloc], py[iterm]);
        // z_ *= std::pow(xvals[iloc], px[iterm]);
        // z_ *= std::pow(yvals[iloc], py[iterm]);
        // zp[iterm] = std::make_pair(z_, std::abs(z_))
      }
      // TODO: figure out why this causes an error
      // std::sort(z.begin(),
                // z.end(),
                // [&zabs](size_t &i1, size_t &i2) {return zabs[i1] < zabs[i2];});

      // TODO: try this alternative
      // std::sort(zp.begin(),
                // zp.end(),
                // [](const auto &a, const auto &b) {return a.second < b.second;});
      // for (size_t iterm = 0; iterm < nCoeff; ++iterm) {
        // z[iterm] = zp[iterm].first;
      // }

      out[ovar][iloc] = static_cast<float>(
        std::accumulate(z.begin(), z.end(), 0.0));
    }
    oops::Log::debug() << "Polynomial2D::compute, completed ovar: " << ovar << std::endl;
  }
  oops::Log::trace() << "Polynomial2D::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & Polynomial2D::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
