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
  oops::Log::trace() << "Polynomial2D::compute, retrieving xvar: " << xvar << std::endl;
  // TODO(JJG): would be nice to make this work for filter-variable resolved variables (like ObsValue)
  // ioda::ObsDataVector<float> xvals(data.obsspace(), xvar.toOopsVariables());
  // ioda::ObsDataVector<float> xvals(data.obsspace(), xvar.toOopsVariables(), xvar.group());
  // std::vector<float> xvals(nlocs);
  std::vector<float> xvalsf;
  data.get(xvar, xvalsf);

  const Variable &yvar = options_.yvar.value();
  oops::Log::trace() << "Polynomial2D::compute, retrieving yvar: " << yvar << std::endl;
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

  oops::Log::trace() << "Polynomial2D::compute, initialize outVars and filterVars" << std::endl;

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
  // Variables outVars(out.varnames());


  // oops::Log::trace() << "Polynomial2D::compute, outVars: " << outVars << std::endl;
  oops::Log::trace() << "Polynomial2D::compute, filterVars: " << filterVars << std::endl;

  // Allocate re-usable placeholder variables
  const std::vector<Polynomial2DTermParameters>
    &allTerms = options_.polynomialTerms.value();
  const std::vector<LinearFitCoefficientsParameters>
    &allVarCoeffs = options_.polynomialCoefficients.value();

  size_t cvar, ovar, fvar;
  double xi, yi, px, py;

  // fill in output values
  // outVars.size() and filterVars.size() are always 1
  // oops::Log::trace() << "Polynomial2D::compute, outVars.size(): " << outVars.size() << std::endl;
  oops::Log::trace() << "Polynomial2D::compute, filterVars.size(): " << filterVars.size() << std::endl;

  // ASSERT(outVars.size() == 1);
  ASSERT(filterVars.size() == 1);

  fvar = 0;
  // for (size_t fvar = 0; fvar < outVars.size(); ++fvar) {
  // for (size_t fvar = 0; fvar < filterVars.size(); ++fvar) {
  // get filter variable
  Variable filterVar = filterVars[fvar];
  oops::Log::trace() << "Polynomial2D::compute, filterVar: " << filterVar << std::endl;
  oops::Log::trace() << "Polynomial2D::compute, filterVar.size(): " << filterVar.size() << std::endl; 
  oops::Log::trace() << "Polynomial2D::compute, filterVar.variable(): " << filterVar.variable() << std::endl;
  oops::Log::trace() << "Polynomial2D::compute, filterVar.channels(): " << filterVar.channels() << std::endl;
  oops::Log::trace() << "Polynomial2D::compute, filterVar.channels().size(): " << filterVar.channels().size() << std::endl;

  // get output variable
  // Variable outVar = outVars[fvar];
  // oops::Log::trace() << "Polynomial2D::compute, outVar: " << outVar << std::endl;
  // oops::Log::trace() << "Polynomial2D::compute, outVar.size(): " << outVar.size() << std::endl; 
  // oops::Log::trace() << "Polynomial2D::compute, outVar.variable(): " << outVar.variable() << std::endl;
  // oops::Log::trace() << "Polynomial2D::compute, outVar.channels(): " << outVar.channels() << std::endl;
  // oops::Log::trace() << "Polynomial2D::compute, outVar.channels().size(): " << outVar.channels().size() << std::endl;

  oops::Log::trace() << "Polynomial2D::compute, out.nvars(): " << out.nvars() << std::endl;
  oops::Log::trace() << "Polynomial2D::compute, out.nlocs(): " << out.nlocs() << std::endl;
  oops::Log::trace() << "Polynomial2D::compute, data.nlocs(): " << data.nlocs() << std::endl;


  for (size_t ovar = 0; ovar < nvars; ++ovar) {
    // oops::Log::trace() << "Polynomial2D::compute, outVar.channels()[" << ovar << "]: " << outVar.channels()[ovar] << std::endl; //duplicated below
  }

  for (size_t ovar = 0; ovar < nvars; ++ovar) {
    oops::Log::trace() << "Polynomial2D::compute, ovar: " << ovar << std::endl;
    // initialize to missing
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      oops::Log::trace() << "Polynomial2D::compute, missing, iloc: " << iloc << std::endl;
      out[ovar][iloc] = missing;
    }
  }
  size_t nVarCoeffs = allVarCoeffs.size();
  // handle scalar and 1-D Variable cases (i.e., channels)
  // for (size_t ovar = 0; ovar < outVar.size(); ++ovar) {
  for (size_t ovar = 0; ovar < nvars; ++ovar) {
  // for (size_t ovar = 0; ovar < filterVar.size(); ++ovar) {
    oops::Log::trace() << "Polynomial2D::compute, ovar: " << ovar << std::endl;

    // oops::Log::trace() << "Polynomial2D::compute, outVar.variable(ovar): " << outVar.variable(ovar) << std::endl;
    // oops::Log::trace() << "Polynomial2D::compute, filterVar.variable(ovar): " << filterVar.variable(ovar) << std::endl;
    // oops::Log::trace() << "Polynomial2D::compute, outVar.channels()[ovar]: " << outVar.channels()[ovar] << std::endl; //duplicated below

    // look for the correct coefficients to match this filter variable and channel (if applicable)
    cvar = nVarCoeffs;
    for (size_t csrch = 0; csrch < nVarCoeffs; ++csrch) {
      const std::string nm = allVarCoeffs[csrch].name.value();
      oops::Log::trace() << "Polynomial2D::compute, allVarCoeffs[csrch].name.value(): " << nm << std::endl;

      if (allVarCoeffs[csrch].name.value() == filterVar.variable()) {
        const int ch = allVarCoeffs[csrch].channel.value().get();
        oops::Log::trace() << "Polynomial2D::compute, allVarCoeffs[csrch].channel.value(): " << ch << std::endl;

        if (allVarCoeffs[csrch].channel.value() != boost::none) {
          // ASSERT(outVar.channels().size() > ovar);
          // if (allVarCoeffs[csrch].channel.value().get() == outVar.channels()[ovar]) {
            // oops::Log::trace() << "Polynomial2D::compute, outVar.channels()[" << ovar <<"]: " << outVar.channels()[ovar] << std::endl;

          ASSERT(filterVar.channels().size() > ovar);
          if (allVarCoeffs[csrch].channel.value() == filterVar.channels()[ovar]) {
            oops::Log::trace() << "Polynomial2D::compute, filterVar.channels()[ovar]: " << filterVar.channels()[ovar] << std::endl;
            cvar = csrch;
            break;
          }
        } else {
          cvar = csrch;
          break;
        }
      }
    }

    oops::Log::trace() << "Polynomial2D::compute, cvar: " << cvar << std::endl;

    // skip filter variables for which there are no coefficients
    if (cvar < 0 || cvar >= nVarCoeffs) {
      // initialize to missing
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        oops::Log::trace() << "Polynomial2D::compute, missing, iloc: " << iloc << std::endl;
        out[ovar][iloc] = missing;
      }
      oops::Log::trace() << "Polynomial2D::compute, out[ovar] initialized to missing" << std::endl;
      continue;
    }

    // initialize coefficients
    oops::Log::trace() << "Polynomial2D::compute, nVarCoeffs: " << nVarCoeffs << std::endl;

    const std::vector<double> coeff = allVarCoeffs[cvar].values.value();
    size_t nTerms = coeff.size();

    oops::Log::trace() << "Polynomial2D::compute, nTerms:" << nTerms << std::endl;

    // xvar and yvar variable indices
    // size_t ixvar = std::min(ovar, xvar.size() - 1);
    // size_t iyvar = std::min(ovar, yvar.size() - 1);

    // allocate intermediate output variables
    std::vector<std::vector<double>> z(nlocs, std::vector<double>(nTerms));

    // evaluate polynomial
    // only use as many terms as there are available coefficients
    for (size_t iterm = 0; iterm < nTerms; ++iterm) {
      // oops::Log::trace() << "Polynomial2D::compute, iterm: " << iterm << std::endl;
      px = static_cast<double>(allTerms[iterm].exponents.value()[0]);
      py = static_cast<double>(allTerms[iterm].exponents.value()[1]);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        // xi = std::pow(xvals[ixvar][iloc], px);
        xi = std::pow(xvals[iloc], px);

        // yi = std::pow(yvals[iyvar][iloc], py);
        yi = std::pow(yvals[iloc], py);

        z[iloc][iterm] = coeff[iterm];
        z[iloc][iterm] *= xi;
        z[iloc][iterm] *= yi;
      }
    }

    oops::Log::trace() << "Polynomial2D::compute, accumulating z" << std::endl;
    // oops::Log::trace() << "Polynomial2D::compute, nlocs: " << nlocs << std::endl;

//    std::vector<double> zabs(nTerms);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      // oops::Log::trace() << "Polynomial2D::compute, accumulate, iloc: " << iloc << std::endl;
      // oops::Log::debug() << "Polynomial2D::compute, ioxyz"
      //   << ", " << iloc
      //   << ", " << ovar
      //   << ", " << xvals[iloc]
      //   << ", " << yvals[iloc]
      //   << ", " << z[iloc]
      //   << std::endl;

//      for (size_t iterm = 0; iterm < nTerms; ++iterm) {
//        // oops::Log::debug() << "Polynomial2D::compute, iterm: " << iterm << std::endl;
//        zabs[iterm] = std::abs(z[iloc][iterm]);
//      }

      // re-order z from smallest to largest magnitude before sum
      // oops::Log::debug() << "Polynomial2D::compute, sort" << std::endl;

//      std::sort(z[iloc].begin(), 
//                z[iloc].end(),
//                [&zabs](const size_t &i1, const size_t &i2) {return zabs[i1] < zabs[i2];});

      // oops::Log::debug() << "Polynomial2D::compute, accumulate" << std::endl;

      out[ovar][iloc] = static_cast<float>(
        std::accumulate(z[iloc].begin(), z[iloc].end(), 0.0));

      // oops::Log::debug() << "Polynomial2D::compute, out"
      //   << ", " << out[ovar][iloc]
      //   << std::endl;

    }
    oops::Log::trace() << "Polynomial2D::compute, completed ovar: " << ovar << std::endl;

  }
  // }
  oops::Log::trace() << "Polynomial2D::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & Polynomial2D::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
