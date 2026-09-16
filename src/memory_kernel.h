#ifndef REMSTATS_MEMORY_KERNEL_H
#define REMSTATS_MEMORY_KERNEL_H

#include <RcppArmadillo.h>
#include <algorithm>

/* CustomKernel

memory = "custom": the weight of a past event is an arbitrary, user-supplied
function of the elapsed time (lag) since the event, tabulated on a grid. The R
side (validate_memory) flattens the table into

    memory_value = c(lag_1, ..., lag_n, weight_1, ..., weight_n)

with the lags sorted increasingly. The weight of a lag is the tabulated weight
at the grid point closest to it (ties go to the smaller lag); lags outside the
grid take the weight of the nearest end point. This mirrors nearest-neighbour
lookup with findInterval() in R, so an R implementation on the same grid gives
identical results.

Unlike the exponential kernel of memory = "decay", a general kernel is not
multiplicative in elapsed time, so the statistics cannot be updated
recursively: every time point has to re-weight the whole relevant past.
*/
struct CustomKernel
{
  arma::vec lags;
  arma::vec wts;

  CustomKernel() {}

  explicit CustomKernel(const arma::vec &memory_value)
  {
    const arma::uword n2 = memory_value.n_elem;
    if (n2 < 2 || (n2 % 2) != 0)
    {
      Rcpp::stop("memory = 'custom' requires memory_value = c(lags, weights) of even length >= 2.");
    }
    const arma::uword n = n2 / 2;
    lags = memory_value.head(n);
    wts = memory_value.tail(n);
    if (!lags.is_sorted("ascend"))
    {
      Rcpp::stop("memory = 'custom': the lag grid in memory_value must be sorted increasingly.");
    }
    if (!lags.is_finite() || !wts.is_finite())
    {
      Rcpp::stop("memory = 'custom': lags and weights in memory_value must be finite.");
    }
  }

  bool empty() const { return lags.n_elem == 0; }

  // Weight at the grid point nearest to 'lag'
  double operator()(double lag) const
  {
    const arma::uword n = lags.n_elem;
    const double *first = lags.memptr();
    const double *last = first + n;
    const double *pos = std::lower_bound(first, last, lag);
    const arma::uword idx = static_cast<arma::uword>(pos - first);
    if (idx == 0)
    {
      return wts(0);
    }
    if (idx == n)
    {
      return wts(n - 1);
    }
    if (lags(idx) == lag)
    {
      return wts(idx);
    }
    // lags(idx - 1) < lag < lags(idx): nearest wins, ties go left
    const double d_left = lag - lags(idx - 1);
    const double d_right = lags(idx) - lag;
    return (d_right < d_left) ? wts(idx) : wts(idx - 1);
  }
};

#endif
