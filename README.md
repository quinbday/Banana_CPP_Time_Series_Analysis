# Banana_CPP_Time_Series_Analysis
# Bananas for Bananas: Forecasting Monthly Banana Prices with SARIMA

*CFRM 586 Midterm — Quinlan Day*

## Overview

This project forecasts the monthly U.S. retail price of bananas from January 2026 through January 2027, comparing several time series models and reporting each one's point forecast alongside a 5th/95th percentile interval. The full historical series, sourced from FRED (`APU0000711211`), runs back to 1990.

## Data

- **Source:** FRED monthly banana price index, price per pound (USD).
- **Sample used for modeling:** restricted to 2010–2026 (see *Restricting the Sample* below for why).
- **Missing value:** October 2025 is missing from the raw series and is imputed as 0.667.
- **Series behavior:** bananas are a classic grocery-store loss leader — their price isn't set to maximize profit on the item itself, but to drive store traffic and cross-selling. That pricing logic, rather than a typical supply/demand price process, explains a lot of what shows up in the data.

## Exploratory Analysis

The full 1990–2026 series shows heavy volatility before 2000 that settles into a much calmer pattern afterward, with a level shift around 2008. Breaking the series into bimonthly segments shows a consistent seasonal rhythm — gridlines every two months line up well with local peaks and troughs across most of the sample.

Looking at sub-periods separately:
- **1990–2010:** bimonthly gridlines track peaks and troughs closely for the first two years, then loosen a bit after 1993 as the extreme early volatility fades.
- **2010–2026:** the same bimonthly rhythm holds consistently across the full window, making this the more stable and analyzable stretch of the data.

In short, the raw data is:
- Broken up well by bimonthly segments.
- Extremely volatile before 2000, and doesn't follow the same trend before vs. after 2010.
- Missing October 2025 (imputed at 0.667).
- Unusual for a commodity price in that it's a loss leader and doesn't follow a typical price distribution.
- Seasonal: it rises from January into February, dips to its lowest point in October, and climbs back up from there.

## Restricting the Sample to 2010–2026

Several diagnostics point to dropping the pre-2010 data rather than using the full 36-year history:

- The STL remainder from 1990–2000 shows a large low-frequency wave pattern, and including data before 2010 heavily violates the normality assumption.
- A Breusch-Pagan test on the full series returns a p-value of roughly 8.727 × 10⁻¹³; restricting to 2010–2026 brings that up to 0.0531, fixing the non-constant variance. A log transform doesn't fix it either — its Breusch-Pagan p-value is even smaller, around 7.932 × 10⁻¹⁵.
- Keeping pre-2010 data also violates the 24- and 36-month stationarity assumptions used later in the analysis.
- A Box-Cox search on the full series suggests a transform of λ = −1, which would flip the price scale (large values become small and vice versa) and make the model very hard to interpret. Restricting the search to λ ∈ [0, 1] and restricting the sample to 2010–2026 together resolve this, returning λ = 1 — i.e., no transform needed.

For monthly data, 36 years of history sounds appealing, but here it actively hurts the analysis: the newer, shorter window gives a cleaner, more stationary series to model. Comparing the STL remainder's ACF/PACF for 1990–2026 vs. 2010–2026 confirms this — the shorter window is far easier to work with, since the STL decomposition strips out more of the seasonality in that subset than it does for the full series. Both windows still show a large first-lag spike, a sign that the series still needs differencing; the PACF benefits especially from dropping the early data, showing much smaller spikes elsewhere.

## Stationarity and Differencing

After first-differencing (d = 1) and seasonal differencing (D = 1, lag 12), the series looks much closer to white noise, with only a few unusual lags remaining. The raw (Box-Cox transformed) series fails both the Augmented Dickey-Fuller (ADF) and Kwiatkowski-Phillips-Schmidt-Shin (KPSS) tests:

| Transformation Stage | ADF p-value | ADF Result | KPSS p-value | KPSS Result | Status |
|---|---|---|---|---|---|
| Raw transformed series | 0.8416 | Non-stationary | ≤ 0.0100 | Non-stationary | Diff needed |
| After d = 1 only | ≤ 0.0100 | Stationary | ≥ 0.1000 | Stationary | Stationary |
| After D = 1 only (lag = 12) | 0.0696 | Non-stationary | 0.0547 | Stationary | Diff needed |
| After d = 1 **and** D = 1 | ≤ 0.0100 | Stationary | ≥ 0.1000 | Stationary | Stationary |

Since the series is stationary at both d = 1 alone and d = 1 combined with D = 1, either differencing scheme is usable going forward.

## Seasonal Structure and Model Family Selection

With d = 1 only, the ACF and PACF are clean, with almost no significant lags outside of lag 1 — the seasonal signal isn't visible yet. Once D = 1 is added, a lag-12 spike appears in the ACF, pointing to a seasonal MA term (Q = 1); the PACF shows the same pattern, suggesting the seasonal AR order can be limited to P = 0 or 1.

Running `auto.arima` with only d = 1 lands on an ARIMA(2,1,2), but its AR and MA coefficients come out as near-mirror-image values — a sign of overfitting. Restricting the search so p and q can't both equal 2 instead returns a simple AR(1) process, reinforcing that the seasonal (D = 1) term is doing real work rather than being redundant with a larger non-seasonal order.

## Model Selection

Comparing non-seasonal ARIMA(p,1,q) candidates by RMSE (overall error) and MASE (trend-capture), rather than seasonal versions, gives:

| Model | RMSE | MASE | AICc |
|---|---|---|---|
| ARIMA(1,1,2) | 0.006892 | 0.347523 | −1357.164 |
| **ARIMA(0,1,2)** | 0.006896 | 0.346917 | −1359.024 |
| ARIMA(1,1,1) | 0.006927 | 0.345308 | −1357.376 |
| ARIMA(1,1,0) | 0.006932 | 0.344661 | −1359.126 |
| ARIMA(0,1,1) | 0.006943 | 0.344591 | −1358.547 |
| ARIMA(0,1,0) | 0.007008 | 0.350573 | −1357.012 |

**ARIMA(0,1,2)** is the pick: its RMSE is barely behind the top candidate, its MASE is solidly competitive, and it has the second-best AICc. For a seasonal counterpart, **SARIMA(1,1,2)(0,1,1)** is chosen — it has the second-lowest MASE among the seasonal fits, and the AR(1,1,2) part lines up with the earlier note about the series being well broken up bimonthly, while the (0,1,1) seasonal term captures the lag-12 spike seen in the d=1, D=1 ACF/PACF.

## Benchmarking

Both chosen models are compared against five standard benchmarks — Mean, Naive, Random Walk with Drift, Seasonal Naive, and AR(1):

| Model | RMSE | MASE |
|---|---|---|
| AR(1) | 0.006986 | 0.349474 |
| Mean Model | 0.026276 | 1.453674 |
| Naive Model | 0.007026 | 0.352255 |
| Random Walk w/ Drift | 0.007018 | 0.353954 |
| Seasonal Naive Model | 0.020219 | 1.000000 |
| Auto ARIMA(0,1,2) | 0.006896 | 0.346917 |
| Auto SARIMA(1,1,2)(0,1,1) | 0.006923 | 0.339635 |

Both selected models beat every benchmark on RMSE. On MASE, SARIMA(1,1,2)(0,1,1) comes out on top, followed by AR(1), then ARIMA(0,1,2).

## Conclusion

- **ARIMA(0,1,2)** is the better choice for overall point-forecast accuracy toward the mean price level. Since it carries no seasonal component, its forecast is a flat point estimate month to month.
- **SARIMA(1,1,2)(0,1,1)** better captures the monthly seasonal shape of the forecast, at a small cost to overall accuracy.

## Forecast Results

12-month-ahead forecasts (point estimate with 5th/95th percentile bounds), starting from January:

**Mean, Naive, Drift, and Seasonal Naive benchmarks**

| Month | Mean – Mean | Mean – Low | Mean – High | Naive – Mean | Naive – Low | Naive – High | Drift – Mean | Drift – Low | Drift – High | SNaive – Mean | SNaive – Low | SNaive – High |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| January | 0.598 | 0.546 | 0.650 | 0.653 | 0.639 | 0.667 | 0.653 | 0.640 | 0.667 | 0.619 | 0.579 | 0.659 |
| February | 0.598 | 0.546 | 0.650 | 0.653 | 0.634 | 0.672 | 0.654 | 0.634 | 0.673 | 0.625 | 0.585 | 0.665 |
| March | 0.598 | 0.546 | 0.650 | 0.653 | 0.629 | 0.677 | 0.654 | 0.630 | 0.678 | 0.635 | 0.595 | 0.675 |
| April | 0.598 | 0.546 | 0.650 | 0.653 | 0.625 | 0.681 | 0.654 | 0.627 | 0.682 | 0.655 | 0.615 | 0.695 |
| May | 0.598 | 0.546 | 0.650 | 0.653 | 0.622 | 0.684 | 0.655 | 0.624 | 0.686 | 0.654 | 0.614 | 0.694 |
| June | 0.598 | 0.546 | 0.650 | 0.653 | 0.619 | 0.687 | 0.655 | 0.621 | 0.689 | 0.657 | 0.617 | 0.697 |
| July | 0.598 | 0.546 | 0.650 | 0.653 | 0.617 | 0.689 | 0.655 | 0.618 | 0.693 | 0.666 | 0.626 | 0.706 |
| August | 0.598 | 0.546 | 0.650 | 0.653 | 0.614 | 0.692 | 0.656 | 0.616 | 0.696 | 0.670 | 0.630 | 0.710 |
| September | 0.598 | 0.546 | 0.650 | 0.653 | 0.612 | 0.694 | 0.656 | 0.614 | 0.698 | 0.667 | 0.627 | 0.707 |
| October | 0.598 | 0.546 | 0.650 | 0.653 | 0.609 | 0.697 | 0.656 | 0.612 | 0.701 | 0.664 | 0.624 | 0.704 |
| November | 0.598 | 0.546 | 0.650 | 0.653 | 0.607 | 0.699 | 0.657 | 0.610 | 0.704 | 0.656 | 0.616 | 0.696 |
| December | 0.598 | 0.546 | 0.650 | 0.653 | 0.605 | 0.701 | 0.657 | 0.608 | 0.706 | 0.653 | 0.613 | 0.693 |

**AR(1), ARIMA(0,1,2), and SARIMA(1,1,2)(0,1,1)**

| Month | AR(1) – Mean | AR(1) – Low | AR(1) – High | ARIMA – Mean | ARIMA – Low | ARIMA – High | SARIMA – Mean | SARIMA – Low | SARIMA – High |
|---|---|---|---|---|---|---|---|---|---|
| January | 0.651 | 0.635 | 0.668 | 0.652 | 0.639 | 0.666 | 0.655 | 0.641 | 0.669 |
| February | 0.649 | 0.627 | 0.673 | 0.652 | 0.634 | 0.670 | 0.657 | 0.639 | 0.675 |
| March | 0.648 | 0.621 | 0.676 | 0.652 | 0.630 | 0.674 | 0.657 | 0.634 | 0.679 |
| April | 0.646 | 0.617 | 0.678 | 0.652 | 0.626 | 0.678 | 0.658 | 0.632 | 0.684 |
| May | 0.644 | 0.612 | 0.680 | 0.652 | 0.623 | 0.681 | 0.658 | 0.628 | 0.687 |
| June | 0.643 | 0.609 | 0.681 | 0.652 | 0.620 | 0.684 | 0.656 | 0.624 | 0.688 |
| July | 0.641 | 0.605 | 0.682 | 0.652 | 0.617 | 0.687 | 0.655 | 0.621 | 0.690 |
| August | 0.640 | 0.602 | 0.683 | 0.652 | 0.614 | 0.689 | 0.656 | 0.619 | 0.693 |
| September | 0.639 | 0.599 | 0.683 | 0.652 | 0.612 | 0.692 | 0.655 | 0.615 | 0.694 |
| October | 0.637 | 0.597 | 0.684 | 0.652 | 0.610 | 0.694 | 0.656 | 0.614 | 0.698 |
| November | 0.636 | 0.595 | 0.684 | 0.652 | 0.608 | 0.696 | 0.655 | 0.611 | 0.699 |
| December | 0.635 | 0.592 | 0.684 | 0.652 | 0.606 | 0.698 | 0.659 | 0.613 | 0.705 |
