# ✈️ Marseille Airport Traffic Forecasting

This project analyses and forecasts **monthly passenger traffic** at Marseille airport using classical time series techniques.  
The goal is to model long-term trends and strong yearly seasonality to produce reliable short-term forecasts.

---

## 🔍 Overview

- Monthly traffic data from **1982 to 2016**
- Clear **upward trend** and strong **seasonal pattern** (s = 12)
- Step-by-step transformation to achieve stationarity:
  - log transformation (stabilize variance)
  - first difference (remove trend)
  - seasonal difference (remove seasonality)
- Fitting and comparing several **SARIMA / SARIMAX models**
- Final models evaluated on a **hold-out validation set**

---

## 🧠 Method Highlights

- Visual inspection of the series (trend, variance, seasonality)
- ACF / PACF analysis to guide SARIMA order selection
- Estimation of:
  - full SARIMA model with seasonal components
  - simplified SARIMA model with fewer parameters
  - SARIMAX model with **dummy variables** to correct outliers
- Diagnostic checks:
  - residual autocorrelation (Ljung–Box test)
  - normality (Shapiro, Jarque–Bera)
  - constant variance (McLeod–Li test)
- Comparison of models using:
  - AIC
  - residual diagnostics
  - coverage of confidence intervals
  - forecast accuracy on a test period

---

## 📊 Results (Summary)

- Two models provide good forecasts:
  - A **SARIMAX model** with dummy variables correcting extreme outliers
  - A more parsimonious **SARIMA (1,1,0)(0,1,1)[12]** model on a stable sub-period (2001–2016)
- One model offers better **fit quality** (AIC, interval coverage),  
  while the other satisfies more classical assumptions (normal and homoscedastic residuals).
- Forecasts capture both the overall trend and the seasonal pattern of passenger traffic.

---

## 🛠 Technologies

- **R**
- `forecast`, `TSA`, `tseries`, `zoo`


