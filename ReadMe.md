
Does ESG Add Value? ESG vs. S&P 500: Risk Modeling with ARMA-GARCH, Copulas & Monte Carlo (2015-2025)
Applied econometrics project (Université Sorbonne Paris Nord, Master Big Data & Business Intelligence, 2026). Does applying ESG filters to a US equity portfolio change its return or its risk profile? We compare an ESG ETF with the S&P 500 over ten years of daily data, from descriptive performance ratios to volatility modeling, tail-dependence analysis and 10,000-scenario Monte Carlo simulations.

Authors: Baptiste Comby, Myriam Aoun, Fanny Lattard

Key findings
ESG (SUSA)	Market (SPY)
Annualized return	13.06%	13.61%
Annualized volatility	17.93%	17.82%
Max drawdown	-32.93%	-33.72%
Sharpe ratio	0.564	0.598
Sortino ratio	0.069	0.071
No ESG alpha: after controlling for market, technology and energy exposure, the intercept is not significant (p = 0.45). ESG does not outperform once risk is accounted for.
A "low-carbon" bias, not a tech bias: the tech factor (XLK) is not significant (p = 0.18), but the energy factor (XLE) is significant and negative (-0.014, p < 0.001). The ESG portfolio differs from the market mainly because of what it excludes (oil & gas).
No safe haven in crashes: a copula analysis gives a lower tail dependence of ≈ 0.89. When the market crashes, the ESG portfolio almost always crashes with it (β ≈ 0.99).
A conditional "shock absorber": in the 10,000 simulated one-year scenarios, ESG has a slightly thinner left tail and a more contained volatility response under stress, with a comparable median outcome.
Stable market relationship: a Rec-CUSUM test (p = 0.419) shows no structural break in ESG's market sensitivity, even through COVID-19.
Takeaway: ESG integration does not create excess return. It mostly reshapes risk, through sector exclusions, with a modest cushioning effect in extreme scenarios.

Data
Assets: iShares MSCI USA ESG Select ETF (SUSA) vs. SPDR S&P 500 ETF (SPY)
Factors and controls: XLK (technology), XLE (energy), ^VIX (volatility index), ^TNX (10-year Treasury yield, used as risk-free rate)
Source: Yahoo Finance, downloaded automatically with quantmod (no data file needed)
Period: January 2015 to 2025, daily adjusted prices; only days when all markets are open are kept
Returns: log-returns for econometric models, simple returns for performance metrics
Methodology
1. Stylized facts and model justification

ADF test (stationarity)
Breusch-Pagan test (heteroscedasticity, which justifies a GARCH model)
Breusch-Godfrey test (autocorrelation, which justifies an ARMA model)
Jarque-Bera test (fat tails, which justifies a Skewed Student-t distribution)
Diagnostic dashboard: return series, density vs. normal, Q-Q plot, ACF of squared returns
2. Performance and resilience

Wealth index, daily returns, drawdown profile
Sharpe, Sortino and Treynor ratios, downside capture
Rolling 12-month volatility, Sharpe ratio and beta; underwater plot
3. Volatility modeling

auto.arima constrained to a maximum order of 2 for parsimony
Final model: ARMA(1,1)-GARCH(1,1) with Skewed Student-t innovations, fitted with the rugarch hybrid solver
Validation: standardized residuals Q-Q plot and Ljung-Box test
4. Monte Carlo simulation

10,000 scenarios × 252 trading days for each portfolio
Outputs: fan charts (50% and 90% intervals), distribution of final prices with the crash zone (losses beyond 20%), forecast volatility, spaghetti plot
5. Dependence and attribution

Copula family selected automatically with VineCopula; lower tail dependence
Factor regression ESG ~ Market + Tech + Energy, with Newey-West (HAC) robust standard errors
6. Structural stability

Rec-CUSUM test with a COVID-19 marker
Repository structure
.
├── README.md
├── esg_vs_spy_analysis.R     # full analysis (steps 1 to 6)
└── paper/
    └── esg_vs_spy_report.pdf # written report (in French, 30 pages)
How to run
Install R (4.x) and RStudio.
Install the packages:
r
install.packages(c("quantmod", "PerformanceAnalytics", "rugarch", "tseries", "copula",
                   "VineCopula", "ggplot2", "reshape2", "forecast", "FinTS", "lmtest",
                   "car", "sandwich", "strucchange", "gt", "dplyr", "scales", "gridExtra"))
Open esg_vs_spy_analysis.R and run it. The data is downloaded from Yahoo Finance automatically.
Note: the 10,000-scenario simulations take a few minutes. Results may differ slightly from the report if you extend the end date, since the data is downloaded live.

Tools
R · quantmod · PerformanceAnalytics · rugarch · VineCopula / copula · forecast · lmtest / sandwich (HAC) · strucchange · ggplot2

Academic project. The written report is in French; the code comments are partly in French.


