## Overview

Here we demonstrate the following goodness-of-fit tests using sample results from tephra inversion modelling.
Open `statistical-tests.html` to see a sample code in R. We also provide sample inputs in this folder.
The html was knitted using `statistical-tests.Rmd`. 

- Kolmogorov-Smirnov (K-S) test
- Cram'er-von Mises (CvM) test
- Anderson-Darling (AD) test
- Shapiro-Wilk test (for Gaussian distributions only)
- Skewness

Each loss function is associated with distributional assumptions on the residuals, hence we can check the appropriateness of the choice of loss function based on the goodness-of-fit of the residuals to these assumed distributions.

## Loading the sample residuals

For this code, we need sample results from different cost functions. Here we utilise results from the following:

- Mean square error (MSE)
- Chi-square error
- Mean square logarithmic error (MSLE)
- Mean absolute error (MAE)
- Mean absolute percentage error (MAPE)

The sample results are generated from multiple inversions, each using a different cost function. Using the best-fit parameters, we predict the tephra load at the observation sites. The results from each inversion are saved in separate tables.

The tables consist of three columns:

- `Observation`: Value of observed data 
- `Prediction`: Value of modelled output
- `Residual`: Calculated as Observation - Prediction

CSV files for these tables are provided in this folder. You can load them in R using:

```
mse <- read.csv("sample_results-mse.csv") # table for Mean square error (MSE)
chi <- read.csv("sample_results-chisquare.csv") # table for Chi-square error
msl <- read.csv("sample_results-msl.csv") # table for Mean square logarithmic error (MSLE)
mae <- read.csv("sample_results-mae.csv") # table for Mean absolute error (MAE)
map <- read.csv("sample_results-map.csv") # table for Mean absolute percentage error (MAPE)
```