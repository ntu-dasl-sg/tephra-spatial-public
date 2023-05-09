# Analysis procedure and notes

Here we outline our analysis.

1. **Prepare data**
* The `data` folder contains the raw observation data.
* We removed the following outliers:
    * Dataset 1: 629785, 9122495     (Load = 7) -- Very small value surrounded by large values
    * Dataset 2: 645089.9   9126336   (Load =  3.50e+02  -- produces large errors
    * Dataset 2: 641221.3   9121778   (Load =  2.10e+02  -- produces large errors
    * Dataset 2: 647377.0   9126551   (Load =   2.24e+02 -- produces large errors

2.  **Prepare input files for inversion**
* We modified the Tephra2 source code (`src/tephra2-source-code/`) to accomodate more cost functions
* We compiled the modified source codes to generate executables for each cost function (`src/executables-costfunctions/`)
* We prepared working folders for the inversion (`results/input for Tephra2/`). These contain the configuration files, wind data, observation data (outliers removed), and the inversion script to run the inversion. 

3.  **Run inversion**
* We use the computing cluster ("Komodo") in Earth Observatory of Singapore to run the inversion.
* The best-fit modelled outputs are saved to `results/best models/`. This contains the values of the best-fit source parameters, and a gridded file containing predicted values at each grid.
* Output csv files containing the predictions at each sampled site are saved to `results/output csv/`. 
* Output figures saved to `results/graphics-maps/` and `graphics-pred_vs_obs`

4. **Check suitability of cost function using goodness of fit tests**
* We provide a sample code of how common goodness-of-fit tests can be done in R. See `stat-tests/statistical-tests.html`

5.  **Create train and test set for analysis`**
* Generate a train and test set from Dataset 1 and 2. The test set is randomly selected from Dataset 1 only with 20-80 test train ratio (resulting to 16 test points). All of Dataset 2 goes to the training set.

6. **Implement the unweighted and weighted kriging techniques**
* Use `src\fusion.R`.
* The results will be in the form of raster maps (.RData file)

7. **Compare the performance of the kriging techniques using LOOCV**
* Use `src\loocv.R`.
