# Analysis procedure

## 1. **Create train and test set with `get_test.rmd`**
* This file generates a train and test set from Dataset 1 and 2.
* The test set is randomly selected from Dataset 1 with 20-80 test train ratio (resulting to 16 test points)
* All of Dataset 2 goes to the training set.
* Raw data comes from `Data` folder.
* Output csv files for Komodo runs saved to `Train and Test` folder. Output figures saved to `Graphics`
* Outliers removed:
    * Dataset 1: 629785, 9122495     (Load = 7) -- Very small value surrounded by large values
    * Dataset 2: 645089.9	9126336	  (Load =  3.50e+02	 -- produces large errors
    * Dataset 2: 641221.3	9121778	  (Load =  2.10e+02	 -- produces large errors
    * Dataset 2: 647377.0	9126551	  (Load =   2.24e+02 -- produces large errors

## 2.  **Run input files in Komodo**
* Use executables from `_Executables_16testpoints` in `t2_inversion` folder in Komodo.

## 3.  **Postprocess output files using `Matlab_batch/batch_postprocess.m`
* This generates `summary.txt` that contains the fit values of each discretized mass-height values.

## 4.  **Get best models using `postprocess_komodo.Rmd`**
* This automatically picks out the folder with the best fit value and copies the output files to a folder called `Best Models`.

## 5. **Generate plots and csv with `process_bestfit.Rmd`**
* This creates maps of the best forward model and a csv of train errors, test errors and ESPs
