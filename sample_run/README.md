# Analysis procedure and notes

Here we demonstrate how we post-process the results of Tephra2 inversion using R.

1.  **Run inversion**
* We use the computing cluster ("Komodo") in Earth Observatory of Singapore to run the inversion.
* Input working folders and executables are saved in this repository in: `inversion/input for Tephra2`.
* Multiple working folders can be post-processed at the same time.
* Transfer the inversion results that you want to post-process in a folder called `ResultsfromKomodo`

2.  **Create train and test set with `get_test.rmd`**
* This file generates a train and test set from Dataset 1 and 2.
* The test set is randomly selected from Dataset 1 only with 20-80 test train ratio (resulting to 16 test points)
* All of Dataset 2 goes to the training set.
* Raw data comes from `Data` folder.
* Output csv files for Komodo runs saved to `Train and Test` folder. Output figures saved to `Graphics`
* Outliers removed:
    * Dataset 1: 629785, 9122495     (Load = 7) -- Very small value surrounded by large values
    * Dataset 2: 645089.9	9126336	  (Load =  3.50e+02	 -- produces large errors
    * Dataset 2: 641221.3	9121778	  (Load =  2.10e+02	 -- produces large errors
    * Dataset 2: 647377.0	9126551	  (Load =   2.24e+02 -- produces large errors

3.  **Get best models using `save_best_fit.Rmd`**
* This automatically picks out the folder with the best fit value and copies the output files to a folder called `Best Models`.

4. **Generate plots and csv with `process_bestfit.Rmd`**
* This creates maps of the best forward model and a csv of train errors, test errors and best-fit ESPs.
