## Results

All the results for all cost functions, and weighting schemes are saved here in the `results` folder. To generate these results, we used batch inversion run using the shell command: `qsub -t 9-11 runInversion_PBS.sh`.

- `best models` - the Tephra2 output (`tephra.out and tephra2.out`) for the run with the best-fit ESPs
- `graphics-maps` - tephra load maps
- `graphics-pred_vs_obs` - plots of predicted vs observed, a typical way of showing the results of inversion
- `output csv` - contains the residuals at each data point, and a summary of training and test errors, and best-fit ESPs
- `kriging_rasters` - Output of the model-data fusion in raster format


## Input files used in the paper:

* `input for Tephra2/unweighted (ds1 and ds2)` - Inversion using points from Dataset 1 and 2. No weights were applied to the training set.
* `input for Tephra2/unweighted (ds1 only)` - Inversion using points from Dataset 1 only. No weights were applied to the training set.
* `input for Tephra2/weighted (ds1 and ds2)` - Inversion using points from Dataset 1 and 2. Weights were applied to the training set.
