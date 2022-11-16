Refer to [S Biasse's guide](https://e5k.github.io/codes/utilities/2018/06/06/inversion/) for more details and troubleshooting on Tephra2 inversion.
Tephra2 source code available [here](https://github.com/ljc-geo/tephra2).

## Input files used in the paper:
* `input for Tephra2/unweighted (ds1 and ds2)` - Inversion using points from Dataset 1 and 2. No weights were applied to the training set.
* `input for Tephra2/unweighted (ds1 only)` - Inversion using points from Dataset 1 only. No weights were applied to the training set.
* `input for Tephra2/weighted (ds1 and ds2)` - Inversion using points from Dataset 1 and 2. Weights were applied to the training set.

## Accommodating new cost functions:
For Tephra2 inversion, we only needed the two executables created by the compilation of Tephra2 (tephra2012_inversion and tephra2-2012). For our paper, we modified these executables and saved them in `input for Tephra2/executables`.
When using a specific cost function, a user can use the corresponding executable based on the filename we provided. (e.g. `tephra2-Chiweighted` and `tephraChiweighted_inversion` is used when inverting using the Chi-square cost function.

Anyone can add more cost functions to the original source code by editing the lines indicated below:

- `Tephra2` → `conf_files` → `tephra2-inversion.conf` - Add option in line 35
- `inversion_src` → `fit_tests.c` — Add formula
- `inversion_src` → `minimizing_func.c` - Line ~170+
- `common_src` → `parameters.h` - Line ~89
- `common_src` → `prototypes.h` - Line 30

Compile the source code using  `make clean` and `make`. Create a new name for the executable.
Incorporate the new name of the executable in the following files:
- `/Tephra2_unwt_chi/inversion_src/makefile`
- `/Tephra2_unwt_chi/forward_src/makefile`
- `/Tephra2_unwt_chi/makefile`
- `inversionConfig.conf` - line 81
- `runInversion_PBS.sh` - line 3 and all instances of loss function name

## Running the inversion

We used the batch inversion run using the shell command:
`qsub -t 9-11 runInversion_PBS.sh`

## Results

All the results for all cost functions, and weighting schemes are saved in the `results` folder.
- `best models` - the Tephra2 output (`tephra.out and tephra2.out`) for the run with the best-fit ESPs
- `graphics-maps` - tephra load maps
- `graphics-pred_vs_obs` - plots of predicted vs observed, a typical way of showing the results of inversion
- `output csv` - contains the residuals at each data point, and a summary of training and test errors, and best-fit ESPs
