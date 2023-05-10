# Inversion and forward estimation with process-based models: an investigation into cost functions, uncertainty-based weights and model-data fusion

This is a code and data repository for our manuscript. Further inquiries can be directed to the corresponding author: Maricar Rabonza (MARICARL001@e.ntu.edu.sg).

## Citation:
```
Rabonza ML, Nguyen M, Biass S, Jenkins S, Taisne B, and Lallemant D (2023) Inversion and forward estimation with process-based models: an investigation into cost functions, uncertainty-based weights and model-data fusion. EMS. In-prep
```

## REPOSITORY OUTLINE

This repository is organised as:

```
tephra-spatial-public/
├── README.md              # overview + Tephra2 quickstart instructions + cost function overview
│ 
├── data/                  # data files used in the project
│   └── README.md          # describes where data came from
│ 
├── src/                   # codes in the project
│   ├── executables/       # modified Tephra2 source code used in the study + compiled executables
│   ├── tephra2-source     # unmodified Tephra2 source code from developers
│   └── fusion/            # R codes used for the model-data fusion and LOOCV
│ 
├── results/               # results of the analysis (data, tables, figures)
│ 
├── doc/                   # documentation
│   ├── stat-tests/        # short report on how we conducted goodness-of-fit tests in R
│   ├── add-costf.md       # instructions to add more cost functions in the Tephra2 source code
│   └── komodo.md          # outline of our analysis using the cluster in EOS Singapore
│ 
├── manuscript/            # files for the journal manuscript submitted to EMS
│   └── README.md          # details of the paper
└── LICENSE                # License
```

## TEPHRA2 OVERVIEW

Our study utilises tephra2 (VERSION 2.0), a tephra dispersion simulation tool used to estimate the mass of tephra that will accumulate at a site or over a region, given explosive eruption conditions. The code was developed by Laura J. Connor and Costanza Bonadonna. Detailed documentation of the software are provided in [Github](https://github.com/geoscience-community-codes/tephra2), [Tephra2 website](https://gscommunitycodes.usf.edu/geoscicommunitycodes/public/tephra2/tephra2.php), and Sebastien Biass' [website](https://e5k.github.io/codes/utilities/2018/06/06/inversion/). 

- Year first available: 2013
- Version used in our study: 2.0, Updated 01-27-2018
- Program language: C
- License: GNU General Public License v3.0
- Availability: https://github.com/geoscience-community-codes/tephra2
- Program size: 963 KB
- Hardware requirements: The inversion model requires MPI (Message Passing Interface) libraries and should be run on a computing cluster with multiple compute nodes. 

## TEPHRA2 QUICK START

**The instructions below are directly derived from the Tephra2 [Github](https://github.com/geoscience-community-codes/tephra2).**

First, download the Tephra2 source code from [Github](https://github.com/geoscience-community-codes/tephra2). A copy of the source code is also provided here in [`/src/tephra2-source/`](https://github.com/ntu-dasl-sg/tephra-spatial-public/tree/main/src/tephra2-source). If interested, instructions for parallelisation are provided [here](https://e5k.github.io/codes/utilities/2018/06/06/inversion/).

### Install the following dependencies:

The following dependencies can be installed via terminal. Visit the corresponding download websites for installation instructions applicable to your computer.

- [gcc](https://gcc.gnu.org)
- [openmpi](https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html)
- [gc](https://www.linuxfromscratch.org/blfs/view/svn/general/gc.html), [gc-devel](https://yum-info.contradodigital.com/view-package/base/gc-devel/), [libatomic_ops](https://github.com/ivmai/libatomic_ops) (for opensuse linux) or [libgc-dev](https://howtoinstall.co/en/libgc-dev), [libgc1c2](https://howtoinstall.co/en/libgc1c2) (for ubuntu linux) or [bdw-gc](https://brewinstall.org/install-bdw-gc-on-mac-with-brew/) (mac - homebrew)

### Compilation

tephra2 is written in C, a compiled language, and must be compiled before it can be executed. See [README](https://github.com/geoscience-community-codes/tephra2/blob/master/README.usage) for quickstart instructions. Specifically, to compile the linux executables, in the top level directory, on the command line type:

```
make
```

This will compile both the forward model and the inversion model. If you do not have openmpi installed the inversion model will fail to compile, but the forward model will still be compiled. This is OK if you only want the forward model.

### Usage

Two configuration files are needed, one for the forward model (tephra2.conf) and another for the inversion model (tephra2-inversion.conf). The program user can change any of the values inside the configuration file; the keywords must not be changed. The forward model can be run on a single computing node. The inversion model requires a computing cluster with multiple compute nodes. 

To run the forward model, at the command line, type:
```
tephra2_2020 tephra2.conf grid_file wind_file > tephra2.out
```
where,
- **tephra2_2020** is the name of the executable
- **tephra2.conf** is the name of the file of configuration parameters (an example)
- **grid_file** is a text file of 3 columns separated by spaces (an example) following the format:
`Easting(m)  Northing(m)  Grid-Elevation(m)`
- **wind_file** is a 3-column text file of wind data (an example) following the format:
`Height(masl)  Wind-Speed(m/s)  Wind-Direction(wind vector azimuth in degrees)`. [Here is a guide](https://github.com/geoscience-community-codes/tephra2/blob/master/plotting_scripts/readme.wind) to download and create wind files for Tephra2 based on NOAA reanalysis data.
- **tephra2.out** is the output file name where the tephra accumulation values will be written, following the format:
`Easting(m)  Northing(m)  Grid-Elevation(m)  Mass(kg/m^2) [weight percent of modeled phi fractions]`

To run the inversion model, at the command line, type:
```
mpirun -np nodes -hostfile machines tephra2-inversion_2020 tepha2-inversion.conf data_file wind_file
```

where,
- **mpirun** is the wrapper script for gcc when using mpi libraries
- **nodes** is the number of cluster compute nodes to use
- **machines** is a text file listing the name of each compute node and the number of cpu cores useable on that node (an example)
- **tephra2-inversion_2020** is the executable name
- **tephra2-inversion.conf** is the file of parameters (an example)
- **data_file** is a text file of tephra accumulation data from chosen area (an example) following the format:
Easting(m)  Northing(m)  Grid-elevation(m)  Mass(kg/m^2)
- **wind_file** is a text file of wind data to use for the inversion (same format as above)

## USING OTHER COST FUNCTIONS WITH TEPHRA2

For our paper, we modified the Tephra2 executables to change the cost functions in the optimisation. The modified executables for the following cost functions are provided in: [`/src/executables/`](https://github.com/ntu-dasl-sg/tephra-spatial-public/tree/main/src/executables).
- Chi-squared error
- Mean absolute error (MAE)
- Mean squared error (MSE)
- Mean absolute percentage error (MAPE)
- Mean squared logarithmic error (MSLE)
- Tokyo-log error

Proceed with running the inversion using the executable that corresponds to the cost function of interest.

To incorporate other cost functions to the original source code, see our guide: [`doc/add-costf.md`](https://github.com/ntu-dasl-sg/tephra-spatial-public/blob/main/doc/add-costf.md)


----

## Author Contributions

MR, MN and DL designed the research. MR wrote the manuscript. SB provided guidance in conducting tephra inversion modelling. All authors discussed the results and commented on the manuscript.

## Funding

This research is supported by the Earth Observatory of Singapore, the National Research Foundation Singapore, and the Singapore Ministry of Education under the NRF-NRFF2018-06 award. 

## Acknowledgements

We are grateful to Edwin Tan for his support in the use of the High-Performance Computing Cluster, Komodo, in the Earth Observatory of Singapore. We thank George Williams for his support and insights related to the tephra load datasets.
