# Inversion and forward estimation with process-based models: an investigation into cost functions, uncertainty-based weights and model-data fusion

This is a code and data repository for our manuscript. Further inquiries can be directed to the corresponding author: Maricar Rabonza (MARICARL001@e.ntu.edu.sg).

## Citation:
```
Rabonza ML, Nguyen M, Biasse S, Jenkins S, Taisne B, and Lallemant D (2023) Inversion and forward estimation with process-based models: an investigation into cost functions, uncertainty-based weights and model-data fusion. EMS. In-prep
```

## REPOSITORY OUTLINE

This repository is organised as:

```
tephra-spatial-public/
├── README.md              # overview + Tephra2 quickstart instructions + cost function overview
├── data/                  # data files used in the project
│   └── README.md          # describes where data came from
├── src/                   # codes in the project
│   ├── executables/       # modified Tephra2 source code used in the study + compiled executables
│   ├── tephra2-source     # unmodified Tephra2 source code from developers
│   └── fusion/            # R codes used for the model-data fusion and LOOCV
├── results/               # results of the analysis (data, tables, figures)
├── doc/                   # documentation
│   ├── stat-tests/        # short report on how we conducted goodness-of-fit tests in R
│   ├── add-costf.md       # instructions to add more cost functions in the Tephra2 source code
│   └── komodo.md          # outline of our analysis using the cluster in EOS Singapore
├── manuscript/            # files for the journal manuscript submitted to EMS
│   └── README.md          # details of the paper
└── LICENSE                # License
```

## TEPHRA2 OVERVIEW

Our study utilises tephra2 (VERSION 2.0), a tephra dispersion simulation tool used to estimate the mass of tephra that will accumulate at a site or over a region, given explosive eruption conditions. The code was developed by Laura J. Connor and Costanza Bonadonna. Detailed documentation of the software are provided in [Github](https://github.com/geoscience-community-codes/tephra2), [Tephra2 website](https://gscommunitycodes.usf.edu/geoscicommunitycodes/public/tephra2/tephra2.php), and Sebastian Biass' [website](https://e5k.github.io/codes/utilities/2018/06/06/inversion/). **The instructions below are directly derived from these sources.**

- Year first available: 2013
- Version used in our study: 2.0, Updated 01-27-2018
- Program language: C
- License: GNU General Public License v3.0
- Availability: https://github.com/geoscience-community-codes/tephra2
- Program size: 963 KB
- Hardware requirements: The inversion model requires MPI (Message Passing Interface) libraries and should be run on a computing cluster with multiple compute nodes. 

### Tephra2 Quick Start

First, download the Tephra2 source code, which is available on [Github](https://github.com/geoscience-community-codes/tephra2). Alternatively, a copy of the source code is provided here in `/src/tephra2-source/`. If interested, instructions for parallelisation are provided [here](https://e5k.github.io/codes/utilities/2018/06/06/inversion/).

#### Install dependencies before compiling:

- gcc
- openmpi
- gc, gc-devel, libatomic_ops (for opensuse linux) or libgc-dev, libgc1c2 (for ubuntu linux) or bdw-gc (mac - homebrew)

#### Compilation

tephra2 is written in C, a compiled language, and must be compiled before it can be executed. The compilation of Tephra2 results to two executables (tephra2012_inversion and tephra2-2012).
```
    To compile type make
```

Two configuration files are needed, one for the forward model (tephra2.conf) and another for the inversion model (tephra2-inversion.conf). The program user can change any of the values inside the configuration file; the keywords must not be changed.

```
    Edit tephra2.conf 
    This is the forward model configuration file.
    
    Edit tephra2-inversion.conf
    This is the inversion configuration file.
    
    Edit the run-inversion script. 
    This script will:
    1)  run the inversion part of tephra2, 
    2)  run PERL/GMT scripts which generate some PNG images of wind field, plume probabilities,
    3)  run the forward tephra2 model with the best-fit results,
    4)  run PERL/GMT scripts to plot contours
``` 

#### Usage

The forward model can be run on a single computing node.
To run the forward model, at the command line, type:

```
tephra2_2020 tephra2.conf grid_file wind_file > tephra2.out
```
where,

- **tephra2_2020** is the name of the executable
- **tephra2.conf** is the name of the file of configuration parameters (an example)
- **grid_file** is is a text file of 3 columns separated by spaces (an example) following the format:
`Easting(m)  Northing(m)  Grid-Elevation(m)`
- **wind_file** is a 3-column text file of wind data (an example) following the format:
`Height(masl)  Wind-Speed(m/s)  Wind-Direction(wind vector azimuth in degrees)`
- **tephra2.out** is the output file name where the tephra accumulation values will be written, following the format:
`Easting(m)  Northing(m)  Grid-Elevation(m)  Mass(kg/m^2) [weight percent of modeled phi fractions]`

To run the inversion and forward model together, run the `run-inversion` script. This requires a computing cluster with multiple compute nodes. 
    
```
    USAGE: nohup sh run-inversion > nohup.out &
    USE: cd tail -f nohup.out
    to view stderr and stdout as it is being saved to this file.
```

Or to run the inversion model, at the command line, type:
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


### Changing the cost functions

For our paper, we modified the Tephra2 executables to change the cost functions in the optimisation. The modified executables for the following cost functions are provided in: `/src/executables/`
- Chi-squared error
- Mean absolute error (MAE)
- Mean squared error (MSE)
- Mean absolute percentage error (MAPE)
- Mean squared logarithmic error (MSLE)
- Tokyo-log error

Proceed with running the inversion using the executable that corresponds to the cost function of interest.

To incorporate other cost functions to the original source code, see our guide: `doc/add-costf.md`


----

## Author Contributions

MR lead the writing, research and analysis of both case studies. All authors contributed to the conceptualisation and design of the study. DL conceived of the idea of celebrating successes in disaster risk reduction using counterfactual analysis. YL and DL provided critical feedback that shaped the research, analysis and manuscript.

## Funding

This project is supported by the National Research Foundation, Prime Minister’s Office, Singapore under the NRF-NRFF2018-06 award, the Earth Observatory of Singapore, the National Research Foundation of Singapore, and the Singapore Ministry of Education under the Research Centers of Excellence initiative. MR is supported by a PhD scholarship from the Earth Observatory of Singapore. 

## Acknowledgements

We are grateful to Edwin Tan for his support in the use of the High-Performance Computing Cluster, Komodo, in the Earth Observatory of Singapore. We thank George Williams for his support and insights related to the tephra load datasets.
