Anyone can add more cost functions to the original source code by editing the lines indicated below:

Edit the configuration files

- `tephra2-inversion.conf` - Add option in line 35
- `inversionConfig.conf` - line 81

Edit the run-inversion script:

- `runInversion_PBS.sh` - line 3 and all instances of loss function name

Edit the source code:

- `inversion_src/fit_tests.c` — Add formula
- `inversion_src/minimizing_func.c` - Line ~170+
- `common_src/parameters.h` - Line ~89
- `common_src/prototypes.h` - Line 30

Compile the source code using  `make clean` and `make`. Create a new name for the executable.

Incorporate the new name of the executable in the following files:
- `inversion_src/makefile`
- `forward_src/makefile`
- `makefile`

