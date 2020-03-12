## Requirements
* Python >= 3.6
    * Required modules
        * numpy (tested with version 1.18.1)
        * scipy (tested with version 1.4.1)
    * Optional modules
        * matplotlib (only needed if plotting the solution, tested with version 3.1.3)
        * numba (reduces setup time substantially, tested with version 0.48.0)
        * tqdm (provides interactive progress bar, tested with version 4.42.1)

In order to build the C and C++ libraries for solving the membrane equations and the extracellular system, respectively, there are some additional dependencies. The program falls back to slower solvers if the libraries haven't been built. Note that the library wrapping ViennaCL must be built in order to use the AMG preconditioner.
The libraries require
* OS: macOS or Linux
* CMake >= 3.9
* Compilers for C and C++ (both GCC and Clang were tested)

## Building the C and C++ libraries

### Downloading ViennaCL
In order to build the C++ library for solving the linear systems with CG and AMG, ViennaCL must be downloaded.
From `c/external`, run
```
./get_viennacl.sh
```

This will download ViennaCL and place it under `c/external/viennacl`.

### Configuring and building
From `c`, run
```
./rebuild_lib.sh
```

This will invoke CMake and GNU Make to build the libraries.


## Usage
The main program is called `run_simulation.py`. For a list of available options, run
```
python run_simulation.py --help
```
