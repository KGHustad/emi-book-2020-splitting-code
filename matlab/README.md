## Requirements
  - One of the following
    - [MATLAB](https://www.mathworks.com/products/matlab.html) (__recommended__) (2018a and 2019b have been tested)
    - [Octave](https://www.gnu.org/software/octave/) (Version 5.1.0 has been tested.)
      - __Note that Octave must be built with support for 64-bit integers in order to run the ischemia simulation with 22x12 cells from the chapter.__
      - Some of the plotting commands may need to be adjusted because of differences between MATLAB and Octave.

## Usage
### Minimal example
We provide a small example program `run_simulation_mini.m` for verifying that the simulation code runs on the system. This program performs a single time step for a mesh with 8x4 cells, where the middle 2x2 cells are ischemic. This minimal example should run in about 30 seconds on a powerful laptop.

### Main programs
There are two main programs:

`run_simulation.m` runs the simulation and stores the result in a file called `solution.mat`. This should be run on a computer with at least 32 GB of memory.

`make_plot.m` loads the result from `solution.mat`, plots it, and stores the result to the file `ischemia.png`.
