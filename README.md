## A Non-Dimensional Scaling Approach to Modeling Material Pyrolysis

This repository contains the data, scripts, and technical documents required to generate the SPyro cone validation cases.

## Instructions
* Clone this repository with submodules
  ```
  git clone --recurse-submodules https://github.com/johodges/materials
  ```
* Pre-process data
  ```
  cd materials/scripts
  python initialize.py
  cd ../..
  ```
* Build latest FDS (example for impi_intel_linux target, but use relevant make target for platform)
  git clone git@github.com:firemodels/fds.git && cd fds
  cd Build/impi_intel_linux
  ./make_fds.sh
  ln -s fds_impi_intel_linux fds
  cd ../..
* Run cases
  ```
  cd materials/scripts
  python evaluate_statistics_new.py
  ```
