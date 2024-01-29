## A Non-Dimensional Scaling Approach to Modeling Material Pyrolysis

This repository contains the data, scripts, and technical documents required to generate the pre-print for this paper submitted to the 2024 Combustion Symposium.

## Instructions
* Clone this repository with submodules
  ```
  git clone --recurse-submodules https://github.com/johodges/nondimensional_spyro
  ```
* Pre-process data
  ```
  cd nondimensional_spyro/scripts
  python fsri_collect_thermophysical_properties.py
  python fsri_collect_cone_data.py
  python process_fsri_database.py
  python process_faa_data.py
  python process_fpl_data.py
  cd ../..
  ```
* Build FDS with updated spyro model
  git clone git@github.com:johodges/fds.git && cd fds
  git checkout nondim_spyro && git checkout fb838f80229163
  cd Build/impi_intel_linux
  ./make_fds.sh
  ln -s fds_impi_intel_linux fds
  cd ../..
* Generate plots
  ```
  python evaluate_faa_polymers.py
  python evaluate_statistics.py
  ```
* Build paper
  ```
  pdflatex hodges2024_pyrolysis.tex
  bibtex hodges2024_pyrolysis
  pdflatex hodges2024_pyrolysis.tex
  pdflatex hodges2024_pyrolysis.tex
  ```
