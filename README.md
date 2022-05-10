# TIFX04-LEMA
Simulation and optimization of a linear electromagnetic accelerator (LEMA) for an applied physics bachelor's thesis project.

## FEM-analysis of accelerator stage
The folder _ComsolAnalysis_ contains code for curve fitting analytical functions to
data from the finite element modeling software COMSOL.
The data was generated from a model (_LEMA_model.mph_) of an accelerator stage.
This is a way to summarise the complex geometry of in a few analytical functions of one or two variables.

## Time dependent simulation

The folder _LEMA_ contains code for running a time dependent simulation of projectile movement using parameters from the COMSOL analysis.

The simulation process is summarised in a single file _one_file_summary.py_. A more modular version of the code can be found in _stage.py_, _projectile.py_ and _experiment.py_. The remaining files are runnable scripts for generating various results form the simulation. In particular, _designplot.py_ generates a plot used for building an optimised accelerator.