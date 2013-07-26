Optimizer optimization
======================

Source code used in "Optimal optimization in comparative methods" by Kurt Michels, Jeremy Beaulieu, Barb Banbury, Naim Matasci and Brian O'Meara.

The script main.R allows the reproduction of the results presented in the article.

fitContinuous.modified.R and fitDiscrete.modified.R are based on fitContinuous and fitDiscrete from the package geiger v. 1.0 by Luke Harmon, respectively.
ace.modified.R is based on ace from the package ape v. 3.0-3 by Emmanuel Paradis
A patched version of the package ouch by Aaron A. King is also included and will be installed in a teporary location.

Dependencies
------------
The following packages are required and are available through CRAN

ape
data.table
mvtnorm
optimx
