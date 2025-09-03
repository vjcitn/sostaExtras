# sostaExtras

Define a simple app to overview some spatialExperiment examples from [sosta](https://bioconductor.org/packages/sosta)  package

## Basic view

use `example(surveyDamond, ask=FALSE)` to start the app after successful installation
via `BiocManager::install("vjcitn/sostaExtras")`

Change the initial 'target' selection to 'INS'

![run1](man/figures/soste1.jpg)

## Cell type view tab

![run2](man/figures/soste2.jpg)

## A different target, filtered by relative 'abundance' level

![run3](man/figures/soste3.jpg)

Note that the qthr control can be used to 'decolor' cells that
have target abundance measure below the selected quartile of
the target
