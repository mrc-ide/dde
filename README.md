# dde

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Travis-CI Build Status](https://travis-ci.org/richfitz/dde.svg?branch=master)](https://travis-ci.org/richfitz/dde)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/richfitz/dde?branch=master&svg=true)](https://ci.appveyor.com/project/richfitz/dde)
[![codecov.io](https://codecov.io/github/richfitz/dde/coverage.svg?branch=master)](https://codecov.io/github/richfitz/dde?branch=master)

This package solves ordinary differential equations (ODEs), delay differential equations (DDEs) and discrete-time *difference* (or recursion) equations, perhaps involving delays.

For all the solvers, the target function can be an R function or a compiled function.

## Ordinary and delay differential equations

ODEs are solved with `DOPRI5` and `DOP853`; the 5th order and 8th order Dormand Prince methods [Dormand Prince](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method), based off of [Ernst Hairer's Fortran implementations](http://www.unige.ch/~hairer/software.html), but implemented in C.  These integrators support dense output.

To solve DDEs, I use the same approach as Hairer in his `RETARD` algorithm, exploiting the dense output of the Dormand Prince solvers.

This is an alternative approach to fitting delay DDE models to using deSolve.  It exists because I have had problems fitting very large DDE systems in deSolve, possibly because the order of interpolation is lower than the order of integration which can cause problems with the overall accuracy of the solution.

By using the dense output, the solution can be computed at any time point the solver has passed to the same accuracy as the solution itself.  This sidesteps the interpolation problem at the cost of a bit more book-keeping.

To store the history without using ever-growing (or just huge) amounts of memory, `dde` uses a [ring buffer](https://github.com/richfitz/ring) to hold the history over time.  This means that the memory required to store the solution does not grow as the total integration length increases (though you still need to pick an amount of memory that scales with the maximum number of steps that span your longest lag at any point in the integration).

These solvers are suitable only for nonstiff problems.

The interface is fairly different to the deSolve interface.  A `deSolve` compatible interface may be provided later (see [this issue](https://github.com/richfitz/dde/issues/2)).

## Discrete time models

Solving discrete time equations is much simpler; you don't have much choice but just to iterate the model.  The package implements this efficiently for compiled models, and also allows models to reference their previous history.
