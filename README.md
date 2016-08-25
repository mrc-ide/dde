# dde

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Travis-CI Build Status](https://travis-ci.org/richfitz/dde.svg?branch=master)](https://travis-ci.org/richfitz/dde)

An R package that implements the `DOPRI5` solver of [Ernst Hairer](http://www.unige.ch/~hairer/software.html), along with the delay differential equation (DDE) form (originally called `RETARD`).

This is an alternative approach to fitting delay DDE models to using deSolve.  It exists because I have had problems fitting very large DDE systems in deSolve, possibly because the order of interpolation is lower than the order of integration which can cause problems with the overall accuracy of the solution.

Hairer addressed this problems by implementing the [Dormand Prince](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method) which has "dense output", which means that the solution can be computed at any time point the solver has passed to the same accuracy as the solution itself.  This sidesteps the interpolation problem at the cost of a bit more book-keeping.

I implemented the same integration algorithm and use a [ring buffer](https://github.com/richfitz/ring) to hold the history over time.  This means that the memory required to store the solution does not grow as the total integration length increases (though you still need to pick an amount of memory that scales with the maximum number of steps that span your longest lag at any point in the integration).

This solver is suitable only for nonstiff problems.

The interface is slightly different to the deSolve interface.  A `deSolve` compatible interface may be provided later.
