# dde

An R package that implements the `DOPRI5` solver of [Ernst Hairer](http://www.unige.ch/~hairer/software.html), along with the delay differential equation (DDE) form (originally called `RETARD`).

This is an alternative approach to fitting delay DDE models to using deSolve.  It exists because I have had problems fitting very large DDE systems in deSolve, possibly because the order of interpolation is lower than the order of integration which can cause problems with the overall accuracy of the solution.
