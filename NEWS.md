## dde 1.1.0

* New batch mode for use with `difeq`, `difeq_batch` to run models with varying parameters (VIMC-1445)

## dde 1.0.1

* Several memory errors fixed [#22](https://github.com/mrc-ide/dde/pull/22)

## dde 1.0.0

* Initial public release
* Validate the derivatives are finite (non-NA, non-infinite) before the first step, avoiding cryptic errors [#13](https://github.com/mrc-ide/dde/issues/13)
* Verbose output and a callback interface to get information about each step (or each evaluation) [#12](https://github.com/mrc-ide/dde/issues/12)

## dde 0.0.6

* Replication interface for `difeq` as `difeq_replicate` [#10](https://github.com/mrc-ide/dde/issues/10)
