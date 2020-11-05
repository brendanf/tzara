
<!-- README.md is generated from README.Rmd. Please edit that file -->

tzara
=====

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/brendanf/tzara.svg?branch=master)](https://travis-ci.com/brendanf/tzara)
[![Codecov test
coverage](https://codecov.io/gh/brendanf/tzara/branch/master/graph/badge.svg)](https://codecov.io/gh/brendanf/tzara?branch=master)
<!-- badges: end -->

To reduce computational complexity,
[dada2](https://benjjneb.github.io/dada2/index.html) only uses
non-singletons as seeds for denoising. For this strategy to work, each
true sequence must be represented by at least two identical reads.
Especially with long amplicons, the probability of two reads having
exactly the same errors is much lower than the probability of being
error-free, so in practice this means that each true sequence must have
two error-free reads. This becomes problematic for rare sequences in
long amplicon libraries.

Tzara (named after [Tristan
Tzara](https://en.wikipedia.org/wiki/Tristan_Tzara), a central figure in
the Dada art movement, who described the [“cut-up”
technique](https://en.wikipedia.org/wiki/Cut-up_technique) for poetry)
implements an alternative method for dealing with long reads: - Cut
reads of the targeted region into smaller domains using hidden Markov
models (via [Hmmer](https://hmmer.org)) or covariance models (via
[Infernal](http://eddylab.org/infernal/)) - Use dada2 to separately
denoise each sub-region/domain - Concatenate denoised sequences for the
different domains which originated in the same read to get denoised
sequences for the full read

There are also methods for deeper chimera checking and attempting to
preserve information from reads which dada2 was unable to successfully
map to a read.

Installation
------------

`tzara` is yet available from CRAN or Bioconductor, but you can install
the latest github version using:

``` r
remotes::install_github("brendanf/tzara")
```

Example
-------

Coming soon!
