Questions in response to the previous submission:

> Why is such a quick update needed? Please explain. Please also re-read the CRAN policies about submission frequency.

This update removes the compile-time dependency of testthat, so that the package can be installed without having to install testthat and all its down-stream dependencies. This was requested by the ggplot2 team as ggplot2 depends on isoband.

This request was made a few days ago and I'm updating the package in response to the request: https://github.com/wilkelab/isoband/issues/19

Since the update significantly reduces load on CRAN mirrors, it seems worthwhile to do it as soon as possible.

I have also fixed the outdated link that was identified in the previous submission.

## Test environment
* ubuntu 20.04, devel and release
* windows, release
* macOS, release

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.
