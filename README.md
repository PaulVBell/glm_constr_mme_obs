# Constraining GLM trends from MMEs using observations

See example.rmd for a worked example of using the hierarchical models and
`R/exeter_example.R` for the same models run on UKCP18 outputs.

## Requirements
The R programming language.

A working Stan installation. This can be done using the [cmdstanr package](https://mc-stan.org/docs/cmdstan-guide/installation.html) using the `install_cmdstan()` function.
Note you will also need a working C++ toolchain.
### R packages
- [tidyverse](https://tidyverse.org/)
- [cmdstanr](https://mc-stan.org/docs/cmdstan-guide/installation.html)
- [bayesplot](https://mc-stan.org/bayesplot/)

### Data licence

UKCP18 data is published under the
[Open Government Licence v3.0](http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).
Redistribution is permitted with attribution — CEDA's registered-user
requirement governs access, not reuse. Cite the source dataset as:

> Met Office Hadley Centre (2019): UKCP Local Projections on a 5km grid over the
> UK for 1980-2080. Centre for Environmental Data Analysis, 3 September 2026.
> https://catalogue.ceda.ac.uk/uuid/e304987739e04cdc960598fa5e4439d0