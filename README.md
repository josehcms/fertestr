
# fertestr

# Demographic tools for the assessment of fertility and parity data.
Date: 2019-12-19

By the end of this initial development phase, this repository will contain indirect fertility estimation functions in a package format. The functions to be implemented in here can be found at the [Tools for Demographic Estimation](https://demographicestimation.iussp.org/) This project is commissioned and financed by the [UN Population Division](http://www.un.org/en/development/desa/population/). [Everton E C Lima](https://twitter.com/notreve81) and I, [José H C Monteiro da Silva](https://josehcms.github.io/), are the authors in direct collaboration with [Patrick Gerland](https://www.researchgate.net/profile/Patrick_Gerland) and [Helena C Castanheira](https://cl.linkedin.com/in/helena-cruz-castanheira-954061105). 

## How to load ```fertestr``` package in R?
```r
install.packages("devtools")

devtools::install_github("josehcms/fertestr")
```

## Note
We are still at **initial development**, so major bugs might be frequent.

## Getting started
Soon we will add some tutorials on the usage of the package. By now, we recommend you to explore the help files which contain some worked examples and description of functions.

```
library(fertestr)
# Brass PF ratio fertility correction:
?fertBrassPF

# Average parities computing with several adjustment and assessment options
?prtyAverage

# El-Badry correction
?prtyElBadry

```
## References
[Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013. Tools for Demographic Estimation. Paris: International Union for the Scientific Study of Population.](https://demographicestimation.iussp.org/)
