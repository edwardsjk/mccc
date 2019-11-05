# mccc

An r package to account for misclassification between censoring and competing events.

mccc can be installed by running

```r
library(devtools)
install_github("edwardsjk/mccc")
```
## Functions

The package contains the following functions.
1. `comp.cens.cor`: A function to estimate a risk function accounting for misclassification between censoring and competing events due to loss to follow-up
2. `bootstrap.cc`: A function to produce bootstrapped standard errors

## Inputs

All functions require inputting the main study data and specifying the time and event type variables. Functions additionally require the user to input a vector of probabilities denoting the probability that participant $i$ should have been censored, given that he or she appeared to be lost. All functions require the user to specify the maximum follow-up time tau. Details on function calls can be found using `?comp.cens.cor` and `?bootstrap.cc`.

## Outputs

`comp.cens.cor` returns a dataframe with a vector of unique event times and the corresponding cumulative incidence estimate. If only the cumulative incidence estimate at the final timepoint is desired, use

```r
ci <- tail(prop<-comp.cens.cor(data=art, tau=2*365.25, t="t", delta="j", p=0.24), n = 1)
```
