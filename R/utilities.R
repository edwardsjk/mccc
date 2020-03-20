# Utilites for accounting for misclassification between censoring and competing events

require(resample)



#' A function to estimate a nonparametric cumulative incidence function accounting for misclassification between competing events and censoring events
#' #'
#' @param data input data source containing main study data, specifically time variable and outcome type
#' @param tau final follow-up time
#' @param t character string denoting name of time variable
#' @param delta character string denoting name of event type variable
#' @param weights character string denoting name of variable containing weights (e.g., survey sampling weights, inverse probability weights, etc)
#' @param p a character string denoting name of the variable containing the probability that person i is in care elsewhere
#' @return data frame with event times and cumulative incidence function for outcome of interest
#' @export
#' @examples
#'  prop<-comp.cens.cor(data=art, tau=2*365.25, t="t", delta="j", p=0.24)


comp.cens.cor<-function(data=data, tau, t="t", delta="delta", weights="NULL",  p="p"){
  data$t <- data[, t]
  data$delta <- data[, delta]
  data$p_ <- data[, p]
  if(weights!="NULL")
  {data$weights <- data[, weights]
  }
  else{data$weights <- 1}

  c.data<-data[!is.na(data$delta),  ]
  c.data$y<-c.data$delta
  c.data$w <- 1
  c.data$p_ <- 1

  #if delta is missing make 2 copies weighted by p
  cen.data <- data[is.na(data$delta), ]
  cen.data$y <- 0
  cen.data$w <- cen.data$p_
  comp.data <- data[is.na(data$delta), ]
  comp.data$y <- 2
  comp.data$w <- 1 - comp.data$p_

  part <- as.data.frame(rbind(c.data, cen.data, comp.data))
  part$new_weight <- part$w * part$weights
  riskfn <- risk(data=part, t="t", delta="y", tau=tau, weights="new_weight")

  return(riskfn)
}

# a helper function which returns a matrix, and takes column vectors as arguments for mean and sd
normv <- function( n , mean , sd , seed, i){
  set.seed(seed+i)
  out <- rnorm( n*length(mean) , mean = mean , sd = sd )
  return( matrix( out , ncol = n , byrow = FALSE ) )
}

#' A helper function to bootstrap risk function correcting for misclassification between censoring and competing events
#' @param B number of bootstrap samples to draw
#' @param bdat1 input data source containing main study data
#' @param tau final follow-up time for which risks and standard errors will be reported
#' @param t character string denoting name of time variable
#' @param delta character string denoting name of event type variable
#' @param weights character string denoting name of variable containing weights (e.g., survey sampling weights, inverse probability weights, etc)
#' @param p a character string denoting name of the variable containing the probability that person i is in care elsewhere
#' @param sep a character string denoting the name of the variable containing the standard error of the estimated probability that person i is in care elsewhere
#' @param seed seed
#' @return standard error of the cumulative incidence (i.e., risk) function at time tau
#' @export
#' @examples
#'  se <- bootstrap.cc(B, data, t = "t", delta = "j", tau, p = "p", sep = "selogitp", seed=123)
#'


bootstrap.cc <- function(B, bdat1, t="t", delta="delta", tau , p="p", sep="sep", seed=456){
  set.seed(seed)
  bdat1$p_ <- bdat1[, p]
  bdat1$sep_ <- bdat1[, sep]
  bdat1$t <- bdat1[, t]
  bdat1$delta <- bdat1[, delta]
  datbi<-samp.bootstrap(nrow(bdat1), B)
  pe <- matrix(nrow = B, ncol = 1)
  gamma <- matrix(nrow = nrow(bdat1), ncol = B)

  for(i in 1:B){ # in each bootstrap sample, ...
    datbii<-bdat1[datbi[,i],] # resample main study data
    gamma <-normv(1, logit(datbii$p_), datbii$sep_, seed, i) #draw gamma (logit(p)) from its asympotitc distribution
    datbii$pbi<-expit(gamma)
    # apply function to correct estimates in each bootstrap sample
    pe[i] <- tail(comp.cens.cor(data=datbii, t="t", delta="delta", tau=tau, p = "pbi")$ci1, n=1)
  }

  se<-apply(pe, 2, sd)

  return(se)
}



#' A helper function to estimate a nonparametric cumulative incidence function to the level needed for plotting
#' @param data input data source containing main study data, specifically time variable and outcome type
#' @param tau final follow-up time
#' @param t character string denoting name of time variable
#' @param delta character string denoting name of event type variable
#' @param weights character string denoting name of variable containing weights (e.g., survey sampling weights, inverse probability weights, etc)
#' @param subvar a character string denoting the name of any subsetting variable (not called by other functions in this package)
#' @param subval the level of subvar of interest
#' @return data frame with event times and cumulative incidence function for outcome of interest
#' @export
#' @examples
#'  r1<-risk(data=part, t="t", delta="y", tau=tau, weights="new_weight")

risk <- function(data = data, t = "t", delta = "y", level=1,  tau=4, weights = "NULL", subvar = "NULL", subval = NA){
  library(survival)
  data$t <- data[, t]
  data$delta <- data[, delta]
  data$status <- as.numeric(data$delta!="0")
  if(weights!="NULL")
  {data$weights <- data[, weights]
  }
  else{data$weights <- 1}
  data <- data[with(data, order(t)),]
  if(subvar!="NULL")
  {
    data$subvar <- data[, subvar]
    ajs <- survfit(Surv(t, status*delta, type='mstate')  ~ 1, data = data, weights=weights, subset = (subvar == subval))
  }
  else{ajs <- survfit(Surv(t, status*delta, type='mstate')  ~ 1, data = data, weights=weights)}

  ci <- ajs$pstate[,level+1]
  times <- data$t[data$t<=tau]
  t1 <-ajs$time
  cif1 <- stepfun(ajs$time, c(0, ci), right=TRUE)
  ci1 <- cif1(times)
  subset <- subval
  risk1<-as.data.frame(cbind(times, ci1, subset))
  risk1$subset<-as.character(risk1$subset)
  return(risk1)
}

### Example usage

# corrected.risk <- tail(comp.cens.cor(data=art, t="t", delta="j", tau=365.25*2, p="p"), n= 1)
# se <- bootstrap.cc(200, art, t="t", delta="j", tau = 365.25*2,  seed=456, p = "p", sep = "sep")
