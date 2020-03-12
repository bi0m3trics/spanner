pielou <- function(X, correction=c("none"), clipregion=NULL)
{
  verifyclass(X, "ppp")
  W <- X$window

  # validate correction argument
  gavecorrection <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(none="none",
                             guard="guard",
                           multi=TRUE))

  # guard correction applied iff 'clipregion' is present
  isguard <- "guard" %in% correction
  askguard <- any(isguard)
  gaveguard <- !is.null(clipregion)
  if(gaveguard)
    clipregion <- as.owin(clipregion)
  if(askguard && !gaveguard) {
    warning("guard correction not performed; clipregion not specified")
    correction <- correction[!isguard]
  } else if(gaveguard && !askguard) 
    correction <- c(correction, "guard")

  result <- pielouCalc(X, correction, clipregion)
  if(length(result) == 1) result <- unname(result)
  return(result)
}

pielou.test <- function (X, ..., correction = "none", clipregion = NULL,
                         alternative = c("two.sided", "less", "greater"), nsim = 999) 
{
  Xname <- short.deparse(substitute(X))
  verifyclass(X, "ppp")
  W <- X$window
  correction <- pickoption("correction", correction, c(none = "none", 
                                                       guard = "guard"))
  switch(correction, none = {
    corrblurb <- "No edge correction"
  }, guard = {
    if (is.null(clipregion)) stop("clipregion not specified")
    clipregion <- as.owin(clipregion)
    corrblurb <- "Guard correction"
  })
  if (missing(alternative) || is.null(alternative)) 
    alternative <- "two.sided"
  alternative <- pickoption("alternative", alternative, c(two.sided = "two.sided", 
                                                          less = "less", clustered = "less", greater = "greater", 
                                                          regular = "greater"))
  altblurb <- switch(alternative, two.sided = "two-sided", 
                     less = "mean squared distance to nearest neighbor less than expected under CSR (clustered)", 
                     greater = "mean squared distance to nearest neighbor greater than expected under CSR (regular)")
  statistic <- pielouCalc(X, correction = correction, clipregion = clipregion, 
                          working = TRUE)
  working <- attr(statistic, "working")
  if (correction == "none") {
    intensity <- with(working, intensity)
    p.value <- switch(alternative, greater = pchisq(statistic * 2 * intensity, 2 * intensity), 
                      less = 1 - pchisq(statistic * 2 * intensity, 2 * intensity), 
                      two.sided = 2 * (1 - pchisq(statistic * 2 * intensity, 2 * intensity)))
    pvblurb <- "X-squared"
  }
  else {
    sims <- numeric(nsim)
    intensity <- working$intensity
    for (i in 1:nsim) {
      Xsim <- rpoispp(npts, win = W)
      sims[i] <- pielouCalc(Xsim, correction = correction, 
                            clipregion = clipregion)
    }
    prob <- mean(sims <= statistic)
    p.value <- switch(alternative, less = prob, greater = 1 - 
                        prob, two.sided = 2 * min(prob, 1 - prob))
    pvblurb <- paste("Monte Carlo test based on", nsim, "simulations of CSR")
  }
  statistic <- as.numeric(statistic)
  names(statistic) <- "Alpha"
  out <- list(statistic = statistic, p.value = round(p.value, 6), alternative = altblurb, 
              method = c("Pielou's test", corrblurb, pvblurb), data.name = Xname)
  class(out) <- "htest"
  return(out)
}

pielouCalc <- function (X, correction = "none", clipregion = NULL, working = FALSE) 
{
  W <- X$window
  area <- area.owin(W)
  npts <- npoints(X)
  intensity <- npts/area
  if (npts == 0) 
    return(NA)
  nncrossX <- nncross(X=runifpoint(npts,win=W), Y=X)
  w_bar <- mean(nncrossX$dist^2)
  statistic <- NULL
  if (working) 
    work <- list(area = area, npts = npts, intensity = intensity, 
                 w_bar = w_bar)
  if ("none" %in% correction) {
    Rnaive <- pi*w_bar*intensity
    statistic <- c(statistic, naive = Rnaive)
  }
  if ("guard" %in% correction && !is.null(clipregion)) {
    ok <- inside.owin(X, , clipregion)
    w_bar_guard <- mean(nncrossX$dist[ok]^2)
    Rguard <- pi*w_bar_guard*intensity
    if (working) 
      work <- append(work, list(w_bar_guard = w_bar_guard))
    statistic <- c(statistic, guard = Rguard)
  }
  if (working) 
    attr(statistic, "working") <- work
  return(statistic)
}

