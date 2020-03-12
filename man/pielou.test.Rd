\name{pielou.test}
\alias{pielou.test}
\title{Pielou's Test}
\description{
  Performs Pielou's test of nonrandomness
  for a spatial point pattern.
}
\usage{
pielou.test(X, ...,
               correction="none",
               clipregion=NULL,
               alternative=c("two.sided", "less", "greater",
                             "clustered", "regular"),
               nsim=999)
}
\arguments{
  \item{X}{
    A spatial point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{Ignored.}
  \item{correction}{
    Character string.
    The type of edge correction to be applied.
    See \code{\link{clarkevans}}
  }
  \item{clipregion}{
    Clipping region for the guard area correction.
    A window (object of class \code{"owin"}).
    See \code{\link{pielou}}
  }
  \item{alternative}{
    String indicating the type of alternative for the
    hypothesis test. Partially matched.
  }
  \item{nsim}{
    Number of Monte Carlo simulations to perform, if a Monte Carlo
    p-value is required.
  }
}
\details{
  This command uses Pielou's  (1959) aggregation index \eqn{\alpha} 
  as the basis for a crude test of clustering or ordering of a point pattern.
  Pielou's index is computed by the function
  \code{\link{pielou}}. See the help for \code{\link{pielou}}
  for information about Pielou's index \eqn{\alpha} and about
  the arguments \code{correction} and \code{clipregion}.

  This command performs a hypothesis test of clustering or ordering of
  the point pattern \code{X}. The null hypothesis is Complete
  Spatial Randomness, i.e.\ a uniform Poisson process. The alternative
  hypothesis is specified by the argument \code{alternative}:
  \itemize{
    \item \code{alternative="less"} or \code{alternative="clustered"}:
    the alternative hypothesis
    is that \eqn{\alpha > 1} corresponding to a clustered point pattern;
    \item \code{alternative="greater"} or \code{alternative="regular"}:
    the alternative hypothesis
    is that \eqn{\alpha < 1} corresponding to a regular or ordered point pattern;
    \item \code{alternative="two.sided"}:
    the alternative hypothesis is that \eqn{\alpha != 1}
    corresponding to a clustered or regular pattern.
  }
  
  Pielou's index \eqn{\alpha} is computed for the data
  as described in \code{\link{pielou}}.

  If \code{correction="none"} and \code{nsim} is missing,
  the \eqn{p}-value for the test is computed using a X-squared 
  distribution  as proposed by Pielou (1959) and referring the
  statistic to the Chi-squared distribution.

  Otherwise, the \eqn{p}-value for the test is computed
  by Monte Carlo simulation of \code{nsim} realisations of
  Complete Spatial Randomness conditional on the
  observed number of points.
}
\value{
  An object of class \code{"htest"} representing the result of the test.
}
\references{
  Pielou, E. C. (1959) The use of point-to-plant distances
  in the study of the pattern of plant populations. \emph{J. Ecol} 
  \bold{47}, 607--613.
  
  Pielou, E. C. (1962) The use of plant-to-neighbor distances
  for the detection of competition. \emph{J. Ecol} \bold{50}, 
  357--367.
}
\author{
  Andrew Sanchez Meador
  \email{andrew.sanchezmeador@nau.edu}
  \url{http://quantitativeecology.org}
}
\seealso{
  \code{\link{pielou}}
}
\examples{
  # Redwood data - clustered
  pielou.test(redwood)
  pielou.test(redwood, alternative="greater")
}
\keyword{spatial}
\keyword{nonparametric}
\keyword{htest}
