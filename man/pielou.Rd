\name{pielou}
\alias{pielou}
\title{Pielou's Index of non-Randomness}
\description{
  Computes Pielou's index
  \eqn{\alpha} for a spatial point pattern.
}
\usage{
pielou(X, correction=c("none", "Donnelly", "cdf"),
              clipregion=NULL)
}
\arguments{
  \item{X}{
    A spatial point pattern (object of class \code{"ppp"}).
  }
  \item{correction}{
    Character vector.
    The type of edge correction(s) to be applied.
  }
  \item{clipregion}{
    Clipping region for the guard area correction.
    A window (object of class \code{"owin"}).
    See Details.
  }
}
\details{
  Pielou's (1959) index of non-randomness \eqn{\alpha} is a crude
  measure of clustering or ordering of a point pattern.
  Pielou (1959) devised an index of non-randomness based on the average distance between
  randomly selected points to the nearest plant and the density of plants per unit area.

  It is based on the average distance between randomly selected points to the 
  nearest neighbour in the pattern of the same intensity.
  A value \eqn{\alpha<1} suggests ordering, while \eqn{\alpha>1} suggests
  clustering.

  Without correction for edge effects, the value of \eqn{\alpha} will be
  positively biased. Edge effects arise because, for a point of \code{X}
  close to the edge of the window, the true nearest neighbour may
  actually lie outside the window. Hence observed nearest neighbour
  distances tend to be larger than the true nearest neighbour distances.

  The argument \code{correction} specifies an edge correction
  or several edge corrections to be applied. It is a character vector
  containing one or more of the options
  \code{"none"} or \code{"guard"}
  (which are recognised by partial matching).
  These edge corrections are:
  \describe{
    \item{"none":}{
      No edge correction is applied. 
    }
    \item{"guard":}{
      Guard region or buffer area method.
      The observed mean nearest neighbour distance
      for the point pattern \code{X}
      is re-defined by averaging only over those points of \code{X}
      that fall inside the sub-window \code{clipregion}.
    }
  }
  If the argument \code{clipregion} is given, then the selected
  edge corrections will be assumed to include \code{correction="guard"}.

  To perform a test based on the Pielou's index,
  see \code{\link{pielou.test}}.
}
\value{
  A numeric value, or a numeric vector with named components
  \item{naive}{\eqn{R} without edge correction}
  \item{guard}{\eqn{R} using guard region}
  (as selected by \code{correction}). 
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
 - built using the code for 
 \code{\link{clarkevans}}
 written by John Rudge
\email{rudge@esc.cam.ac.uk}
with modifications by
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{clarkevans.test}},
  \code{\link{nndist}},
}
\examples{
  # Example of a clustered pattern
  pielou(redwood)

  # Example of an ordered pattern
  pielou(cells)

  # Random pattern
  X <- rpoispp(100)
  pielou(X)

  # How to specify a clipping region
  clip1 <- owin(c(0.1,0.9),c(0.1,0.9))
  clip2 <- erosion(cells$window, 0.1)
  pielou(cells, clipregion=clip1)
  pielou(cells, clipregion=clip2)
}
\keyword{spatial}
\keyword{nonparametric}
