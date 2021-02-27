#' A simple function to download some example MLS data in a ponderosa pine forest north of Flagstaff AZ
#'
 #' @export
getExampleData <- function()
{
  options(timeout=180)
  spannerPath = find.package("spanner", lib.loc = NULL, quiet = FALSE, verbose = getOption("verbose"))
  download.file("http://quantitativeecology.org/files/spanner/Pine_Example.laz", 
                destfile=paste0(spannerPath, "/extdata/Pine_Example.laz"), method="internal", mode="wb")
}