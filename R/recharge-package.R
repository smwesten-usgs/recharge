#'Functions for estimating groundwater recharge from streamflow data
#'
#'This package contains several algorithms designed to estimate groundwater
#' recharge (net infiltration) from daily streamflow data.
#'
#'\tabular{ll}{ Package: \tab recharge\cr
#'Type: \tab Package\cr
#'License: \tab CC0\cr
#'Depends: \tab lubridate, dplyr\cr }
#'This package contains functions that manage or manipipulate
#'hydrologic daily-value data. The functions in this package
#'are listed below, grouped by a general description of the task.
#'
#'Functions to compute baseflow from the streamflow record and
#'functions that summarize those computations.\cr
#'\code{\link{bfi}}\cr
#'\code{\link{hysep}}\cr
#'\code{\link{part}}\cr
#'\code{\link{aggregate.baseflow}}\cr
#'\code{\link{plot.baseflow}}\cr
#'\code{\link{print.baseflow}}\cr
#'
#'Functions that analyze groundwater record, primarily for the
#'estimation of recharge.\cr
#'\code{\link{fall}}\cr
#'\code{\link{rise}}\cr
#'\code{\link{wtf}}\cr
#'\code{\link{aggregate.rise}}\cr
#'\code{\link{confirm.fall}}\cr
#'\code{\link{plot.fall}}\cr
#'\code{\link{plot.rise}}\cr
#'\code{\link{print.fall}}\cr
#'\code{\link{print.rise}}\cr
#'\code{subset.fall}\cr
#'
#'
#'Functions for streamflow recession analysis.\cr
#'\code{\link{join}}\cr
#'\code{\link{recess}}\cr
#'\code{\link{rora}}\cr
#'\code{\link{aggregate.rora}}\cr
#'\code{\link{confirm.recess}}\cr
#'\code{\link{plot.recess}}\cr
#'\code{\link{plot.rora}}\cr
#'\code{\link{print.recess}}\cr
#'
#' @name recharge-package
#' @aliases recharge-package recharge
#' @docType package
#' @author Steve Westenbroek
#'
#' @seealso \code{\link[smwrData:smwrData-package]{smwrData}}
#' @references Need those for the various routines!
#' @keywords package
#' @import lubridate
#' @import dplyr
NULL
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This information is preliminary or provisional and
is subject to revision. It is being provided to meet
the need for timely best science. The information
has not received final approval by the U.S. Geological
Survey (USGS) and is provided on the condition that
neither the USGS nor the U.S. Government shall be held
liable for any damages resulting from the authorized
or unauthorized use of the information.")
}
