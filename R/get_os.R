#' get_os
#'
#' Determine what operating system R is running under.
#'
#' This code is adapted from a post discussing a safe way to
#' detect the operating system that R is running under:
#' https://www.r-bloggers.com/identifying-the-os-from-r/
#'
#' Specifically, this code is designed to return the proper operating
#' system name even if the R function 'Sys.info' is unimplemented.
#'
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
