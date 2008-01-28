.onAttach <- function(...)
{


   # figure out year automatically (probably could be done more elegantly)
   date <- date()
   x <- regexpr("[0-9]{4}", date)
   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

   # echo output to screen
   cat("##\n## On Bayesian analysis of Threshold autoregressive model (BAYSTAR)\n")
   cat("## Copyright (C) 2007-", this.year,
      " Cathy W. S. Chen, Edward M.H. Lin, F.C. Liu, and Richard Gerlach\n", sep="")
   require(mvtnorm, quietly=TRUE)
}