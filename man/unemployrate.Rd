\name{unemployrate}
\alias{unemployrate}
\docType{data}
\title{ U.S. monthly civilian unemployment rate }
\description{
  U.S. monthly civilian unemployment rate, seasonally adjusted and measured in percentage, from January 1948 to March 2004.  The data set is available in Tsay (2005).
}
\usage{data(unemployrate)}
\format{
  unemployrate is a data frame with 675 observations
}
\source{
  The data are obtained from the Bureau of Labor statistics, Department of Labor.
}
\references{
 Tsay, R.S. (2005) Analysis of Financial Time Series. Second ed. Wiley, Hoboken.
}
\examples{
data(unemployrate)
plot.ts(unemployrate)
}
\keyword{datasets}
