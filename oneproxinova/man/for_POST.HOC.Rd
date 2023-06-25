\name{POST.HOC}
\alias{POST.HOC}
\title{
Post-Hoc Analysis using Tukey's Honestly Significant Difference (HSD)
}
\description{
Post-hoc analysis implementing Tukey's Honestly Significant Difference (HSD) test at a specified significance level.
}
\usage{
POST.HOC(data, factor, obs, alpha = 0.05)
}
\arguments{
  \item{data}{
dataset to be used in the data analysis
}
  \item{factor}{
column index or name of the factor/independent variable
}
  \item{obs}{
column index or name of the response variable
}
  \item{alpha}{
level of significance of the multiple comparison of means, default: alpha = 0.05
}
}
\details{
If alpha is not specified, it assumes the default value of 0.05, or a 0.95 confidence level.
Critical value is computed from the Standardized Range Distribution, q, and the t-values
are computed using the point estimates of the mean differences between the two populations and the stadard error.
}
\value{
POST.HOC returns a dataframe of the results of the pairwise mean comparisons at the specified level of significance. The number of the pairwise mean comparisons are computed using the number of factor levels. Significant t-values are marked with "*" on the right of the value.

Null hypothesis and decision rule are printed by default, for reference purposes.
}
\references{
}
\author{
Rhys Allen Abejay, Ashley Nicole Anyayahan, Alexis Mae Roxas, and Daniella Rose Dizon
}
\note{
Significance of the t-value is approximated. Using the R function TukeyHSD to determine the adjusted p-values at each comparison is recommended for better data analysis and interpretation.

The function will not work for factor levels with unequal number of observations and only works for a one-factor(independent) variable.
}

\seealso{
}
\examples{
x <- c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4))
y <- rexp(16, rate = (1/5))

data <- as.data.frame(cbind(x,y))

POST.HOC(data, 1, 2)
POST.HOC(data, 1, 2, alpha = 0.10)

#this will not work
POST.HOC(data, 2, 1)
}
