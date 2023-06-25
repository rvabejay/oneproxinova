\name{SPM}
\alias{SPM}
\title{
Test of Hypothesis on Several Population Means using Analysis of Variance (ANOVA)
}
\description{
Tests of hypothesis on more than two population means using F-Test using Analysis of Variance (ANOVA) technique at a specified significance level.
}
\usage{
SPM(data, factor, obs, alpha = 0.05)
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
}
\value{
SPM returns a dataframe of the results of the Analysis of Variance (ANOVA) given the parameters at the specified level of significance. Signficance of the factor is noted with a "*" on the right side of the p-value.

The decision rule is printed by default. The alternative hypothesis is printed when the factor is significant and the null hypothesis is printed if otherwise. Recommendations on whether to do post-hoc analysis is also printed depending on the significance of the factor.
}
\references{
}
\author{
Rhys Allen Abejay, Ashley Nicole Anyayahan, Jodi Baeza, Hiu Lam Janna Chiu, Mark Louie Lobo, and Norielli Palomares
}
\note{
The function is currently limited to one-way ANOVA (one-factor/independent variable).

Computation of the F-critical value is approximated using Paulson's Approximate Transformation. P-value is also approximated using Smillie and Antsey's Closed Form Appriximation. Function outputs for these values might be slightly different to the R built-in function but similar results are seen on significance tests.

The function will not work for factor levels with unequal number of observations.
}
\seealso{
}
\examples{
x <- c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4))
y <- rexp(16, rate = (1/5))

data <- as.data.frame(cbind(x,y))

SPM(data, 1, 2)
SPM(data, 1, 2, alpha = 0.10)

#this will not work
SPM(data, 2, 1)
}
