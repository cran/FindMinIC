\name{FindMinIC}
\alias{FindMinIC}
\alias{FindMinIC.default}
\alias{FindMinIC.formula}
\alias{fmi}
\alias{fmi.default}
\alias{fmi.formula}
\title{
Find Model with Minimum IC
}
\description{
Evaluates all models in a set of candidates, and ranks them by IC such as AIC.  Either lm or lme can be used for the model.
}
\usage{
# Find the minimum IC
\method{FindMinIC}{default}(coly, candidates = c(""), fixed = c(""), data = list(),
          modeltype = "lm", random = ~1, ic = "AIC", ...)
\method{FindMinIC}{formula}(formula, data=list(), na.action=na.omit, fixed = c(""), random = ~1, ...)

# find the minimum IC, fmi is the shorter name form of FindMinIC
\method{fmi}{default}(coly, candidates = c(""), fixed = c(""), data = list(),
    modeltype = "lm", random = ~1, ic = "AIC", ...)
\method{fmi}{formula}(formula, data=list(), na.action=na.omit, fixed = c(""), random = ~1, ...)
%% TODO: use better name than coly?
}
\arguments{
  \item{formula}{
A formula containing the response variable and terms. All the terms of the formula become candidates for inclusion as covariates.
}
  \item{na.action}{
action to use when data contains NAs.  Options include na.omit, na.exclude, na.fail
}
  \item{coly}{
The name of the column to use for the response variable \code{y} of the model
}
  \item{candidates}{
A list of names of columns that are candidates for inclusion as covariates in the model
}
  \item{fixed}{
A list of names of columns (can be empty) that must always be included in every model
}
  \item{data}{
An object containing the variables for use in the model.
}
  \item{modeltype}{
Currently a choice between \code{"lm"} (the default) and \code{"lme"}.  If a model follows the calling convention of \code{\link{lm}}, it might work here, but it is not guaranteed.
}
  \item{random}{
When \code{modeltype} = \code{"lme"}, use \code{random} the same way as would inside a call to \code{\link{lme}} and to indicate the variable for \code{\link{groupedData}}
}
  \item{ic}{
Type of information criterion to used.  Defaults to \code{"AIC"}.  Other options are \code{"AICc"} or \code{"BIC"}
}
  \item{\dots}{
Extra arguments are passed directly into the call to \code{\link{lm}} or \code{\link{lme}}.
}
}
\details{
FindMinIC tries all possible model combinations of the candidate covariates, while always including the same response variable and fixed variables.  It returns a list of candidate models ranked by IC.  The model combinations include all 2-way interactions among the candidate variables.  Other interactions (like age^2) can be directly included in the candidates or fixed lists.
}
\value{
FindMinIC returns a list of candidate models sorted by information criterion IC.  The first model has the "best" IC.  The list is of class(\code{"cmList"}) while each element of that list is of class(\code{"cm"})
see \code{\link{cmList}} for more details
}
\author{
Nicholas Lange, Tom Fletcher, Kristen Zygmunt
}
\references{
\cite{Burnham, K. P.; Anderson, D. R. (2004), "Multimodel inference: understanding AIC and BIC in Model Selection", Sociological Methods and Research 33: 261-304.}
}
\seealso{
\code{\link{getFirstModel}}
}
\examples{

data(iris)

coly="Sepal.Length"
fixed="Sepal.Width"
candidates=c("Species","-1","Sepal.Width:Species")

results.lm = FindMinIC(coly, candidates, fixed, iris)

# model with lowest IC:
first.model = getFirstModel(results.lm)
print(summary(first.model))

# model with 3rd lowest IC:
third.model = getNthModel(results.lm, 3)
print(summary(third.model))

# list of first 5 models, ordered by AIC
print(summary(results.lm)$table[1:5,])

# list of first 5 models, ordered by BIC
results.bic = FindMinIC(coly, candidates, fixed, iris, ic="BIC")
print(summary(results.bic)$table[1:5,])

fm = FindMinIC(Infant.Mortality ~ ., data = swiss)
summary(fm)

fm2 = FindMinIC(Infant.Mortality ~ Fertility + Agriculture + Education * Catholic,
                data = swiss)
summary(fm2)

# list of first 5 models, ordered by AICc
if (require(nlme)) {
  results.aicc = FindMinIC(distance~age, data=Orthodont, 
                           ic="AICc", model="lme",
                           random= ~ 1 | Subject)
  print(summary(results.aicc))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }
