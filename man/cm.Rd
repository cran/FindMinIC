\name{cm}
\alias{cm}
\alias{cmList}
\alias{getNthModel}
\alias{getFirstModel}
\alias{summaryTable}
\alias{summary.cmList}
\alias{print.summary.cmList}
\alias{formula.cm}
\alias{IC.cm}
\alias{summary.cm}
\alias{print.summary.cm}
\title{
Handling \dQuote{\code{cm}} and \dQuote{\code{cmList}} objects
}
\description{
Methods for working with candidate model (\dQuote{\code{cm}}) and candidate model list (\dQuote{\code{cmList}}) objects.
}
\usage{
getNthModel(object, index)
getFirstModel(object)

summaryTable(object, index, \dots)

\method{summary}{cmList}(object, \dots)

\method{print}{summary.cmList}(x, \dots)

\method{formula}{cm}(x, \dots)

\method{IC}{cm}(object)

\method{summary}{cm}(object, \dots)

\method{print}{summary.cm}(x, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a list of candidate models of class \code{"cmList"} that was returned from a call to \code{\link{fmi}} or \code{\link{FindMinIC}}.
}
  \item{index}{
an index into the \code{cmList}, for instance, using index = 4 in getNthModel will return the 4th best candidate model from the list as ranked by IC.
}
  \item{x}{
a candidate model of class \code{"cm"}
}
  \item{\dots}{
some generic methods such as print and summary can take additional arguments, these are passed in to those methods appropriately.
}
}
\details{
Since \code{\link{FindMinIC}} returns a potentially large list of candidate models, functions such as \code{getNthModel}, \code{getFirstModel}, \code{summaryTable} and \code{summary} can be used to either return specific models or summarize the list of models.  Once a particular candidate model (\code{cm}) has been extracted, functions such as \code{IC}, \code{formula}, and \code{summary} can be used to further understand that particular model.
}
\value{
a cmList is a list containing the following components:
\item{results}{ a list of candidate models of class \code{"cm"}. }
\item{data}{ the data passed in including any changes made for \code{groupedData}. }
\item{best}{ the best candidate model }
\item{modeltype}{ the model type provided to the \code{FindMinIC} call }

a cm object is a candidate model containing the following components:
\item{call}{the call to \code{lm} or \code{lme} that generated the candidate model}
\item{IC}{ the \code{IC} of the model fit }
\item{formula}{ the formula of the model fit }
}

\author{
Nicholas Lange, Tom Fletcher, Kristen Zygmunt
}

\seealso{
\code{\link{FindMinIC}}
}
\examples{
     data(iris)

     coly="Sepal.Length"
     fixed="Sepal.Width"
     candidates=c("Species","-1","Sepal.Width:Species")

     results.lm = FindMinIC(coly, candidates, fixed, iris)

     print(summary(results.lm))

     # best model:
     print(summary(getFirstModel(results.lm)))
     print(summaryTable(results.lm, 1))
     print(summary(results.lm$results[[1]]))

     # 4th best model:
     print(summary(getNthModel(results.lm, 4)))
     print(summary(results.lm$results[[4]]))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ print }
