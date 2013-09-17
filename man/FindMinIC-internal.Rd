\name{FindMinIC-internal}
\alias{splitvars}
\alias{getx}
\alias{gety}
\alias{getIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FindMinIC internal functions
}
\description{
These are internal functions used by FindMinIC.  These should generally \emph{not} be used directly.
}
\usage{
splitvars(fixed)
getx(formula)
gety(formula)
getIC(fit, ictype)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{
a list of strings
}
  \item{formula}{
a formula object
}
  \item{fit}{
a model fit object such as an object returned by lm
}
  \item{ictype}{
type of IC to use, options are "AIC", "AICc", "BIC"
}
}
\details{
splitvars splits the variables in fixed into their component variables.  For instance, splitvars(c("A","B",C*A","-1","A:D","E:F:G")) will return
c("A","B","C","1", "D", "E","F","G"

getx and gety will return the x and y sides of the given formula respectively

getIC returns the AIC, AICc, or BIC for the given model based on which type was passed in
)
}
\value{
see details above

}

\author{
Nicholas Lange, Tom Fletcher, Kristen Zygmunt
}

\seealso{
\code{\link{FindMinIC}} which should be used directly instead of these methods.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
