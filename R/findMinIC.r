# uses lme and groupedData from nlme
# uses as.sets from sets for power set of X

splitvars <- function(fixed) {
  if (length(fixed) == 0) {
    return(fixed)
  }
  tokens = strsplit(fixed,':')
  fixed.terms = c()
  for (tok in tokens) {
    if (length(tok) > 1) {
      fixed.terms = c(fixed.terms, sub("-","",unlist(tok)))
    } else {
      # currently assuming fixed is in form of either a:b:c or a*b*c or a
      # add more steps if other separators are needed
      splittok = strsplit(c(tok),'\\*')
      fixed.terms = c(fixed.terms, sub("-","",unlist(splittok)))
    }
  }
  return(unique(fixed.terms))
}

getIC <- function(fit, ictype) {
  if (substr(ictype,1,3) == "AIC") {
    icval = AIC(fit)
    if (ictype == "AICc") {
      # add the correction
      npar = length(fit$coefficients) + 1
      n = nobs(fit)
      icval = icval + 2 * npar * (npar + 1) / (n - npar - 1)
    }
  } else if (ictype == "BIC") {
    icval = BIC(fit)
  } else {
    warning(paste("Unrecognized IC option", ictype, "using AIC instead"))
    icval = AIC(fit)
  }
  icval
}            

fmi <- function(coly, ...) UseMethod("fmi")

fmi.default <- function(coly, candidates=c(""), fixed=c(""), data=list(), modeltype="lm", group="", ic="AIC", ...) {
  return(FindMinIC.default(coly, candidates, fixed, data, modeltype, group, ic, ...))
}

fmi.formula <- function(formula, data=list(), ...) {
  return(FindMinIC.formula(formula, data, ...))
}

FindMinIC <- function(coly, ...) UseMethod("FindMinIC")

FindMinIC.default <- function(coly, candidates=c(""), fixed=c(""), data=list(), modeltype="lm", group="", ic = "AIC", ...){
  fixed.vars = splitvars(fixed)
  cand.vars = splitvars(candidates)
  if (modeltype == "lme") {
    comb.vars = unique(c(fixed.vars, cand.vars, coly, group))
  } else {
    comb.vars = unique(c(fixed.vars, cand.vars, coly))
  }
  comb.vars = comb.vars[which(comb.vars!="1")]
  tmp.df = subset(data, select=comb.vars)
  
  if (modeltype == "lme") {
    # just using intercept for grouping for now
    # because lme might be less stable/robust for other options
    # TODO contemplate allowing the user to specify the variable(s) for grouping
    gdcall = parse(text = paste("groupedData(",
                     coly, "~ 1 |", group,
                     ", data=tmp.df)"))
    tmp.gds = eval(gdcall)
    tmp.gds = na.omit(tmp.gds)
  } else {
    tmp.gds = na.omit(tmp.df)
  }
  fixed.model = ""
  firstelem = TRUE
  # for the rest of the fixed elements
  for (f in fixed) {
    if (f != "") {
      if (firstelem) {
        fixed.model = paste(fixed.model, f, sep="")
        firstelem = FALSE
      } else {
        fixed.model = paste(fixed.model, "+", f, sep="")
      }
    }
  }

  results = list()

#
# define model space for group analysis
#
  xpnames.set = as.set(candidates)
  Xp = parse(text = 2^xpnames.set) # power set
  M = length(Xp)

  for(m in 1:M){ # all possible models given covariates in candidates
    if(m == 1) model = "1"
    if(m != 1){
      model = NULL
      op = "+"
      xtext = as.character(Xp[[m] ] )[-1] # remove string "list()"
      P = length(xtext)
#
# the following 5 lines could be more well thought out, I think...
#
      if(P == 1) op = ""
      for(p in 1:(P - 1) ){ # weird for P == 1, but it works...
        model = paste(model, xtext[p], op, sep = "")
      }
      if(P > 1) model = paste(model, xtext[P], sep = "")
    }
    if (modeltype == "lme") {
      call = parse(text = paste("withRestarts(lme(",
                     coly,
                     " ~ ",
                     fixed.model,
                     " + ",
                     model,
                     ", random=~1",
                     ", na.action = 'na.omit'",
                     ", method = 'ML'",
                     ", data = tmp.gds",
                     ",...)",
                     ")",
                     sep = "") )
      fit = try(eval(call), silent = TRUE); if(class(fit) == "try-error") next
    
      res = list(call = fit$call, IC = getIC(fit, ic), ictype = ic, formula = formula(fit))
      class(res) = "cm" # for candidate model
      results[[length(results)+1]] = res

    } else {
      # TODO ok to convert modeltype directly into R command?
      call = parse(text = paste(modeltype,
                     "(",
                     coly,
                     " ~ ",
                     fixed.model,
                     " + ",
                     model,
                     ", data = tmp.gds",
                     ",...)",
                     sep = "") )
      fit = try(eval(call), silent = TRUE); if(class(fit) == "try-error") next
      res = list(call = fit$call, IC = getIC(fit, ic), ictype = ic, formula = formula(fit))
      class(res) = "cm" # for candidate model
      results[[length(results)+1]] = res
    }
  }

  if (length(results) > 1) {
    results = results[order(as.numeric(lapply(results, IC.cm)))]
  }
  first = results[[1]]
  if (modeltype == "lme") {
    # for lme, need to do REML instead of ML here
    call = parse(text = paste("lme(",
                   deparse(first$call$fixed, width.cutoff=500),
                   ", random=~1",
                   ", data = tmp.gds",
                   ", na.action = 'na.omit'",
                   ", method = 'REML'",
                   ",...)",
                   sep = ""
                   ) )
  } else {
    call = first$call
  }
  fit = eval(call)
  first$fit = fit
  answer = list(results=results, data=tmp.gds, first=first, modeltype=modeltype)
  class(answer) = "cmList"

  return(answer)
} # end FindMinIC.default

FindMinIC.formula <- function(formula, data=list(), ...)
{
    mf <- model.frame(formula=formula, data=data)
    nms <- names(mf)
    coly <- nms[[1]]
    candidates <- nms[-1]
    est <- FindMinIC.default(coly, candidates, data=data, ...)
    est$call <- match.call()
    est$formula <- formula
    est$terms <- attr(mf, "terms")
    est$data <- mf

    est
}
