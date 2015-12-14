#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2013
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# Author: JKP, IBM SPSS
# Version = 1.0.0

# history
# 16-Oct-2013 - original version

helptext="Compute nonparametric regression estimates for
a continuous or categorical dependent variable.

STATS NONPAR REGR DEPENDENT=varname INDEP=varnames
ID=idvariable
/BANDWIDTH REGTYPE = LOCALCONSTANT* or LOCALLINEAR
ESTIMATORCONT = CVLS* or CVAIC ESTIMATORCAT = CVML* or CVLS or NORMALREF
SCALING = SCALE* or RAW 
CVTYPE = FIXED* or GENERALIZEDNN or ADAPTIVENN
CVKERNELTYPE = GAUSSIAN* or EPANECHNIKOV or UNIFORM
CVORDER = 2* or 4 or 6 or 8
UKERNELTYPE = AITCHISONAITKEN* or LIRACINE
OKERNELTYPE = WANGVANRYZIN* or LIRACINE
NSTART = integer
/OPTIONS MISSING = OMIT* or STOP MAXITER = integer
VALTOLERANCE = number LOCTOLERANCE = number
/OUTPUT SIGTESTS = NO* or YES NVARS = integer
BOOTMETHOD = IID* or WILD or WILDRADEMACHER
BOOTREPS = integer RNSEED = number
PLOTS = RESPONSE* or GRADIENT or NONE INCLUDEPTS = YES* or NO
ERRORBNDS = NONE* or ASYMPTOTIC or BOOTSTRAP
/SAVE DATASET = datasetname RESIDS = NO* or YES PREDICTED = NO* or YES
/HELP

* indicates a default choice.

Note: Statistical details for the technical specifications can be found
in the R np module documentation and associated references.

Note: The computations in this procedure can be very
time consuming.

DEPENDENT and INDEP are the only required specifications.

Example:
STATS NONPAR REGR DEPENDENT=salary INDEP=jobcat jobtime educ.

DEPENDENT and INDEP specify the dependent and independent variables.
The model type is determined by whether the dependent variable 
measurement level is continuous (scale) or categorical (nominal
or ordinal).
Categorical variables are converted to factors.
Some other specifications depend on the measurement level as well.

ID specifies an id variable to be included in the output dataset
if one is created.  This can be useful if merging the dataset
back with the input.

REGTYPE specifies the kernel regression estimator: locally constant or 
locally linear.

ESTIMATORCONT and ESTIMATORCAT specify how to select bandwidths for
continuous or categorical dependent variables, respectively.  Only
one applies.  CVLS specifies least squares cross validation, and
CVAIC specifies expected Kullback-Leibler cross validation.  CVML
is likelihood cross validation, and NORMALREF computes a rule of
thumb validation.

SCALING specifies whether the bandwidths are scale factors or
raw bandwidths.

CVTYPE specifies the bandwidth type for continuous variables:
fixed, generalized nearest neighbors, or adaptive nearest neighbors.

CVKERNELTYPE specifies the continuous kernel type.

CVORDER determines the kernel order.  It does not apply to
the uniform kernel.

OKERNELTYPE specifes the ordinal kernel type.

UKERNELTYPE specifies the nominal kernel type.

The objective function has multiple extrema, and there is the
possibility of finding a local extremum rather than the global
one.  NSTART specifies how many times to restart the estimation
from a different random starting point in order to guard against
this possibility.  The default is the number of independent
variables up to five.

MISSING specifies whether cases with missing values should be
omitted or cause the estimation to stop.

MAXITER specifies the maximum number of iterations and
defaults to 10,000.

VALTOLERANCE specifies the convergence tolerance on the
value of the cross-validation function and defaults to 1.19e-07.

LOCTOLERANCE specifeis the convergence tolerance on the location
of the cross-validation function and defaults to 1.49e-08.

SIGTESTS specifies whether to carry out significance tests
for the independent variables.  These tests are not available
when the dependent variable is categorical.

NVARS indicates which independent variables should
be tested.  The first NVARS independent variables up to
all of them are tested.

BOOTMETHOD specifies the bootstrapping method for the significance
tests.  IID specifies independent, identically distributed draws.
WILD uses a wild bootstrap, and WILDRADEMACHER use a wild
bootstrap with Rademacher variables.  WILD is useful in the
presence of heteroscedastic errors. The Rademacher distribution
has probably 1/2 on -1 and 1, and these variables are combined
with the wild bootstrap.

BOOTREPS specifies the number of bootstrap replications and
defaults to 399.

RNSEED specifies an integer seed for the random number generator
and defaults to 42.  This allows replication of results from run
to run.

PLOTS determines what plots are produced.  The reponse function
or the gradients can be plotted.  Plotting gradients when the
dependent variable is categorical may not be useful with more
than a few independent variables and may fail with more than 5.

INCLUDEPTS applies only to response plots and specifies whether
the data points are included in the plot.

ERRORBNDS specifies whether to include +/- two standard
deviation confidence intervals in the plots, where the
intervals may be asymptotic or come from bootstrapping.

SAVE specifies the creation of a dataset containing casewise
results.

DATASET is the name for a new dataset to hold the results.
The dataset name must not already be in use.  If an ID
variable was specified, it is included in the dataset.
Otherwise the cases are numbered sequentially.

RESIDUALS and PREDICTED specify whether to include variables
holding these values in the dataset.  At least one must 
be specified if DATASET is used.

STATS NONPAR REGR /HELP prints this help and does nothing else.

This procedures uses the R np package by
Tristen Hayfield and Jeffrey S. Racine.
"

npext <- function(dep, indepnp=NULL, id=NULL, regtype="localconstant",
	estimator="cvls", estimatorcat="cvml", scaling="scale", cvtype="fixed", 
  cvkerneltype="gaussian",
	cvorder=2, ukerneltype="aitchisonaitken", okerneltype="wangvanryzin",
	nstart=NULL,
	missingv="omit", maxiter=10000, valtolerance=1.19e-07,
	loctolerance=1.49e-08, sigtests=FALSE, nvars=999,
	bootmethod="iid", bootreps=399, rnseed=42,
	plots="response", includepts=TRUE, errorbnds="none",
	dataset=NULL, resids=FALSE, predicted=FALSE) {
	
	nindepnp = length(indepnp)
	nstart = ifelse(is.null(nstart), min(5, nindepnp), nstart)
	frml = buildfrml(dep, indepnp)
	allargs = as.list(environment())  # for passing external spec to other functions

  setuplocalization("STATS_NONPAR_REGR")

  # A warnings proc name is associated with the regular output
  # (and the same omsid), because warnings/errors may appear in
  # a separate procedure block following the regular output
  procname=gtxt("Nonparametric Regression")
  warningsprocname = gtxt("Nonparametric Regression: Warnings")
  omsid="STATSNPREG"
  warns = Warn(procname=warningsprocname,omsid=omsid)

  tryCatch(library(np, quietly=TRUE), error=function(e){
      warns$warn(gtxtf("The R %s package is required but could not be loaded.","np"),
          dostop=TRUE)
      }
  )
	if (is.null(dataset) && (resids || predicted)) {
		warns$warn(gtxt("A dataset name must be specified if residuals or predicted values are requested"),
		dostop=TRUE)
	}
	dscheck(dataset, warns)
	naaction = list(omit=na.omit, stop=na.fail)[[missingv]] #can't use ifelse here
	regtype = ifelse(regtype == "localconstant", "lc", "ll")
	bwmethod = ifelse(estimator == "cvls", "cv.ls", "cv.aic")
	bwmethodcat = list("cvml"="cv.ml", "cvls"="cv.ls", "normalref"="normal-reference")[[estimatorcat]]
	bwscaling = ifelse(scaling == "scale", TRUE, FALSE)
	bwtype = list("fixed"="fixed", "generalizednn"="generalized_nn", "adaptivenn"="adaptive_nn")[[cvtype]]
	if (bootmethod == "wildrademacher") {
		bootmethod = "wild-rademacher"
	}

	allvars = c(dep, indepnp)

	itmax = maxiter
	ftol = valtolerance
	tol = loctolerance
	
	dta = tryCatch(spssdata.GetDataFromSPSS(allvars, row.label=id, missingValueToNA=TRUE,
		factorMode="labels"), error=function(e) {
      warns$warn(e, dostop=TRUE)
		}
  )
	# remove any factor levels that have a zero count
	# NB: this can't reattach labels to string factors
	dta = cleaned(dta)
	allargs["factorcountnp"] = sum(unlist(lapply(dta[1,2:(nindepnp+1)], is.factor)))
	allargs["ordinalcountnp"] = sum(unlist(lapply(dta[1,2:(nindepnp+1)], is.ordered)))
	if (is.factor(dta[,1])) {
		allargs["estimator"] = estimatorcat
		bwmethod = bwmethodcat
		if (estimatorcat == "normalref" && (cvtype != "fixed" ||
			cvkerneltype != "gaussian")) {
			cvtype = "fixed"
			bwtype="fixed"
			allargs["cvtype"] = cvtype
			cvkerneltype = "gaussian"
			allargs["cvkerneltype"] = cvkerneltype
			warns$warn(gtxt("Normal reference requires a fixed, Gaussian kernel.  Settings have been adjusted"), dostop=FALSE)
		}
	}
	if (!is.null(dataset) && is.factor(dta[,1]) && resids) {
		allargs["resids"] = FALSE
		warns$warn(gtxt("Residuals are not available for categorical dependent variables"),
			dostop=FALSE)
	}
	
	# nonparametric regression
	options(np.messages=FALSE)
	if (!is.factor(dta[,1])) {
		bw = tryCatch({
			npregbw(formula=as.formula(frml),data=dta, na.action=naaction,
				regtype = regtype, bwmethod=bwmethod, bwscaling=bwscaling, 
				bwtype=bwtype, ckertype=cvkerneltype, ckerorder=cvorder,
				ukertype=ukerneltype, okertype=okerneltype, nmulti=nstart,
				ftol=ftol, tol=tol)
			},
			error = function(e) {
				warns$warn(e$message, dostop=TRUE)
			}
		)
			
		res = tryCatch(
      npreg(bws=bw, data=dta, gradients=(plots == "gradient")),
			error=function(e) {
        warns$warn(e$message, dostop=TRUE)
			}
		)
	} else {
		# npcdensbw will not accept a variable holding either a formula string
		# or an as.formula value, but it will work with do.call :-(
		arglist = list(formula=as.formula(frml), data=dta, na.action=naaction,
			bwmethod=bwmethod,
			bwtype=bwtype, ckertype=cvkerneltype, ckerorder=cvorder,
			ukertype=ukerneltype, okertype=okerneltype, nmulti=nstart,
			ftol=ftol, tol=tol
		)

		bw = tryCatch({
			do.call(npcdensbw, arglist)
			},
			error = function(e) {
				warns$warn(e$message, dostop=TRUE)
			}
		)
			
		res = tryCatch(
      npconmode(bws=bw, data=dta, gradients=(plots == "gradient")),
			error=function(e) {
        warns$warn(e$message, dostop=TRUE)
			}
		)
	}
	sigres = NULL
	if (sigtests) {
		if (class(bw) != "rbandwidth") {
			warns$warn(gtxt("Significance tests are not available for a categorical dependent variable"),
				dostop=FALSE)
			sigtests = FALSE
		} else {
			sigindex = min(nvars, nindepnp)
			allargs['sigindex'] = sigindex
			sigres = tryCatch({
				npsigtest(bw, xdat=dta[,2:(nindepnp+1)], ydat=dta[,1], data=dta,
				index=1:sigindex, na.action=naaction,
				boot.method=bootmethod, boot.num=bootreps, random.seed=rnseed)
			},
			error = function(e) {
				warns$warn(e$message, dostop=TRUE)
			}
			)
		}
	}
	

  # print results
  browser()
  displayresults(allargs, allvars, dta, bw, res, sigres, warns)
  if (!is.null(dataset)) {
	  createdataset(res, dta, allargs, warns)
  }

    # clean up workspace
    res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
}



displayresults = function(info, allvars, dta, bw, res, sigres, warns) {

    StartProcedure(gtxt("Nonparametric Regression"), "STATSNPREG")
    
    # summary results
	# input specifications
	lbls = c(gtxt("Dependent Variable"),
		gtxt("Number of Nonparametric Scale Independent Variables"),
		gtxt("Number of Nonparametric Ordinal Independent Variables"),
		gtxt("Number of Nominal Independent Variables"),
		gtxt("Kernel Regression Estimator"),
		gtxt("Bandwidth Selection Method"),
		gtxt("Continuous Variable Bandwidth Type"),
		gtxt("Bandwidth Scaling"),
		gtxt("Continuous Kernel Type"),
		gtxt("Continuous Kernel Order"),
		gtxt("Ordinal Kernel Type"),
		gtxt("Nominal Kernel Type"),
		gtxt("Number of Calculation Starts"),
		gtxt("Cross-Validation Value Tolerance"),
		gtxt("Cross-Validation Location Tolerance"),
		gtxt("Significance Bootstrap Method"),
		gtxt("Missing Value Treatment"),
		gtxt("Output Dataset")
	)
    vals = c(info["dep"],
		info[["nindepnp"]] - info[["factorcountnp"]],
		info["ordinalcountnp"],
		info[["factorcountnp"]] - info[["ordinalcountnp"]],
		info["regtype"],
		ifelse(info[["estimator"]] =="cvls", gtxt("least squares cross validation"), gtxt("Kullback-Leibler cross validation")),
		info[["cvtype"]],
		info["scaling"],
		info["cvkerneltype"],
		info["cvorder"],
		info["okerneltype"],
		info["ukerneltype"],
		info["nstart"],
		info["valtolerance"],
		info["loctolerance"],
		list("iid"=gtxt("iid"), "wild"="wild", "wildrademacher"=gtxt("wild with Rademacher"))[info[["bootmethod"]]],
		info["missingv"],
		ifelse(is.null(info[["dataset"]]), gtxt("--NA--"), info[["dataset"]])
    )

  # settings and result summary
  spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Settings"),
      collabels=c(gtxt("Settings")), templateName="NPSETTINGS", outline=gtxt("Settings"))
     
	scaption = gtxt("Computations done by R package np")
	# bandwidth
  dvcentral = getdepvarstat(dta, bw[['nobs.omit']])
  if (is.factor(dta[[1]])) {
    caption = gtxt("If the dependent variable has multiple modes, the first is reported")
  } else {
    caption = NULL
  }
	lbls = c(gtxt("Number of Cases"),
		gtxt("Valid Cases"),
		gtxt("Omitted Cases"),
    ifelse(is.factor(dta[[1]]), gtxt("Dependent Variable Mode"), gtxt("Dependent Variable Median")),
		gtxt("Objective Function"),
		gtxt("R-Squared"),
		gtxt("Standard Error"),
		gtxt("Mean Absolute Error"),
		gtxt("Correct Classification Percentage")
	)
	if (!is.null(res$ybw)) {
		lbls[length(lbls)+1] = gtxt("Dependent Variable Bandwidth")
	}

	vals = c(bw[['nobs']] + bw[['nobs.omit']],
		bw[['nobs']],
		bw[['nobs.omit']],
    dvcentral,
		round(bw[['fval']], 5),
		ifelse(is.null(res$R2), gtxt("--NA--"), round(res$R2,5)),
		ifelse(is.null(res$MSE), gtxt("--NA--"), round(sqrt(res$MSE),5)),
		ifelse(is.null(res$MAE), gtxt("--NA--"), round(sqrt(res$MAE),5)),
		ifelse(is.null(res$confusion.matrix), gtxt("--NA--"), 
			sum(diag(res$confusion.matrix)) / sum(res$confusion.matrix) * 100)
	)
	if (!is.null(res$ybw)) {
		vals[length(vals)+1] = res$ybw
	}
  spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Bandwidth Summary"),
      collabels=c(gtxt("Statistics")), templateName="NPBWSUMMARY", outline=gtxt("Bandwidth Summary"),
      caption=caption)
	# bandwidth parameters and significance test
	df = buildbwinfo(bw)

	sigcol = rep(NA, info[["nindepnp"]])
	caption = ""
	if (!is.null(sigres)) {
		for (i in 1:length(sigres$P)) {
			sigcol[i] = sigres$P[[i]]
		}
		caption = sprintf(gtxt("Significance test type: %s.  Method: %s.  Bootstrap replications: %s"),
			sigres$ptype, sigres$pmethod, sigres$boot.num)
	}
	df = data.frame(df, sigcol)
	collabels = c(gtxt("Bandwidth"), gtxt("Scale Factor or Lambda"), gtxt("Scale"), gtxt("Sig."))
	spsspivottable.Display(df, title=gtxt("Bandwidth"),
		collabels=collabels, templateName="NPBW",
		outline=gtxt("Bandwidth Parameters"),
		caption=caption
	)

	# confusion matrix for categorical dependent
	if (!is.null(res$confusion.matrix)) {
		condf = data.frame(res$confusion.matrix[, 1:ncol(res$confusion.matrix)])
		names(condf) = dimnames(res$confusion.matrix)$Predicted
		row.names(condf) = dimnames(res$confusion.matrix)$Actual
		spsspivottable.Display(condf, title=gtxt("Confusion Matrix"), 
			templateName="NPCONFUSION",
			outline=gtxt("Confusion Matrix"),
			rowdim=gtxt("Actual"), hiderowdimtitle=FALSE,
			coldim=gtxt("Predicted"), hidecoldimtitle=FALSE,
			format = formatSpec.Count,
      caption = gtxt("Predicted categories with zero count are not shown")
		)
	}

	# plots
	if (!info[['plots']] == "none") {
		plotter(bw, dta, info[['dep']], info[['indepnp']], info[['includepts']], 
			info[['errorbnds']], info[['plots']], warns)
	}
    warns$display(inproc=TRUE)
}

plotter = function(bw, dta, y, x, pts, errorbnds, plots, warns) {
	# plot results each as a separate plot
	# bw is the bandwidth object
	# y and x are the variable names
	# errorbnds is the error bounds setting
	
	# api does not exist in older Statistics versions
	rc=tryCatch(spssRGraphics.SetGraphicsLabel(gtxt("Nonparametric Plot")), error=function(e) {})
	if (class(bw) == "conbandwidth" && plots == "gradient") {
#     if (length(x) >= 4) {  # scale plot to accommodate more variables
#   	  tryCatch({
#         browser()
#   	    temp = paste(tempdir(), "res.png", sep="/")
#   	    size = min(7, length(x)+1) * 100
#   	    png(temp, width=size, height=size, units="px")
#         plot(bw, gradients=(plots=="gradient"))
#         dev.off()
#         spssRGraphics.Submit(temp)
#   	  },
#         error=function(e) {
#           warns$warn(gtxt("Gradient plot coud not be produced"), dostop=FALSE)
#         }
#   	  )
# 	  } else {
      browser()
  		tryCatch(plot(bw, gradients=(plots=="gradient")),
  			error=function(e) {
  				warns$warn(gtxt("Gradient plot could not be produced.  There may be too many variables"),		
            dostop=FALSE)
  			}
  		)
#	  }
	} else {
		if (class(bw) == "conbandwidth") {
			pts = FALSE
			ylab = gtxt("Conditional Density")  # y axis is a probability for these charts
		} else {
			ylab = y
		}
		#get plot data
		pd = plot(bw, plot.errors.method=errorbnds, gradients=(plots=="gradient"), 
			plot.behavior="data", perspective=FALSE)
		thenames = names(pd)
		for (v in 1:length(thenames)) {
			if (regexpr("^r|^cd", thenames[[v]]) == -1) {    # neither response nor gradient
				next
			}
			isgrad = regexpr("^rg", thenames[[v]]) > 0   # indicates gradient
			if (isgrad) {
				if (!(plots == "gradient")) {
					next
				}
				ctitle = gtxt("Gradient")
			} else {
				if (!(plots == "response")) {
					next
				}
				ctitle = gtxt("Response")
			}
			item = pd[[v]]
			if (isgrad) {
				xdata = item$eval[,v]
				plot(xdata, gradients(item), lwd=2, type = "l", xlab=x[v], ylab=y, main=ctitle, perspective=FALSE)
			} else {
				fit = fitted(item)
				if (class(item) == "condensity") {
					lowerci = NULL
					upperci = NULL
					selist = NULL
					fit = item$condens
					ylim = NULL
					if (v < length(thenames)) {  # an x variable
						xdata = item$xeval[,v]
						xlab = x[v]
					} else {  #dep variable
						xdata = item$yeval
						xlab = y
					}
				} else {
					selist = se(item)
					lowerci = fit + selist[,1]
					upperci = fit + selist[,2]
					xdata = item$eval[,v]
					xlab = x[v]
					if (pts) {
						ylim = c(min(dta[,1], lowerci), max(dta[,1], upperci))
					} else {
						ylim = c(min(fit, lowerci), max(fit, upperci))
					}
				}
				suppressWarnings(plot(xdata, fit, lwd=2, type = "l", ylim=ylim, xlab=xlab, ylab=ylab, main=ctitle, perspective=FALSE))
				if (pts) {
					points(dta[,v+1], dta[,1], col="gray")
				}
				if (class(bw) != "conbandwidth" && nrow(selist) > 0) {
					lines(xdata, lowerci, lty=3, lwd=2)
					lines(xdata, upperci, lty=3, lwd=2)
				}
			}
		}
	}
	###rc=tryCatch(spssRGraphics.SetGraphicsLabel(), error=function(e) {})
}

createdataset = function(res, dta, details, warns) {
  # Create residuals and/or predicted values dataset
  # Dataset name is known to be okay, and procedure state is ended
	# res is the bandwidth results
	# dta is the data
	# details is all the input arguments

	if (class(res) == "conbandwidth") {
		details["resids"] = FALSE   # residuals are not available
	}
  if (!is.null(details[["dataset"]])) {
	if (!details[["predicted"]] && !details[["resids"]]) {
		stop(gtxt("An output dataset was specified but no variables were requested or none are available"), 
      call.=FALSE)
	}
	# construct a data frame with the requested variables
	# If the dataset has missing values, the predict and residuals output will have
	# No values for those, but they will carry the ID value, so we have to
	# pick up the ID values from the row names in the generated data frame and
	# copy them to the ID variable for sending back to Statistics
	# n.b. The doc for predict and residuals claims that these values are
	# present but NA, but the facts are otherwise.

	if (details[["predicted"]]) {
		if (class(res) == "conbandwidth") {
			pred = res$condens
		} else {
			pred = fitted(res) 
		}
		theframe = data.frame(pred)
	}
	# residuals flag is turned off above for categorical dependent var models
	if (details[["resids"]]) {
		resids = residuals(res)
		if (!details[["predicted"]]) {
			theframe = data.frame(resids)
		} else {
			theframe[2] = resids
		}
	}
	theframe = data.frame(id=row.names(dta), theframe)
  if (!is.null(details[["id"]])) {  # was an id variable provided
    vardict = spssdictionary.GetDictionaryFromSPSS(details["id"])
    loc = match(details[["id"]], vardict["varName",])
    idtype = as.integer(vardict[["varType", loc]])
    idformat = vardict[["varFormat", loc]]
    idlabel = vardict[["varLabel", loc]]
  	if (idlabel == "") {  # if no label, use the id variable name as the label
  		idlabel = details[["id"]]
  	}
  } else {
      idtype = 0
      idformat = "F10.0"
      idlabel = ""
  }
	dictspec = list(c("ID", idlabel, idtype, idformat, "nominal"))
	if (details[["predicted"]]) {
		if (class(res) == "conmode") {
			vname = "Density"
			vlabel = gtxt("Conditional Density")
		} else {
			vname = "PredictedValues"
			vlabel = "Predicted Values"
		}
		dictspec[2] = list(c(vname, vlabel,
			0, "F8.2", "scale"))
	}
	if (details[["resids"]]) {
		dictspec[length(dictspec)+1] = list(c("Residuals", gtxt("Residuals"),
			0, "F8.2", "scale"))
	}

    dict = spssdictionary.CreateSPSSDictionary(dictspec)
    spssdictionary.SetDictionaryToSPSS(details[["dataset"]], dict)
    spssdata.SetDataToSPSS(details[["dataset"]], theframe)
    spssdictionary.EndDataStep()
  }
}

# This function allows us to use labels as the data option
# Empty levels in a factor confuse np
cleaned = function(dta) {
	# remove any levels that have a zero count
	# dta is a data frame
	
	# String variable values seem to come over unlabelled but
	# with levels including labels and values.  Can't fix that here
	
	prune = function(v) {
		if (is.factor(v)) {
			return(v[,drop=TRUE])
		} else {
			return(v)
		}
	}
	dta = data.frame(lapply(dta, prune), row.names=row.names(dta))
	return(dta)
}

getdepvarstat = function(dta, omitted) {
  # return dependent variable median or mode
  # dta is the data frame with dependent variable first
  # omitted is the count of omitted cases

  if (is.factor(dta[[1]])) { # categorical dependent
    if (omitted == 0) {
      return(levels(dta[[1]])[which.max(table(dta[1]))])
    } else { # find (first) largest factor level count among complete cases
        return(levels(dta[[1]])[which.max(table(dta[complete.cases(dta),][1]))])
    }
  } else { # continuous dependent
        return(median(dta[[1]], na.rm=TRUE))
    }
}

buildfrml = function(y, xnp) {
	# return formula
	# y is the dependent variable, xnp is the nonparametric predictors
	x = paste(xnp, sep="", collapse="+")
	return(paste(y, x, sep="~"))
}
	
buildbwinfo = function(bw) {
	# return data frame of bandwidth information
	
	types = getScaleInfo(bw)
	df = data.frame(cbind(round(bw$bandwidth$x, 5), types, round(bw$sumNum$x, digits=5)), row.names=bw$xnames)
	return(df)
}

getScaleInfo = function(bwobj) {
	# return vector of scale type info
	# parsing the text output from private function in np :-( :-(

	tryCatch({
		lbls = list("Lambda Max", "Scale Factor")
		qqq = np:::genBwScaleStrs(bwobj)
		qqqsplit = strsplit(qqq, "\n")[[1]]
		types = list()
		for (i in 2:length(qqqsplit)) {
			if (length(grep("Scale Factor:", qqqsplit[[i]], ignore.case=TRUE)) > 0) {
				types[i-1] = lbls[[2]]
			} else if (length(grep("Lambda Max:", qqqsplit[[i]], ignore.case=TRUE)) > 0){
				types[i-1] = lbls[[1]]
			} else {
				types[i-1] = ""
			}
		}
		return(types)
		}, error = function(bwobj) {   # in case we can't get the type information			
			return(rep("", length(bwobj$xnames)))
		}
	)
}
	
dscheck = function(alldsspecs, warns) {
	# check dataset validation conditions
	
  if (!is.null(alldsspecs)) {
    alldatasets = spssdata.GetDataSetList()
    if ("*" %in% alldatasets) {
        warns$warn(gtxt("The active dataset must have a name if creating new datasets"),
            dostop=TRUE)
    }
    if (length(intersect(alldsspecs, alldatasets) > 0)) {
        warns$warn(gtxt("One or more specified output dataset names are already in use"),
            dostop=TRUE)
    }
  }
}

# localization initialization
setuplocalization = function(domain) {
  # find and bind translation file names
  # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
  # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

  fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
  browser()
  bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
  if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
      spsspkg.StartProcedure(procname, omsid)
  }
  else {
      spsspkg.StartProcedure(omsid)
  }
}
gtxt <- function(...) {
    return(gettext(...,domain="STATS_NONPAR_REGR"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_NONPAR_REGR"))
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

    if (is.null(msg) || dostop) {
        lcl$display(inproc)  # display messages and end procedure state
        if (dostop) {
            stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
        }
    }
}

    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

    if (lcl$msgnum == 0) {   # nothing to display
        if (inproc) {
            spsspkg.EndProcedure()
        }
    } else {
        if (!inproc) {
            procok =tryCatch({
                StartProcedure(lcl$procname, lcl$omsid)
                TRUE
                },
                error = function(e) {
                    FALSE
                }
            )
        } else {
			procok = TRUE
		}
        if (procok) {  # build and display a Warnings table if we can
            table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
            rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

    for (i in 1:lcl$msgnum) {
        rowcategory = spss.CellText.String(as.character(i))
        BasePivotTable.SetCategories(table,rowdim,rowcategory)
        BasePivotTable.SetCellValue(table,rowcategory, 
            spss.CellText.String(lcl$msglist[[i]]))
    }
    spsspkg.EndProcedure()   # implies display
} else { # can't produce a table
    for (i in 1:lcl$msgnum) {
        print(lcl$msglist[[i]])
    }
}
}
}
return(lcl)
}

Run<-function(args){
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
    spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep", islist=FALSE),
		spsspkg.Template("INDEP", subc="", ktype="existingvarlist", var="indepnp", islist=TRUE),
		spsspkg.Template("ID", subc="", ktype="existingvarlist", var="id", islist=FALSE),
		
    spsspkg.Template("REGTYPE", subc="BANDWIDTH", ktype="str", var="regtype",
			vallist = list("localconstant", "locallinear")),
		spsspkg.Template("ESTIMATORCONT", subc="BANDWIDTH", ktype="str", var="estimator",
			vallist = list("cvls", "cvaic")),
		spsspkg.Template("ESTIMATORCAT", subc="BANDWIDTH", ktype="str", var="estimatorcat",
			vallist = list("cvml", "cvls", "normalref")),
    spsspkg.Template("SCALING", subc="BANDWIDTH",  ktype="str", 
        var="scaling", vallist=list("scale", "raw")),    # must be mapped to TRUE/FALSE
    spsspkg.Template("CVTYPE", subc="BANDWIDTH",  ktype="str", 
        var="cvtype", vallist=list("fixed", "generalizednn", "adaptivenn")),
    spsspkg.Template("CVKERNELTYPE", subc="BANDWIDTH", ktype="str",
        var="cvkerneltype", vallist=list("gaussian", "epanechnikov", "uniform")),
    spsspkg.Template("CVORDER", subc="BANDWIDTH", ktype="int", 
        var="cvorder"),   #2,4,6,8
		spsspkg.Template("UKERNELTYPE", subc="BANDWIDTH", ktype="str", 
        var="ukerneltype", vallist=list("aitchisonaitken", "liracine")),
    spsspkg.Template("OKERNELTYPE", subc="BANDWIDTH", ktype="str", var="okerneltype",
			vallist=list("wangvanryzin", "liracine")),
		spsspkg.Template("NSTART", subc="BANDWIDTH", ktype="int", var="nstart",
			vallist=list(1)),

    spsspkg.Template("MISSING", subc="OPTIONS", ktype="str", 
      var="missingv", vallist = list("omit", "stop")),
		spsspkg.Template("MAXITER", subc="OPTIONS", ktype="int", 
      var="maxiter", vallist = list(1)), 
    spsspkg.Template("VALTOLERANCE", subc="OPTIONS", ktype="float", 
      var="valtolerance", vallist=list(1e-20)),
    spsspkg.Template("LOCTOLERANCE", subc="OPTIONS", ktype="float", 
      var="loctolerance", vallist=list(1e-20)),		
        
    spsspkg.Template("SIGTESTS", subc="OUTPUT", ktype="bool", var="sigtests"),
    spsspkg.Template("NVARS", subc="OUTPUT", ktype="int", var="nvars"),
		spsspkg.Template("BOOTMETHOD", subc="OUTPUT", ktype="str", var="bootmethod",
			vallist=list("iid", "wild", "wildrademacher")),
		spsspkg.Template("BOOTREPS", subc="OUTPUT", ktype="int", var="bootreps"),
		spsspkg.Template("RNSEED", subc="OUTPUT", ktype="float", var="rnseed"),
		spsspkg.Template("PLOTS", subc="OUTPUT", ktype="str", var="plots",
			vallist=list("response", "gradient", "none")),
		spsspkg.Template("INCLUDEPTS", subc="OUTPUT", ktype="bool", var="includepts"),
		spsspkg.Template("ERRORBNDS", subc="OUTPUT", ktype="str", var="errorbnds",
			vallist=list("none", "asymptotic", "bootstrap")),
			
		spsspkg.Template("DATASET", "SAVE", ktype="varname", var="dataset"),
		spsspkg.Template("RESIDUALS", "SAVE", ktype="bool", var="resids"),
		spsspkg.Template("PREDICTED", "SAVE", ktype="bool", var="predicted")	
    ))        
  if ("HELP" %in% attr(args,"names"))
      writeLines(helptext)
  else
      res <- spsspkg.processcmd(oobj,args,"npext")
}

# 
# 
# #################
# npcdensbw <-
#   function(...){
#     args = list(...)
#     if (is(args[[1]],"formula"))
#       UseMethod("npcdensbw",args[[1]])
#     else if (!is.null(args$formula))
#       UseMethod("npcdensbw",args$formula)
#     else
#       UseMethod("npcdensbw",args[[which(names(args)=="bws")[1]]])
#   }
# 
# npcdensbw.formula <-
#   function(formula, data, subset, na.action, call, ...){
# 
#     mf <- match.call(expand.dots = FALSE)
#     m <- match(c("formula", "data", "subset", "na.action"),
#                names(mf), nomatch = 0)
#     mf <- mf[c(1,m)]
#     
#     if(!missing(call) && is.call(call)){
#       ## rummage about in the call for the original formula
#       for(i in 1:length(call)){
#         if(tryCatch(class(eval(call[[i]])) == "formula",
#                     error = function(e) FALSE))
#           break;
#       }
#       mf[[2]] <- call[[i]]
#       
#     }
#     
#     
#     mf[[1]] <- as.name("model.frame")
#     
#     variableNames <- explodeFormula(mf[["formula"]])
#     
#     ## make formula evaluable, then eval
#     varsPlus <- lapply(variableNames, paste, collapse=" + ")
#     mf[["formula"]] <- as.formula(paste(" ~ ", varsPlus[[1]]," + ",
#                                         varsPlus[[2]]),
#                                   env = environment(formula))
#     
#     mf <- eval(mf, parent.frame())
#     
#     ydat <- mf[, variableNames[[1]], drop = FALSE]
#     xdat <- mf[, variableNames[[2]], drop = FALSE]
#     
#     tbw = npcdensbw(xdat = xdat, ydat = ydat, ...)
#     
#     ## clean up (possible) inconsistencies due to recursion ...
#     tbw$call <- match.call(expand.dots = FALSE)
#     environment(tbw$call) <- parent.frame()
#     tbw$formula <- formula
#     tbw$rows.omit <- as.vector(attr(mf,"na.action"))
#     tbw$nobs.omit <- length(tbw$rows.omit)
#     tbw$terms <- attr(mf,"terms")
#     tbw$variableNames <- variableNames
#     
#     tbw
#   }
# 
# npcdensbw.conbandwidth <- 
#   function(xdat = stop("data 'xdat' missing"),
#            ydat = stop("data 'ydat' missing"),
#            bws, bandwidth.compute = TRUE,
#            auto = TRUE,
#            nmulti, remin = TRUE, itmax = 10000,
#            ftol=1.19209e-07, tol=1.49012e-08, small=2.22045e-16,
#            ...){
#     
# 
#     ydat = toFrame(ydat)
#     xdat = toFrame(xdat)
#     
#     if (missing(nmulti)){
#       nmulti <- min(5,(dim(ydat)[2]+dim(xdat)[2]))
#     }
#     
#     if (length(bws$ybw) != dim(ydat)[2])
#       stop(paste("length of bandwidth vector does not match number of columns of", "'ydat'"))
#     
#     if (length(bws$xbw) != dim(xdat)[2])
#       stop(paste("length of bandwidth vector does not match number of columns of", "'xdat'"))
#     
#     if (dim(ydat)[1] != dim(xdat)[1])
#       stop(paste("number of rows of", "'ydat'", "does not match", "'xdat'"))
#     
#     yccon = unlist(lapply(as.data.frame(ydat[,bws$iycon]),class))
#     if ((any(bws$iycon) && !all((yccon == class(integer(0))) | (yccon == class(numeric(0))))) ||
#           (any(bws$iyord) && !all(unlist(lapply(as.data.frame(ydat[,bws$iyord]),class)) ==
#                                     class(ordered(0)))) ||
#           (any(bws$iyuno) && !all(unlist(lapply(as.data.frame(ydat[,bws$iyuno]),class)) ==
#                                     class(factor(0)))))
#       stop(paste("supplied bandwidths do not match", "'ydat'", "in type"))
#     
#     xccon = unlist(lapply(as.data.frame(xdat[,bws$ixcon]),class))
#     if ((any(bws$ixcon) && !all((xccon == class(integer(0))) | (xccon == class(numeric(0))))) ||
#           (any(bws$ixord) && !all(unlist(lapply(as.data.frame(xdat[,bws$ixord]),class)) ==
#                                     class(ordered(0)))) ||
#           (any(bws$ixuno) && !all(unlist(lapply(as.data.frame(xdat[,bws$ixuno]),class)) ==
#                                     class(factor(0)))))
#       stop(paste("supplied bandwidths do not match", "'xdat'", "in type"))
#     
#     ## catch and destroy NA's
#     goodrows <- 1:dim(xdat)[1]
#     rows.omit <- unclass(na.action(na.omit(data.frame(xdat,ydat))))
#     goodrows[rows.omit] <- 0
#     
#     if (all(goodrows==0))
#       stop("Data has no rows without NAs")
#     
#     xdat = xdat[goodrows,,drop = FALSE]
#     ydat = ydat[goodrows,,drop = FALSE]
#     
#     
#     nrow = nrow(ydat)
#     yncol = ncol(ydat)
#     xncol = ncol(xdat)
#     
#     ## at this stage, data to be sent to the c routines must be converted to
#     ## numeric type.
#     
#     ydat = toMatrix(ydat)
#     
#     yuno = ydat[, bws$iyuno, drop = FALSE]
#     ycon = ydat[, bws$iycon, drop = FALSE]
#     yord = ydat[, bws$iyord, drop = FALSE]
#     
#     
#     xdat = toMatrix(xdat)
#     
#     xuno = xdat[, bws$ixuno, drop = FALSE]
#     xcon = xdat[, bws$ixcon, drop = FALSE]
#     xord = xdat[, bws$ixord, drop = FALSE]
#     
#     tbw <- bws
#     
#     if (bandwidth.compute){
#       myopti = list(num_obs_train = nrow,
#                     iMultistart = ifelse(nmulti==0,IMULTI_FALSE,IMULTI_TRUE),
#                     iNum_Multistart = nmulti,
#                     int_use_starting_values = ifelse(all(bws$ybw==0) && all(bws$xbw==0),
#                                                      USE_START_NO, USE_START_YES),
#                     int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
#                     BANDWIDTH_den_extern = switch(bws$type,
#                                                   fixed = BW_FIXED,
#                                                   generalized_nn = BW_GEN_NN,
#                                                   adaptive_nn = BW_ADAP_NN),
#                     itmax=itmax, int_RESTART_FROM_MIN=ifelse(remin,RE_MIN_TRUE,RE_MIN_FALSE), 
#                     int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE), 
#                     bwmethod = switch(bws$method,
#                                       cv.ml = CBWM_CVML,
#                                       cv.ls = CBWM_CVLS,
#                                       cv.ls.np = CBWM_NPLS),        
#                     xkerneval = switch(bws$cxkertype,
#                                        gaussian = CKER_GAUSS + bws$cxkerorder/2 - 1,
#                                        epanechnikov = CKER_EPAN + bws$cxkerorder/2 - 1,
#                                        uniform = CKER_UNI),
#                     ykerneval = switch(bws$cykertype,
#                                        gaussian = CKER_GAUSS + bws$cykerorder/2 - 1,
#                                        epanechnikov = CKER_EPAN + bws$cykerorder/2 - 1,
#                                        uniform = CKER_UNI),
#                     uxkerneval = switch(bws$uxkertype,
#                                         aitchisonaitken = UKER_AIT,
#                                         liracine = UKER_LR),
#                     uykerneval = switch(bws$uykertype,
#                                         aitchisonaitken = UKER_AIT,
#                                         liracine = UKER_LR),
#                     oxkerneval = switch(bws$oxkertype,
#                                         wangvanryzin = OKER_WANG,
#                                         liracine = OKER_LR),
#                     oykerneval = switch(bws$oykertype,
#                                         wangvanryzin = OKER_WANG,
#                                         liracine = OKER_LR),
#                     ynuno = dim(yuno)[2],
#                     ynord = dim(yord)[2],
#                     yncon = dim(ycon)[2],
#                     xnuno = dim(xuno)[2],
#                     xnord = dim(xord)[2],
#                     xncon = dim(xcon)[2],
#                     fast = FALSE,
#                     auto = auto)
#       
#       myoptd = list(ftol=ftol, tol=tol, small=small)
#       
#       if (bws$method != "normal-reference"){
#         myout=
#           .C("np_density_conditional_bw", as.double(yuno), as.double(yord), as.double(ycon),
#              as.double(xuno), as.double(xord), as.double(xcon),
#              as.integer(myopti), as.double(myoptd), 
#              bw = c(bws$xbw[bws$ixcon],bws$ybw[bws$iycon],
#                     bws$ybw[bws$iyuno],bws$ybw[bws$iyord],
#                     bws$xbw[bws$ixuno],bws$xbw[bws$ixord]),
#              fval = double(2),
#              PACKAGE="np" )[c("bw","fval")]
#       } else {
#         nbw = double(yncol+xncol)
#         gbw = bws$yncon+bws$xncon
#         if (gbw > 0){
#           nbw[1:gbw] = (4/3)^0.2
#           if(!bws$scaling)
#             nbw[1:gbw]=nbw[1:gbw]*EssDee(data.frame(xcon,ycon))*nrow^(-1.0/(2.0*bws$cxkerorder+gbw))
#         }
#         myout= list( bw = nbw, fval = c(NA,NA) )
#       }
#       
#       yr = 1:yncol
#       xr = 1:xncol
#       rorder = numeric(yncol + xncol)
#       
#       ## bandwidths are passed back from the C routine in an unusual order
#       ## xc, y[cuo], x[uo]
#       
#       rxcon = xr[bws$ixcon]
#       rxuno = xr[bws$ixuno] 
#       rxord = xr[bws$ixord] 
#       
#       rycon = yr[bws$iycon] 
#       ryuno = yr[bws$iyuno] 
#       ryord = yr[bws$iyord] 
#       
#       
#       ## rorder[c(rxcon,rycon,ryuno,ryord,rxuno,rxord)]=1:(yncol+xncol)
#       
#       tbw <- bws
#       tbw$ybw[c(rycon,ryuno,ryord)] <- myout$bw[yr+bws$xncon]
#       tbw$xbw[c(rxcon,rxuno,rxord)] <- myout$bw[setdiff(1:(yncol+xncol),yr+bws$xncon)]
#       
#       tbw$fval = myout$fval[1]
#       tbw$ifval = myout$fval[2]
#     }
#     
#     ## bandwidth metadata
#     tbw$sfactor <- tbw$bandwidth <- list(x = tbw$xbw, y = tbw$ybw)
#     
#     bwf <- function(i){
#       tbw$bandwidth[[i]][tl[[i]]] <<- (tbw$bandwidth[[i]])[tl[[i]]]*dfactor[[i]]
#     }
#     
#     sff <- function(i){
#       tbw$sfactor[[i]][tl[[i]]] <<- (tbw$sfactor[[i]])[tl[[i]]]/dfactor[[i]]
#     }
#     
#     myf <- if(tbw$scaling) bwf else sff
#     
#     if ((tbw$xnuno+tbw$ynuno) > 0){
#       dfactor <- nrow^(-2.0/(2.0*tbw$cxkerorder+tbw$ncon))
#       dfactor <- list(x = dfactor, y = dfactor)
#       
#       tl <- list(x = tbw$xdati$iuno, y = tbw$ydati$iuno)
#       
#       lapply(1:length(tl), myf)
#     }
#     
#     if ((tbw$xnord+tbw$ynord) > 0){
#       dfactor <- nrow^(-2.0/(2.0*tbw$cxkerorder+tbw$ncon))
#       dfactor <- list(x = dfactor, y = dfactor)
#       
#       tl <- list(x = tbw$xdati$iord, y = tbw$ydati$iord)
#       
#       lapply(1:length(tl), myf)
#     }
#     
#     
#     if (tbw$ncon > 0){
#       dfactor <- nrow^(-1.0/(2.0*tbw$cxkerorder+tbw$ncon))
#       dfactor <- list(x = EssDee(xcon)*dfactor, y = EssDee(ycon)*dfactor)
#       
#       tl <- list(x = tbw$xdati$icon, y = tbw$ydati$icon)
#       
#       lapply(1:length(tl), myf)
#     }
#     
#     tbw <- conbandwidth(xbw = tbw$xbw,
#                         ybw = tbw$ybw,
#                         bwmethod = tbw$method,
#                         bwscaling = tbw$scaling,
#                         bwtype = tbw$type,
#                         cxkertype = tbw$cxkertype,
#                         cxkerorder = tbw$cxkerorder,
#                         uxkertype = tbw$uxkertype,
#                         oxkertype = tbw$oxkertype,
#                         cykertype = tbw$cykertype,
#                         cykerorder = tbw$cykerorder,
#                         uykertype = tbw$uykertype,
#                         oykertype = tbw$oykertype,
#                         fval = tbw$fval,
#                         ifval = tbw$ifval,
#                         nobs = tbw$nobs,
#                         xdati = tbw$xdati,
#                         ydati = tbw$ydati,      
#                         xnames = tbw$xnames,
#                         ynames = tbw$ynames,
#                         sfactor = tbw$sfactor,
#                         bandwidth = tbw$bandwidth,
#                         rows.omit = rows.omit,
#                         bandwidth.compute = bandwidth.compute)
#     
#     tbw
#   }
# 
# nnpcdensbw.NULL <-
#   function(xdat = stop("data 'xdat' missing"),
#            ydat = stop("data 'ydat' missing"),
#            bws, ...){
#     
#     ## maintain x names and 'toFrame'
# 
#     xdat <- toFrame(xdat)
#     
#     ## maintain y names and 'toFrame'
#     ydat <- toFrame(ydat)
#     
#     ## do bandwidths
#     
#     bws = double(ncol(ydat)+ncol(xdat))
#     
#     tbw <- npcdensbw.default(xdat = xdat, ydat = ydat, bws = bws, ...)
#     
#     ## clean up (possible) inconsistencies due to recursion ...
#     mc <- match.call(expand.dots = FALSE)
#     environment(mc) <- parent.frame()
#     tbw$call <- mc
#     
#     tbw
#   }
# 
# npcdensbw.default <-
#   function(xdat = stop("data 'xdat' missing"),
#            ydat = stop("data 'ydat' missing"),
#            bws, 
#            bandwidth.compute = TRUE,
#            auto, nmulti, remin, itmax,
#            ftol, tol, small,
#            ## dummy arguments for conbandwidth() function call
#            bwmethod, bwscaling, bwtype,
#            cxkertype, cxkerorder,
#            cykertype, cykerorder,
#            uxkertype, uykertype,
#            oxkertype, oykertype,
#            ...){
#     
#     ## maintain x names and 'toFrame'
# 
#     xdat <- toFrame(xdat)
#     
#     ## maintain y names and 'toFrame'
#     ydat <- toFrame(ydat)
#     
#     ## first grab dummy args for bandwidth() and perform 'bootstrap'
#     ## bandwidth() call
#     
#     mc.names <- names(match.call(expand.dots = FALSE))
#     margs <- c("bwmethod", "bwscaling", "bwtype", "cxkertype", "cxkerorder",
#                "cykertype", "cykerorder", "uxkertype", "uykertype", "oxkertype",
#                "oykertype")
#     
#     m <- match(margs, mc.names, nomatch = 0)
#     any.m <- any(m != 0)
#     browser()
#     tbw <- eval(parse(text=paste("conbandwidth(",
#                                  "xbw = bws[length(ydat)+1:length(xdat)],",
#                                  "ybw = bws[1:length(ydat)],",
#                                  paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
#                                  ifelse(any.m, ",",""),
#                                  "nobs = nrow(xdat),",
#                                  "xdati = untangle(xdat),",
#                                  "ydati = untangle(ydat),",
#                                  "xnames = names(xdat),",
#                                  "ynames = names(ydat),",
#                                  "bandwidth.compute = bandwidth.compute)")))
#     
#     ## next grab dummies for actual bandwidth selection and perform call
#     
#     mc.names <- names(match.call(expand.dots = FALSE))
#     margs <- c("bandwidth.compute", "auto", "nmulti", "remin", "itmax", "ftol",
#                "tol", "small")
#     m <- match(margs, mc.names, nomatch = 0)
#     any.m <- any(m != 0)
#     
#     tbw <- eval(parse(text=paste("npcdensbw.conbandwidth(xdat=xdat, ydat=ydat, bws=tbw",
#                                  ifelse(any.m, ",",""),
#                                  paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
#                                  ")")))
#     
#     mc <- match.call(expand.dots = FALSE)
#     environment(mc) <- parent.frame()
#     tbw$call <- mc
#     
#     return(tbw)
#   }
# ## No Zero Denominator, used in C code for kernel estimation...
# 
# NZD <- function(a) {
#   sapply(1:NROW(a), function(i) {if(a[i] < 0) min(-.Machine$double.xmin,a[i]) else max(.Machine$double.xmin,a[i])})
# }
# 
# ## Function to test for monotone increasing vector
# 
# is.monotone.increasing <- function(x) {
#   ## Sorted and last value > first value
#   !is.unsorted(x) && x[length(x)] > x[1]
# }
# 
# ## what is a badord? ... an ordered factor of numeric values to treat
# ## them properly one must preserve the numeric value, ie. scale not
# ## just their sorted order Actually, the ord/badord paradigm must go,
# ## in place of levels caching
# 
# matrix.sd <- function(x, na.rm=FALSE) {
#   if(is.matrix(x)) apply(x, 2, sd, na.rm=na.rm)
#   else if(is.vector(x)) sd(x, na.rm=na.rm)
#   else if(is.data.frame(x)) sapply(x, sd, na.rm=na.rm)
#   else sd(as.vector(x), na.rm=na.rm)
# }
# 
# npseed <- function(seed){
#   .C("np_set_seed",as.integer(abs(seed)), PACKAGE = "np")
#   invisible()
# }
# 
# numNotIn <- function(x){
#   while(is.element(num <- rnorm(1),x)){}
#   num
# }
# 
# dlev <- function(x){
#   if(is.ordered(x))
#     x.dlev <- suppressWarnings(as.numeric(levels(x)))
#   if (!is.ordered(x) || any(is.na(x.dlev)))
#     x.dlev <- as.numeric(1:nlevels(x))
#   x.dlev
# }
# 
# isNum <- function(x){
#   return(!any(is.na(suppressWarnings(as.numeric(x)))))
# }
# 
# untangle <- function(frame){
#   if (is.null(frame))
#     return(NULL)
#   
#   iord <- unlist(lapply(frame, is.ordered))
#   iuno <- unlist(lapply(frame, is.factor)) & !iord
#   icon <- unlist(lapply(frame, is.numeric))
#   
#   if(!all(iord | iuno | icon)) 
#     stop("non-allowed data type in data frame")
#   
#   inumord <-
#     suppressWarnings(unlist(lapply(frame,
#                                    function (z) {
#                                      is.ordered(z) && is.numeric(tryCatch(as.numeric(levels(z)), warning = function (y) {
#                                        FALSE }))
#                                    })))
#   
#   all.lev <- lapply(frame, function(y){
#     t.ret <- NULL
#     if (is.factor(y))
#       t.ret <- levels(y)
#     t.ret
#   })
#   
#   all.ulev <- lapply(frame, function (y) {
#     t.ret <- NULL
#     if (is.factor(y))
#       t.ret <- sort(unique(y))
#     t.ret
#   })
#   
#   all.dlev <- lapply(frame, function (y) {
#     t.ret <- NULL
#     if (is.factor(y))
#       t.ret <- dlev(y)
#     t.ret
#   })
#   
#   all.nlev <- lapply(frame, function(y) {
#     t.ret <- NULL
#     if (is.factor(y))
#       t.ret <- nlevels(y)
#     t.ret
#   })
#   
#   list(iord = iord,
#        iuno = iuno,
#        icon = icon,
#        inumord = inumord,
#        all.lev = all.lev,
#        all.ulev = all.ulev,
#        all.dlev = all.dlev,
#        all.nlev = all.nlev)
# }
# 
# validateBandwidth <- function(bws){
#   vari <- names(bws$bandwidth)
#   bchecker <- function(j){
#     v <- vari[j]
#     dati <- bws$dati[[v]]
#     bwv <- bws$bandwidth[[j]]
#     stopifnot(length(bwv) == length(dati$iord))
#     
#     cd <- function(a,b){
#       (a-b)/(a+b+.Machine$double.eps) > 5.0*.Machine$double.eps
#     }
#     
#     vb <- sapply(1:length(bwv), function(i){
#       falg <- (bwv[i] < 0)
#       
#       if (dati$icon[i] && (falg || (!is.finite(bwv[i])))){
#         stop(paste("Invalid bandwidth supplied for continuous",
#                    "variable", bws$varnames[[v]][i], ":",bwv[i]))
#       }
#       
#       if (dati$iord[i] &&
#             (falg || cd(bwv[i],oMaxL(dati$all.nlev[[i]],
#                                      kertype = bws$klist[[v]]$okertype)))){
#         stop(paste("Invalid bandwidth supplied for ordered",
#                    "variable", bws$varnames[[v]][i], ":",bwv[i]))
#       }
#       
#       if (dati$iuno[i] &&
#             (falg || cd(bwv[i],uMaxL(dati$all.nlev[[i]],
#                                      kertype = bws$klist[[v]]$ukertype)))){
#         stop(paste("Invalid bandwidth supplied for unordered",
#                    "variable", bws$varnames[[v]][i], ":",bwv[i]))
#       }
#       return(TRUE)
#     })
#     
#     return(vb)
#   }
#   vbl <- lapply(1:length(vari), bchecker)
#   invisible(vbl)
# }
# 
# explodeFormula <- function(formula){
#   res <- strsplit(strsplit(paste(deparse(formula), collapse=""),
#                            " *[~] *")[[1]], " *[+] *")
#   stopifnot(all(sapply(res,length) > 0))
#   names(res) <- c("response","terms")
#   res
# }
# 
# 
# explodePipe <- function(formula){
#   tf <- as.character(formula)  
#   tf <- tf[length(tf)]
#   
#   eval(parse(text=paste("c(",
#                         ifelse(length(as.character(formula)) == 3,
#                                'strsplit(as.character(formula)[2]," *[+] *"),',""),
#                         'strsplit(strsplit(tf," *[|] *")[[1]]," *[+] *"))')))
# }
# 
# "%~%" <- function(a,b) {
#   all(class(a) == class(b)) && (length(a) == length(b)) &&
#     all(unlist(lapply(a,coarseclass)) == unlist(lapply(b,coarseclass)))
# }
# 
# coarseclass <- function(a) {
#   ifelse(class(a) == "integer", "numeric", class(a))
# }
# 
# toFrame <- function(frame) {
#   if(!is.data.frame(frame)){
#     t.names <- NULL
#     
#     if(!(is.vector(frame) || is.factor(frame) || is.matrix(frame)))
#       stop(deparse(substitute(frame))," must be a data frame, matrix, vector, or factor")
#     
#     if(!is.matrix(frame))
#       t.names <- deparse(eval(substitute(substitute(frame)), envir = parent.frame()))
#     
#     frame <- data.frame(frame, check.names=FALSE)
#     
#     if(!is.null(t.names))
#       names(frame) <- t.names
#   }
#   return(frame)
# }
# 
# 
# cast <- function(a, b, same.levels = TRUE){
#   if(is.ordered(b)){
#     if(same.levels)
#       ordered(a, levels = levels(b))
#     else
#       ordered(a)
#   }   
#   else if(is.factor(b)){
#     if(same.levels)
#       factor(a, levels = levels(b))
#     else
#       factor(a)
#   }
#   else if (coarseclass(b) == "numeric")
#     as.double(a)
#   else if (is.data.frame(b)) {
#     if (dim(a)[2] == dim(b)[2]){
#       r = data.frame(a)
#       for (i in 1: length(b))
#         r[,i] = cast(a[,i],b[,i], same.levels = same.levels)
#       r
#     } else { stop("a could not be cast as b") }
#   }
# }
# 
# subcol <- function(x, v, i){
#   x[,i] = cast(v,x[,i])
#   x
# }
# 
# mcvConstruct <- function(dati){
#   nuno <- sum(dati$iuno)
#   nord <- sum(dati$iord)
#   
#   num.row <- max(sapply(dati$all.lev,length))
#   
#   pad.num <- numNotIn(unlist(dati$all.dlev))
#   
#   mcv <- matrix(data = pad.num, nrow = num.row, ncol = (nuno+nord))
#   attr(mcv, "num.row") <- num.row
#   attr(mcv, "pad.num") <- pad.num
#   
#   cnt <- 0
#   if (nuno > 0)
#     for (i in which(dati$iuno)) 
#       mcv[1:length(dati$all.lev[[i]]), (cnt <- cnt+1)] <- dati$all.dlev[[i]]
#   
#   cnt <- 0
#   if (nord > 0)
#     for (i in which(dati$iord))
#       mcv[1:length(dati$all.lev[[i]]), (cnt <- cnt+1)+nuno] <- dati$all.dlev[[i]]
#   
#   mcv
# }
# 
# ## when admitting new categories, adjustLevels will attempt to catch possible mistakes:
# ## if an unordered variable contains more than one new category, warn
# ## if an ordered, but scaleless variable contains a new category, error
# ## if an ordered, scale-possessing variable contains a new category lacking scale, error
# 
# adjustLevels <- function(data, dati, allowNewCells = FALSE){
#   for (i in which(dati$iord | dati$iuno)){
#     if (allowNewCells){
#       newCats <- setdiff(levels(data[,i]),dati$all.lev[[i]])
#       if (length(newCats) >= 1){
#         if (dati$iuno[i]){
#           if (length(newCats) > 1)
#             warning(paste("more than one 'new' category is redundant when estimating on unordered data.\n",
#                           "training data categories: ", paste(dati$all.lev[[i]], collapse=" "),"\n",
#                           "redundant estimation data categories: ", paste(newCats, collapse=" "), "\n", sep=""))
#           data[,i] <- factor(data[,i], levels = c(dati$all.lev[[i]], newCats))
#         } else {
#           if (dati$inumord[i]){
#             if (!isNum(newCats))
#               stop(paste("estimation data contains a new qualitative category, but training data is\n",
#                          "ordered, and numeric.\n",
#                          "training data categories: ", paste(dati$all.lev[[i]], collapse=" "),"\n",
#                          "conflicting estimation data categories: ", paste(newCats, collapse=" "), "\n", sep=""))
#           } else {
#             stop(paste("estimation beyond the support of training data of an ordered,\n",
#                        "categorical, qualitative variable is not supported.\n"))
#           }
#           
#           data[,i] <- ordered(data[,i], levels = sort(as.numeric(c(dati$all.lev[[i]], newCats))))
#         }
#       } else {
#         data[,i] <- factor(data[,i], levels = dati$all.lev[[i]])
#       }
#     } else {
#       if (!all(is.element(levels(data[,i]), dati$all.lev[[i]])))
#         stop("data contains unknown factors (wrong dataset provided?)")
#       data[,i] <- factor(data[,i], levels = dati$all.lev[[i]])
#     }
#   }
#   
#   data
# }
# 
# toMatrix <- function(data) {
#   tq <- sapply(data, function(y) {
#     if(is.factor(y))
#       y <- dlev(y)[as.integer(y)]
#     y})
#   dim(tq) <- dim(data) ## cover the corner case of single element d.f.
#   tq
# }
# 
# ## this doesn't just strictly check for the response, but does tell you
# ## that evaluating with response fails ... in principle the evaluation
# ## could fail without the response too, but then the calling routine is about
# ## to die a noisy death anyhow ...
# succeedWithResponse <- function(tt, frame){
#   !any(class(try(eval(expr = attr(tt, "variables"),
#                       envir = frame, enclos = NULL), silent = TRUE)) == "try-error")
# }
# 
# ## determine whether a bandwidth
# ## matches a data set
# bwMatch <- function(data, dati){
#   if (length(dati$icon) != ncol(data))
#     stop("bandwidth vector is of improper length")
#   
#   test.dati <- untangle(data)
#   
#   if (any(xor(dati$icon,test.dati$icon)) ||
#         any(xor(dati$iord,test.dati$iord)) ||
#         any(xor(dati$iuno,test.dati$iuno)))
#     stop(paste("supplied bandwidths do not match","data", "in type"))
# }
# 
# uMaxL <- function(c, kertype = c("aitchisonaitken","liracine")){
#   switch(kertype,
#          aitchisonaitken = (c-1.0)/c,
#          liracine = 1.0)
# }
# 
# oMaxL <- function(c, kertype = c("wangvanryzin", "liracine")){
#   switch(kertype,
#          wangvanryzin = 1.0,
#          liracine = 1.0)
# }
# 
# ## tested with: rbandwidth
# ## right now all bandwidth objects have some crusty
# ## vestiges of their evolution, ie. non-list metadata
# ## such as xnames or ynames. The new metadata system is
# ## for the most part list based and facilitates generic
# ## operations
# 
# updateBwNameMetadata <- function(nameList, bws){
#   ## names of 'old' metadata in bw object
#   onames <- names(nameList)
#   lapply(1:length(nameList), function(i) {
#     bws[[onames[i]]] <<- nameList[[i]]
#     bws$varnames[[substr(onames[i],1,1)]] <<- nameList[[i]]
#   })
#   return(bws)
# }
# 
# ## some string utility functions
# 
# pad <- function(s){
#   ifelse(nchar(s) > 0, paste("",s,""), " ")
# }
# 
# rpad <- function(s){
#   ifelse(nchar(s) > 0, paste(s,""), "")
# }
# 
# blank <- function(len){
#   sapply(len, function(nb){
#     paste(rep(' ', times = nb), collapse='')
#   })
# }
# 
# formatv <- function(v){
#   sapply(1:length(v), function(j){ format(v[j]) })
# }
# 
# ## strings used in various report generating functions
# 
# genOmitStr <- function(x){
#   t.str <- ''
#   if(!is.null(x$rows.omit) & !identical(x$rows.omit, NA))
#     t.str <- paste("\nNo. Complete Observations: ", x$nobs,
#                    "No. NA Observations: ", x$nobs.omit,
#                    "\nObservations omitted: ", paste(x$rows.omit, collapse=" "),
#                    "\n")
#   return(t.str)
# }
# 
# ## Estimation-related rgf's
# genGofStr <- function(x){
#   paste(ifelse(is.na(x$MSE),"",paste("\nResidual standard error:",
#                                      format(sqrt(x$MSE)))),
#         ifelse(is.na(x$R2),"",paste("\nR-squared:",
#                                     format(x$R2))), sep="")
# }
# 
# pCatGofStr <- function(x){
#   if(!identical(x$confusion.matrix, NA)){
#     cat("\nConfusion Matrix\n")
#     print(x$confusion.matrix)
#   }
#   
#   if (!identical(x$CCR.overall,NA))
#     cat("\nOverall Correct Classification Ratio: ", format(x$CCR.overall))
#   
#   if (!identical(x$CCR.byoutcome,NA)){
#     cat("\nCorrect Classification Ratio By Outcome:\n")
#     print(x$CCR.byoutcome)
#   }
#   
#   if (!identical(x$fit.mcfadden,NA))
#     cat("\nMcFadden-Puig-Kerschner performance measure: ", format(x$fit.mcfadden))
#   
# }
# 
# genDenEstStr <- function(x){
#   paste("\nBandwidth Type: ",x$ptype,
#         ifelse(is.null(x$log_likelihood) || identical(x$log_likelihood, NA),"",
#                paste("\nLog Likelihood:",
#                      format(x$log_likelihood))),
#         sep="")
# }
# 
# genRegEstStr <- function(x){
#   paste(ifelse(is.null(x$pregtype),"",
#                paste("\nKernel Regression Estimator:",x$pregtype)),
#         ifelse(is.null(x$ptype), "",
#                paste("\nBandwidth Type:",x$ptype)),
#         ifelse(is.null(x$tau), "", paste("\nTau:", x$tau)),
#         sep = "")
# }
# 
# 
# ## bandwidth-related report generating functions
# genBwSelStr <- function(x){
#   paste(ifelse(is.null(x$pregtype),"",paste("\nRegression Type:", x$pregtype)),
#         ifelse(is.null(x$pmethod),"",paste("\nBandwidth Selection Method:",
#                                            x$pmethod)),
#         if (!identical(x$formula,NULL)) paste("\nFormula:",
#                                               paste(deparse(x$formula), collapse="\n")),
#         ifelse(is.null(x$ptype), "",
#                paste("\nBandwidth Type: ",x$ptype, sep="")),
#         ifelse(identical(x$fval,NA),"",
#                paste("\nObjective Function Value: ", format(x$fval),
#                      " (achieved on multistart ", x$ifval, ")", sep="")), sep="")
# }
# 
# genBwScaleStrs <- function(x){
#   ## approach is to take metadata and flatten it so it then can be
#   ## processed into a single string
#   
#   vari <- names(x$sumNum)
#   
#   t.icon <- lapply(vari, function(v){
#     x$dati[[v]]$icon })
#   
#   sumText <- lapply(1:length(vari), function(i) {
#     ifelse(t.icon[[i]],
#            ifelse(x$type == "fixed",
#                   "Scale Factor:",""), "Lambda Max:")
#     
#   })
#   
#   maxNameLen <- max(nchar(unlist(sumText)))
#   print.sumText <- lapply(sumText, '!=', "")
#   
#   sumText <- lapply(1:length(sumText), function(i){
#     paste(blank(maxNameLen - nchar(sumText[[i]])), sumText[[i]], sep="")
#   })
#   
#   t.nchar <- lapply(x$varnames[vari], nchar)
#   
#   maxNameLen <- max(unlist(t.nchar))
#   
#   vatText <- lapply(1:length(t.nchar), function(j){
#     paste("\n", rpad(x$vartitleabb[[vari[j]]]), "Var. Name: ",
#           x$varnames[[vari[j]]],
#           sapply(t.nchar[[j]], function(nc){
#             paste(rep(' ', maxNameLen - nc), collapse='')
#           }), sep='')
#   })
#   
#   return(sapply(1:length(t.nchar), function(j){
#     paste(vatText[[j]], " Bandwidth: ", npFormat(x$bandwidth[[j]]), " ",
#           ifelse(print.sumText[[j]],
#                  paste(sumText[[j]], " ", npFormat(x$sumNum[[j]]), sep=""), ""),
#           sep="", collapse="")
#   }))
# }
# 
# npFormat <- function(x){
#   format(sapply(x,format))
# }
# 
# genBwKerStrs <- function(x){
#   vari <- names(x$klist)
#   
#   ncon <- sapply(vari, function(v){
#     sum(x$dati[[v]]$icon)
#   })
#   
#   nuno <- sapply(vari, function(v){
#     sum(x$dati[[v]]$iuno)
#   })
#   
#   nord <- sapply(vari, function(v){
#     sum(x$dati[[v]]$iord)
#   })
#   
#   cktype <- sapply(vari, function(v){
#     x$klist[[v]]$ckertype
#   })
#   
#   uktype <- sapply(vari, function(v){
#     x$klist[[v]]$ukertype
#   })
#   
#   oktype <- sapply(vari, function(v){
#     x$klist[[v]]$okertype
#   })
#   
#   tt <- ''
#   
#   if(any(ncon > 0)){
#     tt <- paste("\n",
#                 ifelse(length(unique(cktype)) == 1,
#                        paste("\nContinuous Kernel Type:",
#                              x$klist[[vari[1]]]$pckertype),
#                        paste(sapply(1:length(vari), function(v){
#                          ifelse(ncon[v] > 0,
#                                 paste("\nContinuous Kernel Type (",
#                                       x$vartitleabb[[vari[v]]],
#                                       " Var.): ", x$klist[[vari[v]]]$pckertype, sep=""),"")
#                        }), collapse = "")),
#                 sep = "")
#     tt <-
#       paste(tt, paste(sapply(1:length(vari), function(i){
#         ifelse(ncon[i] > 0,
#                paste("\nNo. Continuous", pad(x$vartitle[[vari[i]]]), "Vars.: ",
#                      ncon[i], sep=""), "")
#       }), collapse = ""), sep="")
#   }
#   
#   
#   if(any(nuno > 0)) {
#     tt <- paste(tt, "\n",
#                 ifelse(length(unique(uktype)) == 1,
#                        paste("\nUnordered Categorical Kernel Type:",
#                              x$klist[[vari[1]]]$pukertype),
#                        paste(sapply(1:length(vari), function(i){
#                          ifelse(nuno[i] > 0,
#                                 paste("\nUnordered Categorical Kernel Type (",
#                                       x$vartitleabb[[vari[i]]],
#                                       " Var.): ", x$klist[[vari[i]]]$pukertype, sep=""),"")
#                        }), collapse = "")),
#                 sep = "")
#     tt <-
#       paste(tt, paste(sapply(1:length(vari), function(i){
#         ifelse(nuno[i] > 0,
#                paste("\nNo. Unordered Categorical", pad(x$vartitle[[vari[i]]]), "Vars.: ",
#                      nuno[i], sep=""), "")
#       }), collapse = ""), sep="")
#     
#   }
#   
#   if(any(nord > 0)) {
#     tt <- paste(tt, "\n",
#                 ifelse(length(unique(oktype)) == 1,
#                        paste("\nOrdered Categorical Kernel Type:",
#                              x$klist[[vari[1]]]$pokertype),
#                        paste(sapply(1:length(vari), function(i){
#                          ifelse(nord[i] > 0,
#                                 paste("\nOrdered Categorical Kernel Type (",
#                                       x$vartitleabb[[vari[i]]],
#                                       " Var.): ", x$klist[[vari[i]]]$pokertype, sep=""),"")
#                        }), collapse = "")),
#                 sep = "")
#     tt <-
#       paste(tt, paste(sapply(1:length(vari), function(i){
#         ifelse(nord[i] > 0,
#                paste("\nNo. Ordered Categorical", pad(x$vartitle[[vari[i]]]), "Vars.: ",
#                      nord[i], sep=""), "")
#       }), collapse = ""), sep="")
#     
#   }
#   
#   return(tt)
# }
# 
# genBwKerStrsXY <- function(x){
#   t.str <- ''
#   cnt <- 0
#   
#   if (x$xncon + x$yncon > 0){
#     t.str[cnt <- cnt + 1] <-
#       paste(ifelse(x$pcxkertype == x$pcykertype,
#                    paste("\n\nContinuous Kernel Type:",x$pcxkertype),
#                    paste("\n", ifelse(x$xncon > 0,
#                                       paste("\nContinuous Kernel Type (Exp. Var.):",
#                                             x$pcxkertype), ""),
#                          ifelse(x$yncon > 0,
#                                 paste("\nContinuous Kernel Type (Dep. Var.):",
#                                       x$pcykertype), ""))),
#             ifelse(x$yncon > 0,paste("\nNo. Continuous Dependent Vars.:",x$yncon),""),
#             ifelse(x$xncon > 0,paste("\nNo. Continuous Explanatory Vars.:",x$xncon),""))
#   }
#   
#   if (x$xnuno + x$ynuno > 0){
#     t.str[cnt <- cnt + 1] <-
#       paste(ifelse(x$puxkertype == x$puykertype,
#                    paste("\n\nUnordered Categorical Kernel Type:",x$puxkertype),
#                    paste("\n", ifelse(x$xnuno > 0,
#                                       paste("\nUnordered Categorical Kernel Type (Exp. Var.):",
#                                             x$puxkertype), ""),
#                          ifelse(x$ynuno > 0,
#                                 paste("\nUnordered Categorical Kernel Type (Dep. Var.):",
#                                       x$puykertype), ""))),
#             ifelse(x$ynuno > 0,paste("\nNo. Unordered Categorical Dependent Vars.:",x$ynuno),""),
#             ifelse(x$xnuno > 0,paste("\nNo. Unordered Categorical Explanatory Vars.:",x$xnuno),""))
#   }
#   
#   if (x$xnord + x$ynord > 0){
#     t.str[cnt <- cnt + 1] <-
#       paste(ifelse(x$poxkertype == x$poykertype,
#                    paste("\n\nOrdered Categorical Kernel Type:",x$poxkertype),
#                    paste("\n", ifelse(x$xnord > 0,
#                                       paste("\nOrdered Categorical Kernel Type (Exp. Var.):",
#                                             x$poxkertype), ""),
#                          ifelse(x$ynord > 0,
#                                 paste("\nOrdered Categorical Kernel Type (Dep. Var.):",
#                                       x$poykertype), ""))),
#             ifelse(x$ynord > 0,paste("\nNo. Ordered Categorical Dependent Vars.:",x$ynord),""),
#             ifelse(x$xnord > 0,paste("\nNo. Ordered Categorical Explanatory Vars.:",x$xnord),""))
#   }
#   return(t.str)
# }
# 
# genBwGOFStrs <- function(x) {
#   ###paste("Residual standard error:", sqrt
# }
# 
# ## statistical functions
# RSQfunc <- function(y,y.pred) {
#   y.mean <- mean(y)
#   (sum((y-y.mean)*(y.pred-y.mean))^2)/(sum((y-y.mean)^2)*sum((y.pred-y.mean)^2))
# }
# 
# MSEfunc <- function(y,y.fit) {
#   mean((y-y.fit)^2)
# }
# 
# MAEfunc <- function(y,y.fit) {
#   mean(abs(y-y.fit))
# }
# 
# MAPEfunc <- function(y,y.fit) {
#   jj = which(y != 0)
#   
#   mean(c(abs((y[jj]-y.fit[jj])/y[jj]), as.numeric(replicate(length(y)-length(jj),2))))
# }
# 
# CORRfunc <- function(y,y.fit) {
#   abs(corr(cbind(y,y.fit)))
# }
# 
# SIGNfunc <- function(y,y.fit) {
#   sum(sign(y) == sign(y.fit))/length(y)
# }
# 
# 
# EssDee <- function(y){
#   
#   sd.vec <- apply(as.matrix(y),2,sd)
#   IQR.vec <- apply(as.matrix(y),2,IQR)/(qnorm(.25,lower.tail=F)*2)
#   return(ifelse(sd.vec<IQR.vec|IQR.vec==0,sd.vec,IQR.vec))
#   
# }
# 
# ### holding place for some generic methods
# 
# se <- function(x){
#   UseMethod("se",x)
# }
# 
# gradients <- function(x, ...){
#   UseMethod("gradients",x)
# }
# ### internal constants used in the c backend
# 
# SF_NORMAL = 0
# SF_ARB = 1
# 
# BW_FIXED = 0
# BW_GEN_NN = 1
# BW_ADAP_NN = 2
# 
# IMULTI_TRUE = 1
# IMULTI_FALSE = 0
# 
# RE_MIN_TRUE = 0
# RE_MIN_FALSE = 1
# 
# IO_MIN_TRUE = 1
# IO_MIN_FALSE = 0
# 
# USE_START_NO = 0
# USE_START_YES = 1
# 
# NP_DO_DENS = 1
# NP_DO_DIST = 0
# 
# ## initially making an np-wide option via the 'options' mechanism
# DO_TREE_NO = 0
# DO_TREE_YES = 1
# 
# ##kernel defs
# CKER_GAUSS = 0
# CKER_EPAN  = 4
# CKER_UNI   = 8
# 
# UKER_AIT = 0
# UKER_LR = 1
# 
# OKER_WANG = 0
# OKER_LR = 1
# 
# ##density 
# BWM_CVML = 0
# BWM_CVLS = 1
# BWM_CVML_NP= 2
# 
# ##distribution
# DBWM_CVLS = 0
# 
# ##regression
# BWM_CVAIC = 0
# 
# REGTYPE_LC = 0
# REGTYPE_LL = 1
# 
# ##conditional density/distribution
# CBWM_CVML = 0
# CBWM_CVLS = 1
# CBWM_NPLS = 2
# CBWM_CCDF = 3 # Added 7/2/2010 jracine
# 
# ##conditional distribution
# CDBWM_CVLS = 0
# 
# ##integral operators on kernels
# OP_NORMAL      = 0
# OP_CONVOLUTION = 1
# OP_DERIVATIVE  = 2
# OP_INTEGRAL    = 3
# 
# ALL_OPERATORS = c(OP_NORMAL, OP_CONVOLUTION, OP_DERIVATIVE, OP_INTEGRAL)
# names(ALL_OPERATORS) <- c("normal","convolution", "derivative", "integral")
# 
# ## useful numerical constants of kernel integrals
# int.kernels <- c(0.28209479177387814348, 0.47603496111841936711, 0.62396943688265038571, 0.74785078617543927990,
#                  0.26832815729997476357, 0.55901699437494742410, 0.84658823667359826246, 1.1329342579014329689,
#                  0.5)
# 
# npcdensbw <-
#   function(...){
#     args = list(...)
#     if (is(args[[1]],"formula"))
#       UseMethod("npcdensbw",args[[1]])
#     else if (!is.null(args$formula))
#       UseMethod("npcdensbw",args$formula)
#     else
#       UseMethod("npcdensbw",args[[which(names(args)=="bws")[1]]])
#   }
# 
# npcdensbw.formula <-
#   function(formula, data, subset, na.action, call, ...){
# 
#     mf <- match.call(expand.dots = FALSE)
#     m <- match(c("formula", "data", "subset", "na.action"),
#                names(mf), nomatch = 0)
#     mf <- mf[c(1,m)]
#     
#     if(!missing(call) && is.call(call)){
#       ## rummage about in the call for the original formula
#       for(i in 1:length(call)){
#         if(tryCatch(class(eval(call[[i]])) == "formula",
#                     error = function(e) FALSE))
#           break;
#       }
#       mf[[2]] <- call[[i]]
#       
#     }
#     
#     
#     mf[[1]] <- as.name("model.frame")
#     
#     variableNames <- explodeFormula(mf[["formula"]])
#     
#     ## make formula evaluable, then eval
#     varsPlus <- lapply(variableNames, paste, collapse=" + ")
#     mf[["formula"]] <- as.formula(paste(" ~ ", varsPlus[[1]]," + ",
#                                         varsPlus[[2]]),
#                                   env = environment(formula))
#     
#     mf <- eval(mf, parent.frame())
#     
#     ydat <- mf[, variableNames[[1]], drop = FALSE]
#     xdat <- mf[, variableNames[[2]], drop = FALSE]
#     
#     tbw = npcdensbw(xdat = xdat, ydat = ydat, ...)
#     
#     ## clean up (possible) inconsistencies due to recursion ...
#     tbw$call <- match.call(expand.dots = FALSE)
#     environment(tbw$call) <- parent.frame()
#     tbw$formula <- formula
#     tbw$rows.omit <- as.vector(attr(mf,"na.action"))
#     tbw$nobs.omit <- length(tbw$rows.omit)
#     tbw$terms <- attr(mf,"terms")
#     tbw$variableNames <- variableNames
#     
#     tbw
#   }
# 
# npcdensbw.conbandwidth <- 
#   function(xdat = stop("data 'xdat' missing"),
#            ydat = stop("data 'ydat' missing"),
#            bws, bandwidth.compute = TRUE,
#            auto = TRUE,
#            nmulti, remin = TRUE, itmax = 10000,
#            ftol=1.19209e-07, tol=1.49012e-08, small=2.22045e-16,
#            ...){
#     
#     ydat = toFrame(ydat)
#     xdat = toFrame(xdat)
#     browser()
#     
#     if (missing(nmulti)){
#       nmulti <- min(5,(dim(ydat)[2]+dim(xdat)[2]))
#     }
#     
#     if (length(bws$ybw) != dim(ydat)[2])
#       stop(paste("length of bandwidth vector does not match number of columns of", "'ydat'"))
#     
#     if (length(bws$xbw) != dim(xdat)[2])
#       stop(paste("length of bandwidth vector does not match number of columns of", "'xdat'"))
#     
#     if (dim(ydat)[1] != dim(xdat)[1])
#       stop(paste("number of rows of", "'ydat'", "does not match", "'xdat'"))
#     
#     yccon = unlist(lapply(as.data.frame(ydat[,bws$iycon]),class))
#     if ((any(bws$iycon) && !all((yccon == class(integer(0))) | (yccon == class(numeric(0))))) ||
#           (any(bws$iyord) && !all(unlist(lapply(as.data.frame(ydat[,bws$iyord]),class)) ==
#                                     class(ordered(0)))) ||
#           (any(bws$iyuno) && !all(unlist(lapply(as.data.frame(ydat[,bws$iyuno]),class)) ==
#                                     class(factor(0)))))
#       stop(paste("supplied bandwidths do not match", "'ydat'", "in type"))
#     
#     xccon = unlist(lapply(as.data.frame(xdat[,bws$ixcon]),class))
#     if ((any(bws$ixcon) && !all((xccon == class(integer(0))) | (xccon == class(numeric(0))))) ||
#           (any(bws$ixord) && !all(unlist(lapply(as.data.frame(xdat[,bws$ixord]),class)) ==
#                                     class(ordered(0)))) ||
#           (any(bws$ixuno) && !all(unlist(lapply(as.data.frame(xdat[,bws$ixuno]),class)) ==
#                                     class(factor(0)))))
#       stop(paste("supplied bandwidths do not match", "'xdat'", "in type"))
#     
#     ## catch and destroy NA's
#     goodrows <- 1:dim(xdat)[1]
#     rows.omit <- unclass(na.action(na.omit(data.frame(xdat,ydat))))
#     goodrows[rows.omit] <- 0
#     
#     if (all(goodrows==0))
#       stop("Data has no rows without NAs")
#     
#     xdat = xdat[goodrows,,drop = FALSE]
#     ydat = ydat[goodrows,,drop = FALSE]
#     
#     
#     nrow = nrow(ydat)
#     yncol = ncol(ydat)
#     xncol = ncol(xdat)
#     
#     ## at this stage, data to be sent to the c routines must be converted to
#     ## numeric type.
#     
#     ydat = toMatrix(ydat)
#     
#     yuno = ydat[, bws$iyuno, drop = FALSE]
#     ycon = ydat[, bws$iycon, drop = FALSE]
#     yord = ydat[, bws$iyord, drop = FALSE]
#     
#     
#     xdat = toMatrix(xdat)
#     
#     xuno = xdat[, bws$ixuno, drop = FALSE]
#     xcon = xdat[, bws$ixcon, drop = FALSE]
#     xord = xdat[, bws$ixord, drop = FALSE]
#     
#     tbw <- bws
#     
#     if (bandwidth.compute){
#       myopti = list(num_obs_train = nrow,
#                     iMultistart = ifelse(nmulti==0,IMULTI_FALSE,IMULTI_TRUE),
#                     iNum_Multistart = nmulti,
#                     int_use_starting_values = ifelse(all(bws$ybw==0) && all(bws$xbw==0),
#                                                      USE_START_NO, USE_START_YES),
#                     int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
#                     BANDWIDTH_den_extern = switch(bws$type,
#                                                   fixed = BW_FIXED,
#                                                   generalized_nn = BW_GEN_NN,
#                                                   adaptive_nn = BW_ADAP_NN),
#                     itmax=itmax, int_RESTART_FROM_MIN=ifelse(remin,RE_MIN_TRUE,RE_MIN_FALSE), 
#                     int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE), 
#                     bwmethod = switch(bws$method,
#                                       cv.ml = CBWM_CVML,
#                                       cv.ls = CBWM_CVLS,
#                                       cv.ls.np = CBWM_NPLS),        
#                     xkerneval = switch(bws$cxkertype,
#                                        gaussian = CKER_GAUSS + bws$cxkerorder/2 - 1,
#                                        epanechnikov = CKER_EPAN + bws$cxkerorder/2 - 1,
#                                        uniform = CKER_UNI),
#                     ykerneval = switch(bws$cykertype,
#                                        gaussian = CKER_GAUSS + bws$cykerorder/2 - 1,
#                                        epanechnikov = CKER_EPAN + bws$cykerorder/2 - 1,
#                                        uniform = CKER_UNI),
#                     uxkerneval = switch(bws$uxkertype,
#                                         aitchisonaitken = UKER_AIT,
#                                         liracine = UKER_LR),
#                     uykerneval = switch(bws$uykertype,
#                                         aitchisonaitken = UKER_AIT,
#                                         liracine = UKER_LR),
#                     oxkerneval = switch(bws$oxkertype,
#                                         wangvanryzin = OKER_WANG,
#                                         liracine = OKER_LR),
#                     oykerneval = switch(bws$oykertype,
#                                         wangvanryzin = OKER_WANG,
#                                         liracine = OKER_LR),
#                     ynuno = dim(yuno)[2],
#                     ynord = dim(yord)[2],
#                     yncon = dim(ycon)[2],
#                     xnuno = dim(xuno)[2],
#                     xnord = dim(xord)[2],
#                     xncon = dim(xcon)[2],
#                     fast = FALSE,
#                     auto = auto)
#       
#       myoptd = list(ftol=ftol, tol=tol, small=small)
#       
#       if (bws$method != "normal-reference"){
#         myout=
#           .C("np_density_conditional_bw", as.double(yuno), as.double(yord), as.double(ycon),
#              as.double(xuno), as.double(xord), as.double(xcon),
#              as.integer(myopti), as.double(myoptd), 
#              bw = c(bws$xbw[bws$ixcon],bws$ybw[bws$iycon],
#                     bws$ybw[bws$iyuno],bws$ybw[bws$iyord],
#                     bws$xbw[bws$ixuno],bws$xbw[bws$ixord]),
#              fval = double(2),
#              PACKAGE="np" )[c("bw","fval")]
#       } else {
#         nbw = double(yncol+xncol)
#         gbw = bws$yncon+bws$xncon
#         if (gbw > 0){
#           nbw[1:gbw] = (4/3)^0.2
#           if(!bws$scaling)
#             nbw[1:gbw]=nbw[1:gbw]*EssDee(data.frame(xcon,ycon))*nrow^(-1.0/(2.0*bws$cxkerorder+gbw))
#         }
#         myout= list( bw = nbw, fval = c(NA,NA) )
#       }
#       
#       yr = 1:yncol
#       xr = 1:xncol
#       rorder = numeric(yncol + xncol)
#       
#       ## bandwidths are passed back from the C routine in an unusual order
#       ## xc, y[cuo], x[uo]
#       
#       rxcon = xr[bws$ixcon]
#       rxuno = xr[bws$ixuno] 
#       rxord = xr[bws$ixord] 
#       
#       rycon = yr[bws$iycon] 
#       ryuno = yr[bws$iyuno] 
#       ryord = yr[bws$iyord] 
#       
#       
#       ## rorder[c(rxcon,rycon,ryuno,ryord,rxuno,rxord)]=1:(yncol+xncol)
#       
#       tbw <- bws
#       tbw$ybw[c(rycon,ryuno,ryord)] <- myout$bw[yr+bws$xncon]
#       tbw$xbw[c(rxcon,rxuno,rxord)] <- myout$bw[setdiff(1:(yncol+xncol),yr+bws$xncon)]
#       
#       tbw$fval = myout$fval[1]
#       tbw$ifval = myout$fval[2]
#     }
#     
#     ## bandwidth metadata
#     tbw$sfactor <- tbw$bandwidth <- list(x = tbw$xbw, y = tbw$ybw)
#     
#     bwf <- function(i){
#       tbw$bandwidth[[i]][tl[[i]]] <<- (tbw$bandwidth[[i]])[tl[[i]]]*dfactor[[i]]
#     }
#     
#     sff <- function(i){
#       tbw$sfactor[[i]][tl[[i]]] <<- (tbw$sfactor[[i]])[tl[[i]]]/dfactor[[i]]
#     }
#     
#     myf <- if(tbw$scaling) bwf else sff
#     
#     if ((tbw$xnuno+tbw$ynuno) > 0){
#       dfactor <- nrow^(-2.0/(2.0*tbw$cxkerorder+tbw$ncon))
#       dfactor <- list(x = dfactor, y = dfactor)
#       
#       tl <- list(x = tbw$xdati$iuno, y = tbw$ydati$iuno)
#       
#       lapply(1:length(tl), myf)
#     }
#     
#     if ((tbw$xnord+tbw$ynord) > 0){
#       dfactor <- nrow^(-2.0/(2.0*tbw$cxkerorder+tbw$ncon))
#       dfactor <- list(x = dfactor, y = dfactor)
#       
#       tl <- list(x = tbw$xdati$iord, y = tbw$ydati$iord)
#       
#       lapply(1:length(tl), myf)
#     }
#     
#     
#     if (tbw$ncon > 0){
#       dfactor <- nrow^(-1.0/(2.0*tbw$cxkerorder+tbw$ncon))
#       dfactor <- list(x = EssDee(xcon)*dfactor, y = EssDee(ycon)*dfactor)
#       
#       tl <- list(x = tbw$xdati$icon, y = tbw$ydati$icon)
#       
#       lapply(1:length(tl), myf)
#     }
#     
#     tbw <- conbandwidth(xbw = tbw$xbw,
#                         ybw = tbw$ybw,
#                         bwmethod = tbw$method,
#                         bwscaling = tbw$scaling,
#                         bwtype = tbw$type,
#                         cxkertype = tbw$cxkertype,
#                         cxkerorder = tbw$cxkerorder,
#                         uxkertype = tbw$uxkertype,
#                         oxkertype = tbw$oxkertype,
#                         cykertype = tbw$cykertype,
#                         cykerorder = tbw$cykerorder,
#                         uykertype = tbw$uykertype,
#                         oykertype = tbw$oykertype,
#                         fval = tbw$fval,
#                         ifval = tbw$ifval,
#                         nobs = tbw$nobs,
#                         xdati = tbw$xdati,
#                         ydati = tbw$ydati,      
#                         xnames = tbw$xnames,
#                         ynames = tbw$ynames,
#                         sfactor = tbw$sfactor,
#                         bandwidth = tbw$bandwidth,
#                         rows.omit = rows.omit,
#                         bandwidth.compute = bandwidth.compute)
#     browser()
#     
#     tbw
#   }
# 
# npcdensbw.NULL <-
#   function(xdat = stop("data 'xdat' missing"),
#            ydat = stop("data 'ydat' missing"),
#            bws, ...){
#     
#     ## maintain x names and 'toFrame'
#     browser()
#     xdat <- toFrame(xdat)
#     
#     ## maintain y names and 'toFrame'
#     ydat <- toFrame(ydat)
#     
#     ## do bandwidths
#     
#     bws = double(ncol(ydat)+ncol(xdat))
#     
#     tbw <- npcdensbw.default(xdat = xdat, ydat = ydat, bws = bws, ...)
#     
#     ## clean up (possible) inconsistencies due to recursion ...
#     mc <- match.call(expand.dots = FALSE)
#     environment(mc) <- parent.frame()
#     tbw$call <- mc
#     
#     tbw
#   }
# 
# npcdensbw.default <-
#   function(xdat = stop("data 'xdat' missing"),
#            ydat = stop("data 'ydat' missing"),
#            bws, 
#            bandwidth.compute = TRUE,
#            auto, nmulti, remin, itmax,
#            ftol, tol, small,
#            ## dummy arguments for conbandwidth() function call
#            bwmethod, bwscaling, bwtype,
#            cxkertype, cxkerorder,
#            cykertype, cykerorder,
#            uxkertype, uykertype,
#            oxkertype, oykertype,
#            ...){
#     
#     ## maintain x names and 'toFrame'
#     xdat <- toFrame(xdat)
#     
#     ## maintain y names and 'toFrame'
#     ydat <- toFrame(ydat)
#     
#     ## first grab dummy args for bandwidth() and perform 'bootstrap'
#     ## bandwidth() call
#     
#     mc.names <- names(match.call(expand.dots = FALSE))
#     margs <- c("bwmethod", "bwscaling", "bwtype", "cxkertype", "cxkerorder",
#                "cykertype", "cykerorder", "uxkertype", "uykertype", "oxkertype",
#                "oykertype")
#     
#     m <- match(margs, mc.names, nomatch = 0)
#     any.m <- any(m != 0)
#     
#     tbw <- eval(parse(text=paste("conbandwidth(",
#                                  "xbw = bws[length(ydat)+1:length(xdat)],",
#                                  "ybw = bws[1:length(ydat)],",
#                                  paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
#                                  ifelse(any.m, ",",""),
#                                  "nobs = nrow(xdat),",
#                                  "xdati = untangle(xdat),",
#                                  "ydati = untangle(ydat),",
#                                  "xnames = names(xdat),",
#                                  "ynames = names(ydat),",
#                                  "bandwidth.compute = bandwidth.compute)")))
#     
#     ## next grab dummies for actual bandwidth selection and perform call
#     
#     mc.names <- names(match.call(expand.dots = FALSE))
#     margs <- c("bandwidth.compute", "auto", "nmulti", "remin", "itmax", "ftol",
#                "tol", "small")
#     m <- match(margs, mc.names, nomatch = 0)
#     any.m <- any(m != 0)
#     
#     tbw <- eval(parse(text=paste("npcdensbw.conbandwidth(xdat=xdat, ydat=ydat, bws=tbw",
#                                  ifelse(any.m, ",",""),
#                                  paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
#                                  ")")))
#     
#     mc <- match.call(expand.dots = FALSE)
#     environment(mc) <- parent.frame()
#     tbw$call <- mc
#     
#     return(tbw)
#   }
# conbandwidth <-
#   function(xbw,
#            ybw,
#            bwmethod = c("cv.ml","cv.ls","normal-reference", "cv.ls.np", "manual"),
#            bwscaling = FALSE,
#            bwtype = c("fixed","generalized_nn","adaptive_nn"),
#            cxkertype = c("gaussian", "epanechnikov","uniform"), 
#            cxkerorder = c(2,4,6,8),
#            uxkertype = c("aitchisonaitken","liracine"),
#            oxkertype = c("wangvanryzin","liracine"),
#            cykertype = c("gaussian", "epanechnikov","uniform"), 
#            cykerorder = c(2,4,6,8),
#            uykertype = c("aitchisonaitken"),
#            oykertype = c("wangvanryzin"),
#            fval = NA,
#            ifval = NA,
#            nobs = NA,
#            xdati, ydati,
#            xnames = character(length(xbw)),
#            ynames = character(length(ybw)),
#            sfactor = NA, bandwidth = NA,
#            rows.omit = NA, bandwidth.compute = TRUE,...){
#     
#     if (missing(xbw) | missing(ybw))
#       stop("improper invocation of conbandwidth constructor: 'bw' or i[cuo]* missing")
#     
#     xndim = length(xbw)
#     yndim = length(ybw)
#     browser()
#     
#     bwmethod = match.arg(bwmethod)
#     bwtype = match.arg(bwtype)
#     
#     cxkertype = match.arg(cxkertype)
#     cykertype = match.arg(cykertype)
#     
#     if(missing(cxkerorder))
#       cxkerorder = 2
#     else if (cxkertype == "uniform")
#       warning("ignoring kernel order specified with uniform kernel type")
#     else {
#       kord = eval(formals()$cxkerorder) 
#       if (!any(kord == cxkerorder))
#         stop("cxkerorder must be one of ", paste(kord,collapse=" "))
#     }
#     
#     if (bwmethod == "normal-reference" && (cxkertype != "gaussian" || bwtype != "fixed")){    
#       warning("normal-reference bandwidth selection assumes gaussian kernel with fixed bandwidth")
#       bwtype = "fixed"
#       cxkertype = "gaussian"
#     }
#     
#     if(missing(cykerorder))
#       cykerorder = 2
#     else if (cykertype == "uniform")
#       warning("ignoring kernel order specified with uniform kernel type")
#     else {
#       kord = eval(formals()$cykerorder) 
#       if (!any(kord == cykerorder))
#         stop("cykerorder must be one of ", paste(kord,collapse=" "))
#     }
#     
#     if (bwmethod == "normal-reference" && (cykertype != "gaussian" || bwtype != "fixed")){    
#       warning("normal-reference bandwidth selection assumes gaussian kernel with fixed bandwidth")
#       bwtype = "fixed"
#       cykertype = "gaussian"
#     }
#     
#     if (cxkerorder != cykerorder & bwscaling)
#       stop("scale factors with different order kernels for dependent and explanatory variables is unsupported")
#     
#     uxkertype = match.arg(uxkertype)
#     uykertype = match.arg(uykertype)
#     
#     oxkertype = match.arg(oxkertype)
#     oykertype = match.arg(oykertype)
#     
#     pxorder = switch( cxkerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )
#     pyorder = switch( cykerorder/2, "Second-Order", "Fourth-Order", "Sixth-Order", "Eighth-Order" )
#     
#     dati <- list(x = xdati, y = ydati)
#     
#     if (!identical(sfactor,NA)){
#       ## using the new model for generically accessing bandwidth objects
#       
#       okertype <- list(x = oxkertype, y = oykertype)
#       ukertype <- list(x = uxkertype, y = uykertype)
#       
#       scaleOrMax <- function(i, j) {
#         if (dati[[j]]$icon[i])
#           return((sfactor[[j]])[i])
#         
#         if (dati[[j]]$iord[i])
#           return(oMaxL(dati[[j]]$all.nlev[[i]], kertype = okertype[[j]]))
#         
#         if (dati[[j]]$iuno[i])
#           return(uMaxL(dati[[j]]$all.nlev[[i]], kertype = ukertype[[j]]))
#       }
#       
#       sumNum <- list(x = NA, y = NA)
#       sumNum[] <- lapply(1:length(dati), function(i) {
#         sapply(1:length(dati[[i]]$icon), scaleOrMax, j = i)
#       })
#     } else {
#       sumNum <- NA
#     }
#     
#     if (length(rows.omit) == 0)
#       rows.omit <- NA
#     
#     mybw = list(
#       xbw=xbw,
#       ybw=ybw,
#       method = bwmethod,
#       pmethod = switch( bwmethod,
#                         cv.ml = "Maximum Likelihood Cross-Validation",
#                         cv.ls = "Least Squares Cross-Validation",
#                         cv.ls.np = "Least Squares Cross-Validation (block algorithm)",
#                         "normal-reference" = "Normal Reference"),
#       fval = fval,
#       ifval = ifval,
#       scaling = bwscaling,
#       pscaling = ifelse(bwscaling, "Scale Factor(s)", "Bandwidth(s)"),
#       type = bwtype,
#       ptype = switch( bwtype,
#                       fixed = "Fixed",
#                       generalized_nn = "Generalized Nearest Neighbour",
#                       adaptive_nn = "Adaptive Nearest Neighbour" ),
#       cxkertype = cxkertype,
#       cykertype = cykertype,
#       cxkerorder = cxkerorder,
#       cykerorder = cykerorder,
#       pcxkertype = switch(cxkertype,
#                           gaussian = paste(pxorder,"Gaussian"),
#                           epanechnikov =  paste(pxorder,"Epanechnikov"),
#                           uniform = "Uniform"),
#       pcykertype = switch(cykertype,
#                           gaussian = paste(pyorder,"Gaussian"),
#                           epanechnikov =  paste(pyorder,"Epanechnikov"),
#                           uniform = "Uniform"),
#       uxkertype = uxkertype,
#       uykertype = uykertype,
#       puxkertype = switch( uxkertype,
#                            aitchisonaitken = "Aitchison and Aitken",
#                            liracine = "Li and Racine"),
#       puykertype = switch( uykertype,
#                            aitchisonaitken = "Aitchison and Aitken"),
#       oxkertype = oxkertype,
#       oykertype = oykertype,
#       poxkertype = switch( oxkertype,
#                            wangvanryzin = "Wang and Van Ryzin",
#                            liracine = "Li and Racine"),
#       poykertype = switch( oykertype,
#                            wangvanryzin = "Wang and Van Ryzin"),
#       nobs = nobs,
#       xndim = xndim,
#       yndim = yndim,
#       ndim = xndim + yndim,
#       xncon = sum(xdati$icon),
#       xnuno = sum(xdati$iuno),
#       xnord = sum(xdati$iord),
#       yncon = sum(ydati$icon),
#       ynuno = sum(ydati$iuno),
#       ynord = sum(ydati$iord),
#       ncon = sum(c(xdati$icon, ydati$icon)),
#       ixcon = xdati$icon,
#       ixuno = xdati$iuno,
#       ixord = xdati$iord,
#       iycon = ydati$icon,
#       iyuno = ydati$iuno,
#       iyord = ydati$iord,
#       xnames = xnames,
#       ynames = ynames,
#       xdati = xdati,
#       ydati = ydati,
#       xmcv = mcvConstruct(xdati),
#       ymcv = mcvConstruct(ydati),
#       sfactor = sfactor,
#       bandwidth = bandwidth,
#       sumNum = sumNum,
#       dati = dati, 
#       varnames = list(x = xnames, y = ynames),
#       vartitle = list(x = "Explanatory", y = "Dependent"),
#       vartitleabb = list(x = "Exp.", y = "Dep."),
#       rows.omit = rows.omit,
#       nobs.omit = ifelse(identical(rows.omit,NA), 0, length(rows.omit)))
#     
#     mybw$klist = list(
#       x =
#         list(ckertype = cxkertype,
#              pckertype = mybw$pcxkertype,
#              ukertype = uxkertype,
#              pukertype = mybw$puxkertype,
#              okertype = oxkertype,
#              pokertype = mybw$poxkertype),
#       y =
#         list(ckertype = cykertype,
#              pckertype = mybw$pcykertype,
#              ukertype = uykertype,
#              pukertype = mybw$puykertype,
#              okertype = oykertype,
#              pokertype = mybw$poykertype))
#     
#     if(!bandwidth.compute)
#       mybw$pmethod <- "Manual"
#     
#     
#     class(mybw) = "conbandwidth"
#     browser()
#     if(!any(is.na(mybw$bandwidth)))
#       validateBandwidth(mybw)
#     browser()
#     mybw
#   }
# 
# print.conbandwidth <- function(x, digits=NULL, ...){
#   cat("\nConditional density data (",x$nobs," observations, ",
#       (x$xndim+x$yndim)," variable(s))",
#       "\n(", x$yndim, " dependent variable(s), and ", x$xndim, " explanatory variable(s))\n\n",
#       sep="")
#   print(matrix(x$ybw,ncol=x$yndim,dimnames=list(paste("Dep. Var. ",x$pscaling,":",sep=""),x$ynames)))
#   
#   print(matrix(x$xbw,ncol=x$xndim,dimnames=list(paste("Exp. Var. ",x$pscaling,":",sep=""),x$xnames)))
#   
#   cat(genBwSelStr(x))
#   cat(genBwKerStrsXY(x))
#   
#   cat("\n\n")
#   if(!missing(...))
#     print(...,digits=digits)
#   invisible(x)
# }
# 
# plot.conbandwidth <- function(...) { npplot(...) }
# 
# summary.conbandwidth <- function(object, ...) {
#   cat("\nConditional density data (",object$nobs," observations, ",
#       (object$xndim+object$yndim)," variable(s))",
#       "\n(", object$yndim, " dependent variable(s), and ", object$xndim, " explanatory variable(s))\n",
#       sep="")
#   
#   cat(genOmitStr(object))
#   cat(genBwSelStr(object))
#   
#   cat(paste("\n", genBwScaleStrs(object), sep=""))
#   cat(genBwKerStrs(object))
#   
#   cat("\n\n")
# }
# 
# predict.conbandwidth <- function(...) { eval(npcdens(...), envir = parent.frame()) }
# 
