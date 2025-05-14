install_packages <- function() {
  is_windows <- .Platform$OS.type == "windows"
  
  # Define target_lib only on Windows
  if (is_windows) {
    get_libpath <- function() {
      libpaths <- .libPaths()
      target_lib <- libpaths[grepl("one.*extensions$", libpaths)]
      
      if (length(target_lib) == 0) {
        stop("No valid lib path found with 'one/extensions'")
      }
      
      return(target_lib)
    }
    
    libpath <- get_libpath()
    
    load_or_install <- function(pkg, libpath) {
      is_installed <- suppressWarnings(requireNamespace(pkg, quietly = TRUE, lib.loc = libpath))
      if (!is_installed) {
        cat(paste0("Installing '", pkg, "' to ", libpath, "\n"))
        install.packages(pkg, lib = libpath, dependencies = TRUE, type = "binary")
        if (pkg == "Matrix") assign("was_installed", TRUE, envir = .GlobalEnv)
      }
    }
    
    load_or_install("Matrix", libpath)
    load_or_install("np", libpath)
    
  } else {
    # For macOS or other platforms
    load_or_install <- function(pkg) {
      is_installed <- suppressWarnings(requireNamespace(pkg, quietly = TRUE))
      if (!is_installed) {
        cat(paste0("Installing '", pkg, "' to default library\n"))
        install.packages(pkg, dependencies = TRUE)  # default lib, default type
        if (pkg == "Matrix") assign("was_installed", TRUE, envir = .GlobalEnv)
      }
    }
    
    load_or_install("Matrix")
    load_or_install("np")
  }

  # Show restart note if Matrix was just installed
  if (exists("was_installed", envir = .GlobalEnv)) {
    cat("\nNOTE: Please restart SPSS to reload packages.\n")
    cat("      This is necessary because 'Matrix' and 'np' was just installed or updated.\n")
  }
}

# Run the installation
install_packages()


# Main nonparametric regression function with debug statements
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
	
	# Debug: Log function entry and input parameters
	#cat("DEBUG: Entering npext function\n")
	#cat("DEBUG: dep =", dep, "\n")
	#cat("DEBUG: indepnp =", paste(indepnp, collapse=", "), "\n")
	#cat("DEBUG: id =", id, "\n")
	
	nindepnp = length(indepnp)
	nstart = ifelse(is.null(nstart), min(5, nindepnp), nstart)
	frml = buildfrml(dep, indepnp)
	
	# Debug: Log formula and nindepnp
	#cat("DEBUG: Formula created:", frml, "\n")
	#cat("DEBUG: nindepnp =", nindepnp, "\n")
	
	allargs = as.list(environment())  # for passing external spec to other functions

  # Debug: Before calling setuplocalization
  #cat("DEBUG: Calling setuplocalization with domain STATS_NONPAR_REGR\n")
  tryCatch({
    setuplocalization("STATS_NONPAR_REGR")
    #cat("setuplocalization completed successfully\n")
  }, error = function(e) {
    cat("Error in setuplocalizationJune 06, 2025 ERROR: setuplocalization failed:", e$message, "\n")
    stop(e$message)
  })

  # A warnings proc name is associated with the regular output
  # (and the same omsid), because warnings/errors may appear in
  # a separate procedure block following the regular output
  procname=gtxt("Nonparametric Regression")
  warningsprocname = gtxt("Nonparametric Regression: Warnings")
  omsid="STATSNPREG"
  warns = Warn(procname=warningsprocname,omsid=omsid)

  # Debug: Before loading np package
  #cat("DEBUG: Attempting to load np package\n")
  #cat("DEBUG: Current .libPaths():", paste(.libPaths(), collapse=", "), "\n")
  tryCatch(library(np, quietly=TRUE), error=function(e){
      #cat("Failed to load np package:", e$message, "\n")
      warns$warn(gtxtf("The R %s package is required but could not be loaded.","np"),
          dostop=TRUE)
      }
  )
  #cat("DEBUG: np package loaded successfully\n")
  
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
	
	# Debug: Before fetching data from SPSS
	#cat("DEBUG: Fetching data from SPSS for variables:", paste(allvars, collapse=", "), "\n")
	dta = tryCatch(spssdata.GetDataFromSPSS(allvars, row.label=id, missingValueToNA=TRUE,
		factorMode="labels"), error=function(e) {
      cat("Error fetching SPSS data:", e$message, "\n")
      warns$warn(e, dostop=TRUE)
		}
  )
	# Debug: Log data frame summary
	#cat("DEBUG: Data fetched, dimensions:", dim(dta), "\n")
	
	# remove any factor levels that have a zero count
	# NB: this can't reattach labels to string factors
	dta = cleaned(dta)
	# Debug: Log cleaned data frame
	#cat("DEBUG: Data cleaned, dimensions:", dim(dta), "\n")
	
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
		# Debug: Before computing bandwidth
		#cat("DEBUG: Computing bandwidth for continuous dependent variable\n")
		bw = tryCatch({
			npregbw(formula=as.formula(frml),data=dta, na.action=naaction,
				regtype = regtype, bwmethod=bwmethod, bwscaling=bwscaling, 
				bwtype=bwtype, ckertype=cvkerneltype, ckerorder=cvorder,
				ukertype=ukerneltype, okertype=okerneltype, nmulti=nstart,
				ftol=ftol, tol=tol)
			},
			error = function(e) {
				cat("Error in npregbw:", e$message, "\n")
				warns$warn(e$message, dostop=TRUE)
			}
		)
		# Debug: Log bandwidth object
		#cat("DEBUG: Bandwidth computed, class:", class(bw), "\n")
			
		res = tryCatch(
      npreg(bws=bw, data=dta, gradients=(plots == "gradient")),
			error=function(e) {
        cat("Error in npreg:", e$message, "\n")
        warns$warn(e$message, dostop=TRUE)
			}
		)
	} else {
		# Debug: Before computing bandwidth for categorical variable
		#cat("DEBUG: Computing bandwidth for categorical dependent variable\n")
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
				cat("Error in npcdensbw:", e$message, "\n")
				warns$warn(e$message, dostop=TRUE)
			}
		)
		# Debug: Log bandwidth object
		#cat("DEBUG: Bandwidth computed, class:", class(bw), "\n")
			
		res = tryCatch(
      npconmode(bws=bw, data=dta, gradients=(plots == "gradient")),
			error=function(e) {
        cat("Error in npconmode:", e$message, "\n")
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
			# Debug: Before significance tests
			#cat("DEBUG: Running significance tests, sigindex =", sigindex, "\n")
			sigres = tryCatch({
				npsigtest(bw, xdat=dta[,2:(nindepnp+1)], ydat=dta[,1], data=dta,
				index=1:sigindex, na.action=naaction,
				boot.method=bootmethod, boot.num=bootreps, random.seed=rnseed)
			},
			error = function(e) {
				cat("Error in npsigtest:", e$message, "\n")
				warns$warn(e$message, dostop=TRUE)
			}
			)
		}
	}
	
  # print results
  # browser()  # Commented out to avoid pausing execution
  displayresults(allargs, allvars, dta, bw, res, sigres, warns)
  if (!is.null(dataset)) {
	  createdataset(res, dta, allargs, warns)
  }

  # clean up workspace
  res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
  # Debug: Log function exit
  #cat("DEBUG: Exiting npext function\n")
}

# Localization initialization with debug statements
setuplocalization = function(domain) {
  # find and bind translation file names
  # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
  # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

  # Debug: Log entry and domain
  #cat("DEBUG: Entering setuplocalization, domain =", domain, "\n")
  #cat("DEBUG: Current .libPaths():", paste(.libPaths(), collapse=", "), "\n")
  
  # Debug: Log file path construction
  target_file = paste(domain, ".R", sep="")
  #cat("DEBUG: Searching for file:", target_file, "\n")
  
  fpath = Find(file.exists, file.path(.libPaths(), target_file))
  # Debug: Log file path result
  #if (is.null(fpath)) {
    #cat("File not found:", target_file, "\n")
  #} else {
    #cat("File found at:", fpath, "\n")
  #}
  
  # Debug: Log localization directory
  loc_dir = file.path(dirname(fpath), domain, "lang")
  #cat("DEBUG: Localization directory:", loc_dir, "\n")
  #cat("DEBUG: Localization directory exists:", file.exists(loc_dir), "\n")
  
  # Debug: Before binding text domain
  #cat("DEBUG: Binding text domain:", domain, "to", loc_dir, "\n")
  result = bindtextdomain(domain, loc_dir)
  # Debug: Log result of bindtextdomain
  #cat("DEBUG: bindtextdomain result:", result, "\n")
}

# The rest of the script remains unchanged
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
		info["cvtype"],
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
		gtxt("*).Objective Function"),
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
		row.names(condf) = dimnames(res$collision_matrix)$Actual
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
      # browser()  # Commented out
  		tryCatch(plot(bw, gradients=(plots=="gradient")),
  			error=function(e) {
  				warns$warn(gtxt("Gradient plot could not be produced.  There may be too many variables"),		
            dostop=FALSE)
  			}
  		)
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