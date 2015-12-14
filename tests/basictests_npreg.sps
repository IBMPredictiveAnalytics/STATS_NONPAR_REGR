﻿get file="c:/spss22/samples/english/employee data.sav".
dataset name emp.
compute salbinary = salary gt 40000.
compute salbinaryrev = 1-salbinary.
compute id = -id.
variable level salbinary salbinaryrev(ordinal).
exec.
STATS NONPAR REGR DEPENDENT=salary   indep=jobtime prevexp minority
/bandwidth nstart=1.

STATS NONPAR REGR DEPENDENT=salary   indep=jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN NSTART=1
/OUTPUT SIGTESTS=YES NVARS=999
 BOOTMETHOD=IID BOOTREPS=9  RNSEED=42 PLOTS=GRADIENT.

STATS NONPAR REGR DEPENDENT=salary   indep=jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN NSTART=1
/OUTPUT SIGTESTS=YES NVARS=999
 BOOTMETHOD=IID BOOTREPS=9  RNSEED=42 PLOTS=RESPONSE ERRORBNDS =asymptotic.

STATS NONPAR REGR DEPENDENT=salary   indep=jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN NSTART=1
/OUTPUT SIGTESTS=YES NVARS=999
 BOOTMETHOD=IID BOOTREPS=9  RNSEED=42 PLOTS=RESPONSE.

STATS NONPAR REGR DEPENDENT=salbinary   indep=jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN NSTART=1
/OUTPUT SIGTESTS=YES NVARS=999
 BOOTMETHOD=IID BOOTREPS=9  RNSEED=42 PLOTS=GRADIENT.

STATS NONPAR REGR DEPENDENT=salary   indep=jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN NSTART=1
/OUTPUT SIGTESTS=YES NVARS=999
 BOOTMETHOD=IID BOOTREPS=9  RNSEED=42 PLOTS=RESPONSE errorbnds=asymptotic.

STATS NONPAR REGR DEPENDENT=salary   indep=jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN NSTART=1
/OUTPUT SIGTESTS=YES NVARS=999
 BOOTMETHOD=IID BOOTREPS=9  RNSEED=42 PLOTS=RESPONSE errorbnds=asymptotic.

STATS NONPAR REGR DEPENDENT=salary   indep=jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN NSTART=1
/OUTPUT SIGTESTS=YES NVARS=999
 BOOTMETHOD=IID BOOTREPS=9  RNSEED=42 PLOTS=RESPONSE includepts=no errorbnds=asymptotic.


variable level educ(scale).
BEGIN PROGRAM R.
dta <- spssdata.GetDataFromSPSS(factorMode="labels")
res = summary(lm(salary ~ educ + jobcat + prevexp + jobtime, data = dta))
print(res)
spsspivottable.Display(coef(res), templateName="SPSSLM", title="My Regression",
caption = sprintf("R-Squared: %.3f Residual SE: %.3f", res$r.squared, res$sigma) )
end program.


STATS NONPAR REGR DEPENDENT=salary   indep=jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=GENERALIZEDNN SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN  nstart=1
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=NO.


STATS NONPAR REGR DEPENDENT=salary   indep=jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=GENERALIZEDNN SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN nstart=1
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=NO PLOTS=RESPONSE.

DATASET ACTIVATE emp.
STATS NONPAR REGR DEPENDENT=salbinary   indep=jobcat jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
NSTART=1
/OUTPUT SIGTESTS=NO PLOTS=RESPONSE.

STATS NONPAR REGR DEPENDENT=salbinary   indep=jobcat jobtime prevexp minority 
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
NSTART=1
/OUTPUT SIGTESTS=NO  PLOTS=RESPONSE.

STATS NONPAR REGR DEPENDENT=salbinary   indep=jobcat
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
NSTART=1
/OUTPUT SIGTESTS=NO  PLOTS=RESPONSE.

STATS NONPAR REGR DEPENDENT=salbinaryrev   indep=jobcat
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
NSTART=1
/OUTPUT SIGTESTS=NO  PLOTS=RESPONSE.

DATASET ACTIVATE emp.
STATS NONPAR REGR DEPENDENT=salary   INDEP=educ jobtime 
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=YES NVARS=999 BOOTMETHOD=WILD BOOTREPS=399 RNSEED=42 PLOTS=RESPONSE.

STATS NONPAR REGR DEPENDENT=salary   INDEP=educ jobtime
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVLS CVTYPE=GENERALIZEDNN SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=YES NVARS=999 BOOTMETHOD=WILD BOOTREPS=399 RNSEED=42 PLOTS=GRADIENT.
*fixed vs generalizednn vs adaptive.
STATS NONPAR REGR DEPENDENT=salary   INDEP= jobtime
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVLS CVTYPE=fixed SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=YES NVARS=999 BOOTMETHOD=WILD BOOTREPS=399 RNSEED=42 PLOTS=RESPONSE.
STATS NONPAR REGR DEPENDENT=salary   INDEP= jobtime
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVLS CVTYPE=GENERALIZEDNN SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=YES NVARS=999 BOOTMETHOD=WILD BOOTREPS=399 RNSEED=42 PLOTS=RESPONSE.

* dataset creation.
DATASET ACTIVATE emp.
STATS NONPAR REGR DEPENDENT=salary   INDEP=educ jobtime ID=id
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=NO PLOTS=NONE ERRORBNDS=NONE
/SAVE DATASET=results PREDICTED=YES RESIDUALS=YES.

STATS NONPAR REGR DEPENDENT=salary   INDEP=educ jobtime
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS CVTYPE=FIXED SCALING=RAW CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=NO PLOTS=NONE ERRORBNDS=NONE
/SAVE DATASET=resultsresonly RESIDUALS=YES.

* parameter variations.
DATASET ACTIVATE emp.
STATS NONPAR REGR DEPENDENT=salary   INDEP=jobcat jobtime ID=id
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVAIC CVTYPE=ADAPTIVENN SCALING=SCALE CVKERNELTYPE=UNIFORM     
CVORDER=4 OKERNELTYPE=LIRACINE UKERNELTYPE=LIRACINE NSTART=1
/OUTPUT SIGTESTS=YES NVARS=1 BOOTMETHOD=IID BOOTREPS=100 RNSEED=42 PLOTS=RESPONSE 
    ERRORBNDS=ASYMPTOTIC
/SAVE  RESIDUALS=YES PREDICTED=YES.

STATS NONPAR REGR DEPENDENT=salary   INDEP=jobcat jobtime ID=id
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVAIC CVTYPE=ADAPTIVENN SCALING=SCALE CVKERNELTYPE=UNIFORM     
CVORDER=4 OKERNELTYPE=LIRACINE UKERNELTYPE=LIRACINE NSTART=1
/OUTPUT SIGTESTS=NO PLOTS=RESPONSE ERRORBNDS=NONE
/SAVE RESIDUALS=YES PREDICTED=YES.

STATS NONPAR REGR DEPENDENT=salary   INDEP=jobcat jobtime ID=id
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVAIC CVTYPE=ADAPTIVENN SCALING=SCALE CVKERNELTYPE=UNIFORM     
CVORDER=4 OKERNELTYPE=LIRACINE UKERNELTYPE=LIRACINE NSTART=1
/OUTPUT SIGTESTS=NO PLOTS=RESPONSE ERRORBNDS=NONE
/SAVE RESIDUALS=YES PREDICTED=YES.

DATASET ACTIVATE emp.
STATS NONPAR REGR DEPENDENT=salbinary   INDEP=educ jobtime ID=id
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVAIC ESTIMATORCAT=NORMALREF  CVTYPE=ADAPTIVENN 
    SCALING=SCALE CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=NO PLOTS=RESPONSE ERRORBNDS=NONE.


get file="%er%\stats_nonpar_regr\tests\cps71.sav".
dataset name cps.
DATASET ACTIVATE cps.
STATS NONPAR REGR DEPENDENT=logwage   INDEP=age ID=id
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS ESTIMATORCAT=CVML  CVTYPE=FIXED SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=YES NVARS=999 BOOTMETHOD=WILD BOOTREPS=399 RNSEED=42 PLOTS=RESPONSE ERRORBNDS=NONE.    
   

STATS NONPAR REGR DEPENDENT=logwage   INDEP=age ID=id
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVLS ESTIMATORCAT=CVML  CVTYPE=FIXED SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=YES NVARS=999 BOOTMETHOD=WILD BOOTREPS=399 RNSEED=42 PLOTS=RESPONSE 
    ERRORBNDS=ASYMPTOTIC.

STATS NONPAR REGR DEPENDENT=logwage   INDEP=age ID=id
/BANDWIDTH REGTYPE=LOCALLINEAR ESTIMATORCONT=CVLS ESTIMATORCAT=CVML  CVTYPE=FIXED SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OUTPUT SIGTESTS=YES NVARS=999 BOOTMETHOD=WILD BOOTREPS=399 RNSEED=42 PLOTS=RESPONSE 
    ERRORBNDS=BOOTSTRAP.

DATASET ACTIVATE cps.
STATS NONPAR REGR DEPENDENT=logwage   INDEP=age ID=id
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS ESTIMATORCAT=CVML  CVTYPE=FIXED SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OPTIONS MISSING=OMIT MAXITER=10000 VALTOLERANCE=.0000000119 LOCTOLERANCE=.00000000149
/OUTPUT SIGTESTS=NO PLOTS=RESPONSE ERRORBNDS=NONE
INCLUDEPTS=NO.

get file="%er%\stats_nonpar_regr\tests\italy.sav".
dataset name italy.

DATASET ACTIVATE italy.
STATS NONPAR REGR DEPENDENT=gdp   INDEP=year ID=ID
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS ESTIMATORCAT=CVML  CVTYPE=FIXED SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=4 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OPTIONS MISSING=STOP MAXITER=10000 VALTOLERANCE=.0000000119 LOCTOLERANCE=.00000000149
/OUTPUT SIGTESTS=NO PLOTS=RESPONSE ERRORBNDS=ASYMPTOTIC
INCLUDEPTS=YES.
variable level year(ordinal).

STATS NONPAR REGR DEPENDENT=gdp   INDEP=year ID=ID
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS ESTIMATORCAT=CVML  CVTYPE=FIXED SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=4 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OPTIONS MISSING=STOP MAXITER=10000 VALTOLERANCE=.0000000119 LOCTOLERANCE=.00000000149
/OUTPUT SIGTESTS=NO PLOTS=GRADIENT ERRORBNDS=ASYMPTOTIC
INCLUDEPTS=YES.
variable level year(nominal).




DATASET ACTIVATE emp.
STATS NONPAR REGR DEPENDENT=jobcat   INDEP=educ gender jobtime ID=id
/BANDWIDTH REGTYPE=LOCALCONSTANT ESTIMATORCONT=CVLS ESTIMATORCAT=CVML  CVTYPE=FIXED SCALING=RAW 
    CVKERNELTYPE=GAUSSIAN 
CVORDER=2 OKERNELTYPE=WANGVANRYZIN UKERNELTYPE=AITCHISONAITKEN 
/OPTIONS MISSING=OMIT MAXITER=10000 VALTOLERANCE=.0000000119 LOCTOLERANCE=.00000000149
/OUTPUT SIGTESTS=NO PLOTS=RESPONSE ERRORBNDS=ASYMPTOTIC
INCLUDEPTS=NO /SAVE DATASET=results RESIDUALS=NO PREDICTED=YES.