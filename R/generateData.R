#' Generate simulated data
#'
#' \code{generateData} generates a set of simulated clinical trials' data
#'
#' This is the core data simulation code. See vignettes for additional details.
#' This function is primarily called by \link{generateSimulatedResults}.
#'
#' @param modelparam a datatable with a single entry per named column:
#' \itemize{
#'  \item{\code{N}}{  Number of simulated participants}
#'  \item{\code{c.bm}}{  Correlation between the biomarker and the biologic response}
#'  \item{\code{carryover_t1half}}{  Halflife of the carryover effect}
#'  \item{\code{c.tv}}{  Autocorrelation for the tv factor across timepoints}
#'  \item{\code{c.pb}}{ Autocorrelation for the pb factor across timepoints}
#'  \item{\code{c.br}}{  Autocorrelation for the br factor across timepoints}
#'  \item{\code{c.cf1t}}{  Correlation between different factors at a single timepoint}
#'  \item{\code{c.cfct}}{  Correlation between different factors at different timepoints}
#' }
#' @param respparam Data table with 5 columns and 3 rows. The first row lists the 3 factors
#'   defining the response trajectories (tv, pb, and br). The rest of the columns are the
#'   input parameters for the \link{modgompertz} function:
#' \itemize{
#'  \item{\code{max}}{  Maximum response value}
#'  \item{\code{disp}}{  Displacement}
#'  \item{\code{rate}}{  Rate}
#'  \item{\code{sd}}{  Standard deviation}
#' }
#' @param blparam Data table with 3 colums and 3 rows. First row lists the varaibles being
#'   defined, which consist of \code{BL} for the baseline value of the outcome measure, and
#'   \code{bm} for the biomarker values. Remaining columns:
#'   \itemize{
#'     \item{\code{m}}{  Mean}
#'     \item{\code{sd}}{  Standard deviation}
#'   }
#' @param trialdesign The information defining a single path through the clinical trial
#'   design, as defined by the \code{output$trialpaths} output of \link{buildtrialdesign}
#' @param empirical Should the correlations be empirical? (Usually \code{FALSE} unless
#'   you're doing something odd for testing purposes)
#' @param makePositiveDefinite Should the covariance matrix be forced to be positive definite?
#'   (Usually yes, although you may want to turn this off at times to test the impact of
#'   this step on your covariance sturcture)
#' @param seed Randomization seed (defaults to NA)
#' @param lambda_cor Rate of BM-BR correlation decay during off-drug
#'   periods (per week). Controls how fast the biomarker's predictive
#'   power decays after drug discontinuation. The correlation decays
#'   as exp(-lambda_cor * tsd). Defaults to NA, which computes the
#'   pharmacokinetically consistent value ln(2) / carryover_t1half,
#'   matching the drug elimination rate. Set to 0 to disable
#'   correlation decay (original publication behavior). Set to a
#'   positive value to override with a custom decay rate.
#' @param verbose Set to \code{TRUE} if you want chatty outputs; defaults to \code{FALSE}
#' @returns A \code{dat} file that contains both the total symptom scores at each timepoint
#'   and also all the individual factors that were used to generate those total scores
#' @examples
#'   #See vignettes for examples of how to use this
#' @export


generateData<-function(modelparam,respparam,blparam,trialdesign,empirical,makePositiveDefinite,seed=NA,lambda_cor=NA,verbose=FALSE){

  # Compute lambda_cor from carryover half-life if not specified
  if(is.na(lambda_cor)){
    if(modelparam$carryover_t1half > 0){
      lambda_cor <- log(2) / modelparam$carryover_t1half
    } else {
      lambda_cor <- 0
    }
  }

  # I. Turn the trial design information into something easier to use
  d<-data.table(trialdesign)
  d$t_wk_cumulative<-cumulative(d$t_wk)
  d[,onDrug:=(tod>0)]
  nP<-dim(trialdesign)[1]

  # II. Now set up a whole host of variables that we're going to want to track - essentially our baseline
  #     parameters for each subject ("bm","BL"), and the three modeled factors for each stage of the trial.

  # Set up the variable names
  cl<-c("tv","pb","br")
  labels<-c(c("bm","BL"),
            paste(trialdesign$timeptname,cl[1],sep="."),
            paste(trialdesign$timeptname,cl[2],sep="."),
            paste(trialdesign$timeptname,cl[3],sep="."))

  # Set up vectors with the standard deviations and means
  sds<-c(blparam[cat=="bm"]$sd,blparam[cat=="BL"]$sd)
  sds<-c(sds,rep(respparam[cat=="tv"]$sd,nP))
  sds<-c(sds,rep(respparam[cat=="pb"]$sd,nP)*trialdesign$e)
  sds<-c(sds,rep(respparam[cat=="br"]$sd,nP))
  means<-c(blparam[cat=="bm"]$m,blparam[cat=="BL"]$m)
  for (c in cl){
    rp<-respparam[cat==c]
    if(c=="tv"){
      means<-c(means,modgompertz(d$t_wk_cumulative,rp$max,rp$disp,rp$rate))
    }
    if(c=="pb"){
      means<-c(means,modgompertz(d$tpb,rp$max,rp$disp,rp$rate)*trialdesign$e)
    }
    if(c=="br"){
      brmeans<-modgompertz(d$tod,rp$max,rp$disp,rp$rate)
      if(nP>1){
        for(p in 2:nP){
          if(!d[p]$onDrug){
            if(d[p]$tsd>0){
              brmeans[p]<-brmeans[p]+brmeans[p-1]*(1/2)^(d$tsd[p]/modelparam$carryover_t1half)
            }
          }
        }
      }
      means<-c(means,brmeans)
    }
  }

  # Turn this into a matrix of correlations
  correlations<-diag(length(labels))
  rownames(correlations)<-labels
  colnames(correlations)<-labels
  for(c in cl){
    # AR(1) autocorrelations across time: rho^|t_i - t_j|
    if(nP>1){
      rho<-modelparam[[paste("c",c,sep=".")]]
      for(p in 1:(nP-1)){
        for(p2 in (1+p):nP){
          n1<-paste(trialdesign$timeptname[p],c,sep=".")
          n2<-paste(trialdesign$timeptname[p2],c,sep=".")
          time_gap<-abs(d$t_wk_cumulative[p2] - d$t_wk_cumulative[p])
          correlations[n1,n2]<-rho^time_gap
          correlations[n2,n1]<-rho^time_gap
        }
      }
    }
    # cross-factor correlations (same time and different time)
    for(c2 in setdiff(cl,c)){
      for(p in 1:nP){
        n1<-paste(trialdesign$timeptname[p],c,sep=".")
        n2<-paste(trialdesign$timeptname[p],c2,sep=".")
        correlations[n1,n2]<-modelparam$c.cf1t
        correlations[n2,n1]<-modelparam$c.cf1t
      }
      for(p in 1:(nP-1)){
        for(p2 in (1+p):nP){
          n1<-paste(trialdesign$timeptname[p],c,sep=".")
          n2<-paste(trialdesign$timeptname[p2],c2,sep=".")
          time_gap<-abs(d$t_wk_cumulative[p2] - d$t_wk_cumulative[p])
          correlations[n1,n2]<-modelparam$c.cfct * modelparam[[paste("c",c,sep=".")]]^time_gap
          correlations[n2,n1]<-modelparam$c.cfct * modelparam[[paste("c",c,sep=".")]]^time_gap
        }
      }
    }
    # biomarker-BR correlation with independent decay
    for(p in 1:nP){
      n1<-paste(trialdesign$timeptname[p],"br",sep=".")
      if(d[p]$onDrug){
        correlations[n1,'bm']<-modelparam$c.bm
        correlations['bm',n1]<-modelparam$c.bm
      } else if(d[p]$tsd>0 && lambda_cor>0){
        decay<-exp(-lambda_cor * d$tsd[p])
        correlations[n1,'bm']<-modelparam$c.bm * decay
        correlations['bm',n1]<-modelparam$c.bm * decay
      }
    }
  }
  # Turn correlation matrix into covariance matrix
  sigma<-outer(sds,sds)*correlations
  #  should be faster: outer(S,S)*R where S is SDs and R is correlation matrix
  #  previous: sigma<-correlations*sds%o%sds

  # Check and enforce positive definiteness:
  if(makePositiveDefinite){
    if(!is.positive.definite(sigma)){
      if(verbose){
        evals<-eigen(sigma, symmetric=TRUE, only.values=TRUE)$values
        neg_evals<-evals[evals<0]
        warning(sprintf(
          "Non-PD covariance matrix corrected (min eigenvalue: %.4f, %d negative). Consider constraining correlation parameters.",
          min(evals), length(neg_evals)))
      }
      sigma<-make.positive.definite(sigma, tol=1e-3)
    }
  }

  # Set seed:
  #if(is.na(seed)){seed<-round(rnorm(1,10000,10000))}
  #set.seed(seed)

  # Make the componant data!
  dat<-mvrnorm(n=modelparam$N,mu=means,Sigma=sigma,empirical=empirical)
  dat<-data.table(dat)
  setnames(dat,names(dat),labels)

  # TESTING
  #tdat1<-dat[, lapply(.SD, mean)]
  #plot(c(d$t_wk_cumulative),tdat1[,.(OL.br,BD.br,UD.br,COd.br,COp.br)])


  #tdat2<-dat[, lapply(.SD, mean)]
  #plot(c(d$t_wk_cumulative),tdat2[,.(OL.br,BD.br,UD.br,COd.br,COp.br)])

  # Turn it into actual fake data... calc intervals separately, because will use
  dat[,ptID:=1:modelparam$N]
  for(p in 1:nP){
    n<-trialdesign$timeptname[p]
    evalstring<-paste(
      "dat[,D_",n,":=sum(",paste(paste(n,cl,sep='.',collapse=','),sep=','),"),by='ptID']",
      sep="")
    eval(parse(text=evalstring))
  }
  for(p in 1:nP){
    n<-trialdesign$timeptname[p]
    evalstring<-paste(
      "dat[,",n,":=sum(BL,-",paste(paste(n,cl,sep='.',collapse=',-'),sep=','),"),by='ptID']",
      sep="")
    eval(parse(text=evalstring))
  }

  # TESTING
  #tdat3<-dat[, lapply(.SD, mean)]
  #plot(c(0,d$t_wk_cumulative),tdat3[,.(BL,OL,BD,UD,COd,COp)])
  #plot(c(d$t_wk_cumulative),tdat3[,.(D_OL,D_BD,D_UD,D_COd,D_COp)])

  return(dat)
}
