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
#'  \item{\code{c.bm}}{  Biomarker-BR parameter. Under MVN architecture
#'    (default), this is the correlation between biomarker and biological
#'    response. Under mean moderation architecture, this is the regression
#'    coefficient scaling the BR component by centered biomarker value.}
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
#' @param cached_sigma Pre-built sigma matrix (list with sigma, means,
#'   labels, nP, cl, trialdesign components). When provided, skips
#'   sigma construction entirely. Use \link{buildSigma} to create.
#' @param dgp_architecture DGP architecture for the biomarker-treatment
#'   interaction. \code{"mvn"} (default) uses differential correlation
#'   in the covariance structure (Architecture B). \code{"mean_moderation"}
#'   uses direct mean scaling of the BR component by the biomarker
#'   (Architecture A). See the DGP architecture white paper for details.
#' @returns A \code{dat} file that contains both the total symptom scores at each timepoint
#'   and also all the individual factors that were used to generate those total scores
#' @examples
#'   #See vignettes for examples of how to use this
#' @export


generateData<-function(modelparam,respparam,blparam,trialdesign,empirical,makePositiveDefinite,seed=NA,lambda_cor=NA,verbose=FALSE,cached_sigma=NULL,dgp_architecture="mvn"){

  dgp_architecture<-match.arg(dgp_architecture, c("mvn", "mean_moderation"))

  if(!is.null(cached_sigma)){
    sigma<-cached_sigma$sigma
    means<-cached_sigma$means
    labels<-cached_sigma$labels
    nP<-cached_sigma$nP
    cl<-cached_sigma$cl
    chol_sigma<-cached_sigma$chol_sigma
  } else {
    built<-buildSigma(modelparam,respparam,blparam,trialdesign,makePositiveDefinite,lambda_cor,verbose,
                      dgp_architecture=dgp_architecture)
    sigma<-built$sigma
    means<-built$means
    labels<-built$labels
    nP<-built$nP
    cl<-built$cl
    chol_sigma<-built$chol_sigma
    trialdesign<-built$trialdesign
  }

  # Draw participants using pre-computed Cholesky factor
  # Equivalent to mvrnorm(n, mu, Sigma) but avoids
  # recomputing the Cholesky decomposition on each call
  n<-modelparam$N
  p<-length(means)
  Z<-matrix(rnorm(n * p), nrow=n, ncol=p)
  dat<-Z %*% chol_sigma + matrix(means, nrow=n, ncol=p, byrow=TRUE)
  dat<-data.table(dat)
  setnames(dat,names(dat),labels)

  # Architecture A: additive biomarker moderation of BR
  # Each participant's BR is shifted proportionally to their
  # standardized biomarker value when on drug. The shift magnitude
  # is c.bm * sigma_br per SD of biomarker, matching the
  # conditional expectation under Architecture B's MVN model:
  #   E[BR | bm, on_drug] = mu_BR + c.bm * (sigma_br/sigma_bm) * (bm - mu_bm)
  # This produces comparable effect sizes between architectures
  # for the same c.bm value.
  if(dgp_architecture == "mean_moderation"){
    beta_bm<-modelparam$c.bm
    bm_mean<-blparam[cat=="bm"]$m
    bm_sd<-blparam[cat=="bm"]$sd
    br_sd<-respparam[cat=="br"]$sd
    bm_z<-(dat$bm - bm_mean) / bm_sd

    d<-data.table(trialdesign)
    d[,onDrug:=(tod>0)]
    tnames<-if(is.data.frame(trialdesign)) trialdesign$timeptname else trialdesign$timeptnames
    if(is.null(tnames)) tnames<-trialdesign$timeptname

    for(tp in 1:nP){
      if(d[tp]$onDrug){
        br_col<-paste(tnames[tp], "br", sep=".")
        dat[,(br_col):=get(br_col) + beta_bm * bm_z * br_sd]
      }
    }
  }

  # Compute outcome columns using vectorized operations
  dat[,ptID:=1:modelparam$N]
  tnames<-if(is.data.frame(trialdesign)) trialdesign$timeptname else trialdesign$timeptnames
  if(is.null(tnames)) tnames<-trialdesign$timeptname
  for(p in 1:nP){
    n<-tnames[p]
    comps<-paste(n,cl,sep=".")
    dat[,(paste0("D_",n)):=rowSums(.SD),.SDcols=comps]
  }
  for(p in 1:nP){
    n<-tnames[p]
    comps<-paste(n,cl,sep=".")
    dat[,(n):=BL - rowSums(.SD),.SDcols=comps]
  }

  return(dat)
}


#' Build sigma matrix for a trial path
#'
#' \code{buildSigma} constructs the covariance matrix, mean vector,
#' and labels for a single trial path. The result can be cached and
#' passed to \link{generateData} via the cached_sigma parameter to
#' avoid redundant matrix construction across replicates.
#'
#' @param modelparam Model parameters (see \link{generateData})
#' @param respparam Response parameters (see \link{generateData})
#' @param blparam Baseline parameters (see \link{generateData})
#' @param trialdesign Trial path specification
#' @param makePositiveDefinite Force positive definiteness?
#' @param lambda_cor Correlation decay rate (NA for auto)
#' @param verbose Print diagnostics?
#' @param dgp_architecture DGP architecture (\code{"mvn"} or
#'   \code{"mean_moderation"}). Under mean moderation, the BM-BR
#'   correlation block is omitted from the covariance matrix.
#' @return List with sigma, means, labels, nP, cl, trialdesign
#' @export
buildSigma<-function(modelparam,respparam,blparam,trialdesign,makePositiveDefinite=TRUE,lambda_cor=NA,verbose=FALSE,dgp_architecture="mvn"){

  # Compute lambda_cor from carryover half-life if not specified
  if(is.na(lambda_cor)){
    if(modelparam$carryover_t1half > 0){
      lambda_cor <- log(2) / modelparam$carryover_t1half
    } else {
      lambda_cor <- 0
    }
  }

  d<-data.table(trialdesign)
  d$t_wk_cumulative<-cumulative(d$t_wk)
  d[,onDrug:=(tod>0)]
  nP<-dim(trialdesign)[1]

  cl<-c("tv","pb","br")
  labels<-c(c("bm","BL"),
            paste(trialdesign$timeptname,cl[1],sep="."),
            paste(trialdesign$timeptname,cl[2],sep="."),
            paste(trialdesign$timeptname,cl[3],sep="."))

  sds<-c(blparam[cat=="bm"]$sd,blparam[cat=="BL"]$sd)
  sds<-c(sds,rep(respparam[cat=="tv"]$sd,nP))
  sds<-c(sds,rep(respparam[cat=="pb"]$sd,nP)*trialdesign$e)
  sds<-c(sds,rep(respparam[cat=="br"]$sd,nP))
  means<-c(blparam[cat=="bm"]$m,blparam[cat=="BL"]$m)
  for (cc in cl){
    rp<-respparam[cat==cc]
    if(cc=="tv"){
      means<-c(means,modgompertz(d$t_wk_cumulative,rp$max,rp$disp,rp$rate))
    }
    if(cc=="pb"){
      means<-c(means,modgompertz(d$tpb,rp$max,rp$disp,rp$rate)*trialdesign$e)
    }
    if(cc=="br"){
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

  # Build correlation matrix
  nLabels<-length(labels)
  correlations<-diag(nLabels)
  rownames(correlations)<-labels
  colnames(correlations)<-labels

  # Pre-compute time gaps for AR(1)
  weeks<-d$t_wk_cumulative

  for(cc in cl){
    rho<-modelparam[[paste("c",cc,sep=".")]]
    # AR(1) autocorrelations
    if(nP>1){
      for(p in 1:(nP-1)){
        for(p2 in (1+p):nP){
          n1<-paste(trialdesign$timeptname[p],cc,sep=".")
          n2<-paste(trialdesign$timeptname[p2],cc,sep=".")
          tg<-abs(weeks[p2] - weeks[p])
          ar1_val<-rho^tg
          correlations[n1,n2]<-ar1_val
          correlations[n2,n1]<-ar1_val
        }
      }
    }
    # Cross-factor correlations
    for(c2 in setdiff(cl,cc)){
      for(p in 1:nP){
        n1<-paste(trialdesign$timeptname[p],cc,sep=".")
        n2<-paste(trialdesign$timeptname[p],c2,sep=".")
        correlations[n1,n2]<-modelparam$c.cf1t
        correlations[n2,n1]<-modelparam$c.cf1t
      }
      if(nP>1){
        for(p in 1:(nP-1)){
          for(p2 in (1+p):nP){
            n1<-paste(trialdesign$timeptname[p],cc,sep=".")
            n2<-paste(trialdesign$timeptname[p2],c2,sep=".")
            tg<-abs(weeks[p2] - weeks[p])
            cfval<-modelparam$c.cfct * rho^tg
            correlations[n1,n2]<-cfval
            correlations[n2,n1]<-cfval
          }
        }
      }
    }
    # BM-BR correlation with decay (Architecture B / MVN only)
    if(cc=="br" && dgp_architecture=="mvn"){
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
  }

  # Covariance matrix
  sigma<-outer(sds,sds)*correlations

  # PD check (skip eigenvalue computation unless verbose)
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

  # Pre-compute Cholesky factor for efficient MVN draws
  chol_sigma<-tryCatch(chol(sigma), error=function(e){
    chol(make.positive.definite(sigma, tol=1e-3))
  })

  list(sigma=sigma, means=means, labels=labels, nP=nP, cl=cl,
       trialdesign=trialdesign, chol_sigma=chol_sigma)
}


#' Validate parameter grid for positive definiteness
#'
#' Tests all parameter combinations for PD before simulation.
#' Returns only valid combinations.
#'
#' @param trialdesigns List of trial designs
#' @param modelparams Model parameter grid
#' @param respparam Response parameters
#' @param blparam Baseline parameters
#' @param lambda_cor Correlation decay rate
#' @return List with valid_params (filtered grid) and n_rejected
#' @export
validateParameterGrid<-function(trialdesigns,modelparams,respparam,blparam,lambda_cor=NA){
  n_total<-0
  n_rejected<-0
  rejected<-character()

  for(iTD in 1:length(trialdesigns)){
    td<-trialdesigns[[iTD]][[2]]
    for(iP in 1:length(td)){
      for(iM in 1:nrow(modelparams)){
        n_total<-n_total+1
        mp<-modelparams[iM,]
        result<-tryCatch({
          s<-buildSigma(mp,respparam,blparam,td[[iP]],
                        makePositiveDefinite=FALSE,lambda_cor=lambda_cor)
          is.positive.definite(s$sigma)
        }, error=function(e) FALSE)
        if(!result){
          n_rejected<-n_rejected+1
          td_name<-trialdesigns[[iTD]][[1]]$name_shortform
          rejected<-c(rejected, sprintf(
            "%s path %d: c.bm=%.2f co=%.1f",
            td_name, iP, mp$c.bm, mp$carryover_t1half))
        }
      }
    }
  }

  if(n_rejected > 0){
    cat(sprintf("PD validation: %d/%d rejected (%.1f%%)\n",
                n_rejected, n_total, 100*n_rejected/n_total))
    for(r in rejected) cat("  ", r, "\n")
  } else {
    cat(sprintf("PD validation: all %d combinations valid\n", n_total))
  }

  list(n_total=n_total, n_rejected=n_rejected, rejected=rejected)
}
