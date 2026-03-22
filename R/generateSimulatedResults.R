#' Generate Simulated Results
#'
#' \code{generateSimulatedResults} rotates through a parameter space, creating
#'   and analyizing simulated data and collating the results
#'
#' This is a wrapper for \link{generateData} and \link{lme_analysis} that
#' systematically works its way through a set of parameters and creates
#' a specified number of repetitions of simulated data and results of
#' analyizing that data for each set of parameter specifications.
#'
#' @param trialdesigns A list of trial designs of the form generated
#'   by \link{buildtrialdesigns}
#' @param respparamsets A list of options for \code{respparam} as
#'   described in \link{generateData}
#' @param blparamsets A list of options for \code{blparam} as
#'   described in \link{generateData}
#' @param censorparams A data.table where each row is a set of options for
#'   \code{censorparam} as described in \link{censordata}
#' @param modelparams A data.table where each row is a set of options for
#'   \code{modelparams} as described in \link{generateData}
#' @param simparam The options that control the logistics of the simulations:
#'   \itemize{
#'     \item{\code{Nreps}}{  Number of repititions for a given set of parameters}
#'     \item{\code{progressiveSave=TRUE}}{  Do you want to save chunks as you go?
#'       Helpful if the run is going to take a long time, and you don't want the
#'       run to be interrupted 8 hours in to a 9 hour run and lose everything.}
#'     \item{\code{basesavename}}{  If you are using progressivesave, what do you
#'       want the base of the filename to be?  Will have _save1, _save2, etc appended}
#'     \item{\code{nRep2save}}{  How many parameter sets to you want to run at
#'       a time before saving?}
#'     \item{\code{saveunit2start}}{  If you got through part of your parameter
#'       space but not all of it, you can start partway through. If starting
#'       from the beginning, set to 1}
#'     \item{\code{savedir}}{  What directory do you want to save to?}
#'   }
#' @param analysisparams The options used by \link{lme_analysis}. Has the
#'   following components:
#'   \itemize{
#'     \item{\code{useDE=TRUE}}{  Binary: should the LME model use expectnacy
#'       information if available?}
#'     \item{\code{t_random_slope=FALSE}}{  Binary: should the LME model use
#'       random slopes, as well as random intercepts, for participants?}
#'   }
#' @param rawdataout=FALSE Binary - output the raw data for all
#'   simulations?  This is memory intensive for large runs.
#' @param lambda_cor Correlation decay rate (NA for auto, see
#'   \link{generateData})
#' @return Returns a list with three named parts:
#'   \itemize{
#'     \item{\code{results}}{  A large data table that tells you the parameters used
#'       and the analysis results for every simulated trial, one row per simualated trial.}
#'     \item{\code{parameterselections}}{  A list desribing the parameterspace used.
#'       This is important because some of the information in \code{out$results} is
#'       the index of the parameter set selection, and you will need this additional
#'       information to interpret it. For example, \code{trialdesigns}, \code{respparamset},
#'       \code{blparamset}, and \code{modelparamset} are all indexed.}
#'     \item{\code{rawdataout}}{  If you have rawdataout set to \code{TRUE}, you will
#'       also get (1) \code{drawdataout$precensor}, a list of simulated trial data of
#'       a length corresponding to the number of parameter permulations, with the order
#'       parallel to that of the main output results. Repititions of a particular
#'       parameter set are mixed, identified by \code{irep}. And, \code{rawdataout$postcensor},
#'       which has the same format but with the psotcensoring data}
#'   }
#' @examples
#'   # See vignettes for examples
#' @export

generateSimulatedResults<-function(trialdesigns,respparamsets,blparamsets,
                                   censorparams,modelparams,simparam,analysisparams,
                                   rawdataout=FALSE,lambda_cor=NA,n_cores=1){

  # defaults
  if(missing(analysisparams)) analysisparams<-list(useDE=TRUE,
                                                   t_random_slope=FALSE,
                                                   full_model_out=FALSE)

  # Auto-detect cores if requested
  if(n_cores < 1) n_cores <- max(1, parallel::detectCores() - 1)

  # Hold original directory in case we change it
  initialdirectory<-getwd()

  # Setup parameter grid
  tic("Total time:")
  VPGmaster<-expand.grid(
    trialdesign=1:length(trialdesigns),
    respparamset=1:length(respparamsets),
    blparamset=1:length(blparamsets),
    modelparamset=1:dim(modelparams)[1]
  )
  nV<-dim(VPGmaster)[1]
  irep<-1
  totrep<-nV*simparam$Nreps
  totparamsets<-nV

  # --- Optimization 1: Pre-build and cache sigma matrices ---
  cat("Caching sigma matrices...\n")
  tic("Sigma cache built")
  sigma_cache<-list()
  for(iR in 1:nV){
    td<-trialdesigns[[VPGmaster[iR,"trialdesign"]]][[2]]
    rp<-respparamsets[[VPGmaster[iR,"respparamset"]]]$param
    bp<-blparamsets[[VPGmaster[iR,"blparamset"]]]$param
    mp<-modelparams[VPGmaster[iR,"modelparamset"],]
    nP<-length(td)
    for(iP in 1:nP){
      cache_key<-paste(VPGmaster[iR,"trialdesign"],
                       VPGmaster[iR,"respparamset"],
                       VPGmaster[iR,"blparamset"],
                       VPGmaster[iR,"modelparamset"],
                       iP, sep="_")
      if(is.null(sigma_cache[[cache_key]])){
        sigma_cache[[cache_key]]<-buildSigma(
          mp, rp, bp, td[[iP]],
          makePositiveDefinite=TRUE,
          lambda_cor=lambda_cor, verbose=FALSE)
      }
    }
  }
  toc()
  cat(sprintf("Cached %d unique sigma matrices\n\n",
              length(sigma_cache)))

  # Progressive save setup
  if(simparam$progressiveSave==FALSE){
    nLargeLoops<-1
    LLstarts<-1
    LLstops<-nV
  }else{
    nLargeLoops<-ceiling(nV/simparam$nRep2save)
    if(nLargeLoops>1){
      LLstarts<-c(1,1+simparam$nRep2save*(1:(nLargeLoops-1)))
      LLstops<-c(LLstarts[2:nLargeLoops]-1,nV)
    }else{
      LLstarts<-1
      LLstops<-nV
    }
  }

  for(iLL in simparam$saveunit2start:nLargeLoops){
    tic(paste("***Progressive Save Unit ",(iLL+1-simparam$saveunit2start)," of ",(nLargeLoops+1-simparam$saveunit2start),
              " complete",sep=""))
    VPG<-VPGmaster[LLstarts[iLL]:LLstops[iLL],]
    iparamset<-LLstarts[iLL]

    if(rawdataout){d_out<-list(precensor=list(),postcensor=list())}

    # Define per-parameter-set worker function
    run_one_paramset<-function(iR, VPG, trialdesigns, respparamsets,
                               blparamsets, modelparams, simparam,
                               analysisparams, censorparams,
                               sigma_cache, nV, iLL_offset){
      td<-trialdesigns[[VPG[iR,"trialdesign"]]][[2]]
      pp<-list()
      pp$respparam<-respparamsets[[VPG[iR,"respparamset"]]]$param
      pp$blparam<-blparamsets[[VPG[iR,"blparamset"]]]$param
      pp$modelparam<-modelparams[VPG[iR,"modelparamset"],]

      nP<-length(td)
      Ns<-pp$modelparam$N%/%nP
      Ns<-Ns+c(rep(1,pp$modelparam$N%%nP),rep(0,nP-pp$modelparam$N%%nP))
      NNs<-Ns*simparam$Nreps

      pp$modelparam$N<-NNs[[1]]
      ck<-paste(VPG[iR,"trialdesign"], VPG[iR,"respparamset"],
                VPG[iR,"blparamset"], VPG[iR,"modelparamset"],
                1, sep="_")
      dat<-generateData(pp$modelparam,pp$respparam,pp$blparam,td[[1]],FALSE,TRUE,
                        cached_sigma=sigma_cache[[ck]])
      dat[,path:=1]
      dat[,replicate:=rep(1:simparam$Nreps,Ns[1])]
      if(nP>1){
        for(iP in 2:nP){
          pp$modelparam$N<-NNs[[iP]]
          ck<-paste(VPG[iR,"trialdesign"], VPG[iR,"respparamset"],
                    VPG[iR,"blparamset"], VPG[iR,"modelparamset"],
                    iP, sep="_")
          dat2<-generateData(pp$modelparam,pp$respparam,pp$blparam,td[[iP]],FALSE,TRUE,
                             cached_sigma=sigma_cache[[ck]])
          dat2[,path:=iP]
          dat2[,replicate:=rep(1:simparam$Nreps,Ns[iP])]
          dat<-rbind(dat,dat2)
        }
      }
      pp$modelparam$N<-sum(Ns)

      results_list<-list()
      ridx<-0
      for(iAP in 1:dim(analysisparams)[1]){
        ETanalysisout<-lme_analysis(td,dat,analysisparams[iAP,])
        for(iS in 1:simparam$Nreps){
          analysisout<-lme_analysis(td,dat[replicate==iS],analysisparams[iAP,])
          ridx<-ridx+1
          results_list[[ridx]]<-cbind(VPG[iR,],as.data.table(pp$modelparam),
                   data.table(censorparamset=0,use_DE=analysisparams[iAP,]$useDE,
                              t_random_slope=analysisparams[iAP,]$t_random_slope,
                              irep=((iR-1)*simparam$Nreps+iS),frac_NA=0,
                              ETbeta=ETanalysisout$beta,ETbetaSE=ETanalysisout$betaSE,
                              beta=analysisout$beta,betaSE=analysisout$betaSE,
                              p=analysisout$p,issingular=analysisout$issingular,
                              warning=analysisout$warning))
        }
        nocensoringflag<-FALSE
        if(length(censorparams)<2){
          if(is.na(censorparams)) nocensoringflag<-TRUE
        }
        if(!nocensoringflag){
          for(iC in 1:dim(censorparams)[1]){
            datc<-censordata(dat,td[[1]],censorparams[iC,])
            for(iS in 1:simparam$Nreps){
              frac_NA<-sum(is.na(datc[replicate==iS]))/(pp$modelparam$N*length(td[[1]]$t_wk))
              analysisout<-lme_analysis(td,datc[replicate==iS],analysisparams[iAP,])
              ridx<-ridx+1
              results_list[[ridx]]<-cbind(VPG[iR,],as.data.table(pp$modelparam),
                       data.table(censorparamset=iC,use_DE=analysisparams[iAP,]$useDE,
                                  t_random_slope=analysisparams[iAP,]$t_random_slope,
                                  irep=((iR-1)*simparam$Nreps+iS),frac_NA=frac_NA,
                                  ETbeta=ETanalysisout$beta,ETbetaSE=ETanalysisout$betaSE,
                                  beta=analysisout$beta,betaSE=analysisout$betaSE,
                                  p=analysisout$p,issingular=analysisout$issingular,
                                  warning=analysisout$warning))
            }
          }
        }
      }
      rbindlist(results_list)
    }

    # Run parameter sets (parallel or sequential)
    nR<-dim(VPG)[1]
    use_parallel<-(n_cores > 1) && requireNamespace("furrr", quietly=TRUE)

    if(use_parallel){
      future::plan(future::multisession, workers=n_cores)
      cat(sprintf("Parallel: %d workers for %d parameter sets\n",
                  n_cores, nR))

      # Wrap worker to carry function definitions via closure
      make_worker<-function(genData_fn, lmeAnal_fn, censor_fn,
                            worker_fn){
        function(iR){
          # Inject functions into worker's scope
          environment(worker_fn)$generateData<-genData_fn
          environment(worker_fn)$lme_analysis<-lmeAnal_fn
          environment(worker_fn)$censordata<-censor_fn
          worker_fn(iR, VPG, trialdesigns, respparamsets,
                    blparamsets, modelparams, simparam,
                    analysisparams, censorparams,
                    sigma_cache, nV, LLstarts[iLL])
        }
      }
      worker<-make_worker(generateData, lme_analysis,
                          censordata, run_one_paramset)

      out_parts<-furrr::future_map(1:nR, worker,
        .options=furrr::furrr_options(
          seed=TRUE,
          packages=c("data.table","nlme","MASS","corpcor")),
        .progress=TRUE)

      out<-rbindlist(out_parts)
      future::plan(future::sequential)
    } else {
      out_parts<-vector("list", nR)
      for(iR in 1:nR){
        tic(paste("Parameter set ",iR," of ",nR,
                  " (", iparamset," of ",nV," total)",sep=""))
        out_parts[[iR]]<-run_one_paramset(
          iR, VPG, trialdesigns, respparamsets,
          blparamsets, modelparams, simparam,
          analysisparams, censorparams,
          sigma_cache, nV, LLstarts[iLL])
        iparamset<-iparamset+1
        toc()
      }
      out<-rbindlist(out_parts)
    }

    # Organize the output and return it
    if(rawdataout){
      outpt<-list(results=out,
                  parameterselections=list(trialdesigns=trialdesigns,
                                           respparamsets=respparamsets,
                                           blparamsets=blparamsets,
                                           censorparams=censorparams,
                                           modelparams=modelparams,
                                           analysisparams=analysisparams,
                                           simparam=simparam),
                  rawdata=d_out)
    }else{
      outpt<-list(results=out,
                  parameterselections=list(trialdesigns=trialdesigns,
                                           respparamsets=respparamsets,
                                           blparamsets=blparamsets,
                                           censorparams=censorparams,
                                           modelparams=modelparams,
                                           analysisparams=analysisparams,
                                           simparam=simparam))
    }

    if(simparam$progressiveSave){
      setwd(simparam$savedir)
      saveRDS(outpt,paste(simparam$basesavename,iLL,sep="_save"))
    }
    toc()
  } # end iLL
  toc()
  setwd(initialdirectory)
  if(!simparam$progressiveSave) return(outpt)
}
