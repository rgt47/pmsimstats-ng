# Short helper functions

#' Cumulative
#'
#' \code{cumulative} turns a vector of intervals into a cumulative running time
#'
#' @export
cumulative<-function(x){
  I<-length(x)
  y<-rep(NA,I)
  y[1]<-x[1]
  if(I>1){
    for(i in 2:I){
      y[i]<-x[i]+y[i-1]
    }
  }
  return(y)
}

#' Modified Gompertz
#'
#' \code{modgompertz} takes response parameters and returns change in symptoms as a function of time
#'
#' This is a slight modification of a standard gompertz function that is designed to model
#' how symptoms can change as a function of time.
#'
#' @param t time - no specific units, just make sure it's consistant across your modeling
#' @param maxr maximum response that can be attained (ceiling)
#' @param disp vertical displacement
#' @param rate rate
#' @return A numeric value indicating how much "improvement" there has been since baseline
#'   as of point \code{t} in time. A negative value indicates worsening rather than improvement.
#' @export
modgompertz<-function(t,maxr,disp,rate){
  # Guard: maxr=0 implies y identically zero; otherwise the final
  # rescale divides by (maxr - vert_offset) = 0 and returns NaN.
  if (isTRUE(all.equal(maxr, 0))) return(rep(0, length(t)))
  y<-maxr*exp(-disp*exp(-rate*t))
  vert_offset<-maxr*exp(-disp*exp(-rate*0))
  y<-y-vert_offset
  y<-y*(maxr/(maxr-vert_offset)) # return the assymptotic max to what it was
  return(y)
}

# -----------------------------------------------------------------
# Trajectory family dispatch
# -----------------------------------------------------------------
#
# Three additional monotone-saturating families share modGompertz's
# vertical-offset convention, that is, y(0)=0 and y(t -> infty) = maxr.
# Each family is a 3-parameter form (matching modGompertz's three free
# parameters: maxr, disp/rate-equivalent, rate/breakpoint-equivalent).
# Calibration of the secondary parameters to match cross-family asymp-
# tote and rise time is the caller's responsibility (see, e.g.,
# analysis/scripts/gompertz-evaluation/01-architecture-a-trajectory-
# sweep.R for the canonical week-8 anchor calibration).
#
# Logistic: y(t) = maxr / (1 + exp(-log_rate*(t - t0))), with vertical
# offset subtracted and asymptote rescaled so y(0)=0 and y(infty)=maxr,
# in direct analogy to modGompertz. Free parameters: maxr, log_rate, t0.
#
# Hyperbolic tangent: y(t) = (maxr/2)*(1 + tanh(rate*(t - t0))), with
# vertical offset adjustment. Free parameters: maxr, rate, t0. The
# dispersion-equivalent here is t0, which controls the inflection
# location and therefore the rise sharpness for fixed rate.
#
# Piecewise linear with breakpoint: a linear ramp from 0 to maxr over
# [0, t_breakpoint], then a slow post-breakpoint slope. Free parameters:
# maxr, t_breakpoint, post_slope. The smooth limit (post_slope=0) is a
# pure ramp-then-flat function. y(0)=0 already, no offset adjustment.

#' Logistic trajectory (vertical-offset adjusted)
#'
#' Logistic curve with the same vertical-offset convention as
#' \link{modgompertz}: y(0)=0 and y(t->infty)=maxr.
#'
#' @param t time
#' @param maxr asymptotic maximum response
#' @param log_rate steepness parameter (analog of Gompertz rate)
#' @param t0 inflection-point time (analog of Gompertz dispersion)
#' @return numeric vector of length(t)
#' @export
shape_logistic <- function(t, maxr, log_rate, t0) {
  if (isTRUE(all.equal(maxr, 0))) return(rep(0, length(t)))
  y <- maxr / (1 + exp(-log_rate * (t - t0)))
  vert_offset <- maxr / (1 + exp(-log_rate * (0 - t0)))
  y <- y - vert_offset
  y <- y * (maxr / (maxr - vert_offset))
  y
}

#' Hyperbolic-tangent trajectory (vertical-offset adjusted)
#'
#' Hyperbolic tangent curve with the same vertical-offset convention
#' as \link{modgompertz}: y(0)=0 and y(t->infty)=maxr.
#'
#' @param t time
#' @param maxr asymptotic maximum response
#' @param rate steepness (analog of Gompertz rate)
#' @param t0 inflection-point time (analog of Gompertz dispersion)
#' @return numeric vector of length(t)
#' @export
shape_hyperbolic_tangent <- function(t, maxr, rate, t0) {
  if (isTRUE(all.equal(maxr, 0))) return(rep(0, length(t)))
  y <- (maxr / 2) * (1 + tanh(rate * (t - t0)))
  vert_offset <- (maxr / 2) * (1 + tanh(rate * (0 - t0)))
  y <- y - vert_offset
  y <- y * (maxr / (maxr - vert_offset))
  y
}

#' Piecewise-linear breakpoint trajectory
#'
#' A linear ramp from 0 at t=0 to maxr at t=t_breakpoint, then a
#' (typically mild) post-breakpoint slope. y(0)=0 by construction.
#'
#' @param t time
#' @param maxr value attained at the breakpoint
#' @param t_breakpoint week at which the rise saturates
#' @param post_slope post-breakpoint slope (default 0 = flat)
#' @return numeric vector of length(t)
#' @export
shape_piecewise_linear_breakpoint <- function(t, maxr, t_breakpoint,
                                              post_slope = 0) {
  if (isTRUE(all.equal(maxr, 0))) return(rep(0, length(t)))
  ramp <- maxr * pmin(1, t / t_breakpoint)
  tail <- post_slope * pmax(0, t - t_breakpoint)
  ramp + tail
}

#' Trajectory family dispatch
#'
#' Single entry point that returns a trajectory of the requested
#' family. Each family has the modGompertz vertical-offset convention
#' (y(0)=0 and y(t->infty)=maxr) except piecewise_linear_breakpoint,
#' which is exactly zero at t=0 by construction.
#'
#' @param family one of "gompertz", "logistic", "hyperbolic_tangent",
#'   "piecewise_linear_breakpoint"
#' @param t time vector
#' @param maxr asymptotic maximum response
#' @param p2,p3 family-specific shape parameters (see Details)
#' @details
#' For "gompertz", \code{p2 = disp} and \code{p3 = rate}.
#' For "logistic", \code{p2 = log_rate} and \code{p3 = t0}.
#' For "hyperbolic_tangent", \code{p2 = rate} and \code{p3 = t0}.
#' For "piecewise_linear_breakpoint", \code{p2 = t_breakpoint} and
#'   \code{p3 = post_slope}.
#' @return numeric vector of length(t)
#' @export
trajectoryShape <- function(family, t, maxr, p2, p3) {
  switch(family,
    gompertz                    = modgompertz(t, maxr, p2, p3),
    logistic                    = shape_logistic(t, maxr, p2, p3),
    hyperbolic_tangent          = shape_hyperbolic_tangent(t, maxr, p2, p3),
    piecewise_linear_breakpoint = shape_piecewise_linear_breakpoint(
                                    t, maxr, p2, p3),
    stop(sprintf("Unknown trajectory family: '%s'", family))
  )
}

#' Reknit simulated results
#'
#' \code{reknitsimresults} recombines the results of simulations that have the
#' exact same parameter space set, but were run as separate chunks - e.g., you
#' ran 100 repititions, realized you needed more, and ran 150 more to have 250
#' total. Does not work with rawdataout.
#'
#' @param basesavename The basesavename as described in \link{generateSimulatedResults}
#' @param savedir The directory the files are saved in
#' @return The output is a file just like you would have received from
#'   \link{generateSimulatedResults} if it had all been run at once.
#' @export
reknitsimresults<-function(basesavename,savedir){
  if(missing(savedir)){
    initialdirectory<-getwd()
  }else{
    initialdirectory<-getwd()
    setwd(savedir)
  }
  sr<-vector(mode="list",length=length(basesavename))
  for(iSN in 1:length(basesavename)){
    simresults<-readRDS(paste(basesavename[1],"1",sep="_save"))
    nextfilenum<-2
    while(nextfilenum>0){
      nextfile<-paste(basesavename[1],nextfilenum,sep="_save")
      if(file.exists(nextfile)){
        newresults<-readRDS(nextfile)
        simresults$results<-rbind(simresults$results,newresults$results)
        if(length(names(simresults))>2){
          warning("Haven't implemented adding rawdata together yet")
        }
        nextfilenum<-nextfilenum+1
      }else{
        nextfilenum<-0
      }
    }
    sr[[iSN]]<-simresults
  }
  if(length(sr)>1){
    # test the parameter settings are the same
    for(iSN in 2:length(sr)){
      if(!identical(sr[[1]]$parameterselections$trialdesigns,
          sr[[iSN]]$parameterselections$trialdesigns)){
        warning("Trial design settings don't match")
      }
      if(!identical(sr[[1]]$parameterselections$respparamsets,
                    sr[[iSN]]$parameterselections$respparamsets)){
        warning("Response paramset settings don't match")
      }
      if(!identical(sr[[1]]$parameterselections$blparamsets,
                    sr[[iSN]]$parameterselections$blparamsets)){
        warning("Baseline paramset settings don't match")
      }
    }
    # Merge (inefficiently, if had saved differently could use rbindlist...)
    simresults<-sr[[1]]
    for(iSN in 2:length(sr)){
      simresults$results<-rbind(simresults$results,sr[[iSN]]$results)
    }
  }else{
    simresults<-sr[[1]]
  }
  return(simresults)
}
