repel.sim.capt <- function(fit = NULL, traps = NULL, mask = NULL, popn = NULL,
                     infotypes = character(0), detfn = "hn",
                     pars = NULL, ss.opts = NULL, cue.rates = NULL, survey.length = NULL,
                     freq.dist = "edf", sound.speed = 330, test.detfn = FALSE,
                     first.only = FALSE, keep.locs = FALSE, keep.ids = FALSE, ...){
  arg.names <- names(as.list(environment()))
  extra.args <- list(...)
  if (any(names(extra.args) == "call.freqs")){
    if (!missing(cue.rates)){
      stop("The argument `cue.rates' has replaced `call.freqs'; use only the former.")
    }
    warning("The argument `call.freqs' is deprecated; please rename to `cue.rates' instead.")
    cue.rates <- extra.args[["call.freqs"]]
  }
  ## Some error checking.
  if (any(infotypes == "ss")){
    stop("Signal strength information is simulated by setting argument 'detfn' to \"ss\".")
  }
  if (!missing(ss.opts) & detfn != "ss"){
    warning("The argument 'ss.opts' is being ignored, as 'detfn' is not \"ss\".")
  }
  if (keep.ids & is.null(cue.rates)){
    warning("IDs are necessarily different as only one call is simulated from each individual.")
  }
  ## Warnings for ignored parameters. Is there a neater way of doing
  ## this? The 'missing()' function is annoying.
  if (!missing(mask) & !missing(fit)){
    warning("Argument 'mask' is being ignored as 'fit' was provided.")
  }
  if (!missing(infotypes) & !missing(fit)){
    warning("Argument 'infotypes' is being ignored as 'fit' was provided.")
  }
  if (!missing(detfn) & !missing(fit)){
    warning("Argument 'detfn' is being ignored as 'fit' was provided.")
  }
  if (!missing(traps) & !missing(fit)){
    warning("Argument 'traps' is being ignored as 'fit' was provided.")
  }
  if (!missing(pars) & !missing(fit)){
    warning("Argument 'pars' is being ignored as 'fit' was provided.")
  }
  if (!missing(ss.opts) & !missing(fit)){
    warning("Argument 'ss.opts' is being ignored as 'fit' was provided.")
  }
  if (!missing(cue.rates) & !missing(fit)){
    warning("Argument 'cue.rates' is being ignored as 'fit' was provided.")
  }
  if (!missing(sound.speed) & !missing(fit)){
    warning("Argument 'sound.speed' is being ignored as 'fit' was provided.")
  }
  ## Grabbing values from fit if required.
  if (!is.null(fit)){
    traps <- get.traps(fit)
    mask <- get.mask(fit)
    infotypes <- fit$infotypes
    detfn <- fit$args$detfn
    pars <- get.par(fit, "fitted", as.list = TRUE)
    ss.opts <- fit$args$ss.opts
    cue.rates <- fit$args$cue.rates
    survey.length <- fit$args$survey.length
    sound.speed <- fit$args$sound.speed
    ## Setting up correct arguments for a simulating from a first-call model.
    if (!is.null(ss.opts$lower.cutoff)){
      cue.rates <- Inf
      first.only <- TRUE
    }
  }
  ## Setting up logical indicators for additional information types.
  supp.types <- c("bearing", "dist", "ss", "toa", "mrds")
  sim.types <- supp.types %in% infotypes
  names(sim.types) <- supp.types
  sim.bearings <- sim.types["bearing"]
  sim.dists <- sim.types["dist"]
  sim.toas <- sim.types["toa"]
  sim.mrds <- sim.types["mrds"]
  sim.ss <- ifelse(detfn == "ss", TRUE, FALSE)
  cutoff <- ss.opts$cutoff
  lower.cutoff <- ss.opts$lower.cutoff
  het.source <- ss.opts$het.source
  directional <- ss.opts$directional
  ss.link <- ss.opts$ss.link
  ## Sorting out inf calls stuff.
  if (any(cue.rates == Inf)){
    if (!first.only){
      warning("Setting 'first.only' to 'TRUE' as 'cue.rates' is Inf.")
      first.only <- TRUE
    }
    if (is.null(lower.cutoff)){
      #stop("A lower cutoff must be specified if individuals call until they are detected.")
    }
    inf.calls <- TRUE
  } else {
    inf.calls <- FALSE
  }
  ## Sorting out directional calling stuff.
  if (sim.ss){
    if (is.null(cutoff)){
      stop("For signal strength models, the 'cutoff' component of 'ss.opts' must be specified.")
    }
    if (is.null(directional)){
      if ("b2.ss" %in% names(pars)){
        directional <- TRUE
      } else {
        directional <- FALSE
        pars$b2.ss <- 0
      }
    } else if (directional & !("b2.ss" %in% names(pars))){
      stop("Parameter 'b2.ss' must be specified for a directional calling model.")
    } else if (!directional & "b2.ss" %in% names(pars)){
      if (pars$b2.ss != 0){
        warning("Parameter 'b2.ss' in 'pars' is being ignored as the 'directional' component of 'ss.opts' is 'FALSE'.")
        pars$b2.ss <- 0
      }
    }
    if (is.null(het.source)){
      if ("sigma.b0.ss" %in% names(pars)){
        het.source <- TRUE
      } else {
        het.source <- FALSE
        pars$sigma.b0.ss <- 0
      }
    } else if (het.source & !("sigma.b0.ss" %in% names(pars))){
      stop("Parameter 'sigma.b0.ss' must be specified for a model with heterogeneity in source strengths'.")
    } else if (!het.source & "sigma.b0.ss" %in% names(pars)){
      if (pars$sigma.b0.ss != 0){
        warning("Parameter 'sigma.b0.ss' in 'pars' is being ignores ad the 'het.source' component of 'ss.opts' is 'FALSE'.")
        pars$sigma.b0.ss <- 0
      }
    }
    if (is.null(ss.link)){
      ss.link <- "identity"
    }
    ## Setting b2.ss to 0 if model is not directional.
    if (!directional){
      pars$b2.ss <- 0
    }
    ## Setting sigma.b0.ss if model does not have heterogeneity in source strengths.
    if (!het.source){
      pars$sigma.b0.ss <- 0
    }
  }
  ## Working out required parameters.
  suppar.names <- c("kappa", "alpha", "sigma.toa")[sim.types[c("bearing", "dist", "toa")]]
  if (sim.ss){
    if (ss.link == "identity"){
      detfn <- "ss"
    } else if (ss.link == "log"){
      detfn <- "log.ss"
    } else if (ss.link == "spherical"){
      stop("Simulation for spherical spreading models is not yet implemented.")
    } else {
      stop("The argument 'ss.link' must be either \"identity\" or \"log\"")
    }
  }
  detpar.names <- switch(detfn,
                         hn = c("g0", "sigma"),
                         hr = c("g0", "sigma", "z"),
                         th = c("shape", "scale"),
                         lth = c("shape.1", "shape.2", "scale"),
                         ss = c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss"),
                         log.ss = c("b0.ss", "b1.ss", "sigma.ss"))
  par.names <- c("D", detpar.names, suppar.names)
  if (!identical(sort(par.names), sort(names(pars)))){
    msg <- paste("The following must be named components of the list 'pars': ",
                 paste(par.names, collapse = ", "), ".", sep = "")
    stop(msg)
  }
  ## Grabbing detection function parameters.
  detpars <- pars[detpar.names]
  ## Specifies the area in which animal locations can be generated.
  core <- data.frame(x = range(mask[, 1]), y = range(mask[, 2]))
  ## Simulating population.
  if (is.null(cue.rates)){
    ## Indicates which individual is being detected.
    individual <- 1:nrow(popn)
  } else {
    D <- pars$D
    if (!first.only){
      ## This is super messy, but it's scaling D from call
      ## density to animal density.
      D <- D/(mean(cue.rates)*survey.length)
    }
    if (is.null(popn)){
      popn <- as.matrix(sim.popn(D = D, core = core, buffer = 0))
    }
    n.a <- nrow(popn)
    if (freq.dist == "edf"){
      if (length(cue.rates) == 1){
        freqs <- rep(cue.rates*survey.length, n.a)
      } else {
        freqs <- sample(cue.rates*survey.length, size = n.a, replace = TRUE)
      }
    } else if (freq.dist == "norm"){
      if (diff(range(cue.rates)) == 0){
        freqs <- rep(unique(cue.rates)*survey.length, n.a)
      } else {
        freqs <- rnorm(n.a, mean(cue.rates)*survey.length, sd(cue.rates)*survey.length)
      }
    } else {
      stop("The argument 'freq.dist' must be either \"edf\" or \"norm\"")
    }
    ## Rounding frequencies up and down at random, depending
    ## on which integer is closer.
    which.integers <- floor(freqs) == freqs
    for (i in (1:n.a)[!which.integers]){
      prob <- freqs[i] - floor(freqs[i])
      freqs[i] <- floor(freqs[i]) + rbinom(1, 1, prob)
    }
    ## Indicates which individual is being detected.
    if (!first.only){
      if (n.a == 0){
        individual <- numeric(0)
      } else {
        individual <- rep(1:n.a, times = freqs)
      }
      popn <- popn[individual, , drop=FALSE]
    } else {
      individual <- seq_along(numeric(n.a))
    }
  }
  n.popn <- nrow(popn)
  ## Calculating distances.
  dists <- crossdist(popn[, 1], popn[, 2], traps[, 1], traps[, 2])
  n.traps <- nrow(traps)
  ## Creating empty bincapt if no animals in population.
  if (n.popn == 0){
    captures <- numeric(0)
    bin.capt <- matrix(0, nrow = 0, ncol = n.traps)
    out <- list(bincapt = bin.capt)
    if (sim.ss){
      out$ss <- bin.capt
    }
  } else {
    ## Calculating detection probabilities and simulating captures.
    if (!sim.ss){
      det.probs <- ascr:::calc.detfn(dists, detfn, detpars, ss.link)
      if (first.only){
        ## If only first calls are required, simulate each call separately.
        full.bin.capt <- matrix(0, nrow = n.a, ncol = n.traps)
        for (i in 1:n.a){
          det <- FALSE
          j <- 1
          while (!det & j <= freqs[i]){
            ind.bin.capt <- as.numeric(runif(n.traps) < det.probs[i, ])
            if (sum(ind.bin.capt) > 0){
              full.bin.capt[i, ] <- ind.bin.capt
              det <- TRUE
            }
            j <- j + 1
          }
        }
      } else {
        full.bin.capt <- matrix(as.numeric(runif(n.popn*n.traps) < det.probs),
                                nrow = n.popn, ncol = n.traps)
      }
      captures <- which(apply(full.bin.capt, 1, sum) > 0)
      bin.capt <- full.bin.capt[captures, , drop=FALSE]
      out <- list(bincapt = bin.capt)
    } else {
      if (ss.link == "identity"){
        inv.ss.link <- identity
      } else if (ss.link == "log"){
        inv.ss.link <- exp
      } else {
        stop("Argument 'ss.link' must be \"identity\" or \"log\".")
      }
      pars$cutoff <- cutoff
      detpars$cutoff <- cutoff
      ## Simulating animal directions and calculating orientations
      ## to traps.
      if (pars$b2.ss != 0){
        if (!is.null(cue.rates) & !first.only){
          warning("Call directions are being generated independently.")
        }
        popn.dirs <- runif(n.popn, 0, 2*pi)
        popn.bearings <- t(bearings(traps, popn))
        popn.orientations <- abs(popn.dirs - popn.bearings)
      } else {
        popn.orientations <- 0
      }
      ## Expected received strength at each microphone for each call.
      ss.mean <- inv.ss.link(pars$b0.ss - (pars$b1.ss - pars$b2.ss*(cos(popn.orientations) - 1)/2)*dists)
      ## Random error at each microphone.
      sigma.mat <- matrix(pars$sigma.b0.ss^2, nrow = n.traps, ncol = n.traps)
      diag(sigma.mat) <- diag(sigma.mat) + pars$sigma.ss^2
      if (first.only){
        if (pars$sigma.b0.ss > 0){
          stop("Simulation of first call data for situations with heterogeneity in source signal strengths is not yet implemented.")
          ## Though note that everything is OK for directional calling.
        }
        if (inf.calls){
          log.det.probs <- pnorm(cutoff, ss.mean, pars$sigma.ss, FALSE, TRUE)
          log.evade.probs <- pnorm(cutoff, ss.mean, pars$sigma.ss, TRUE, TRUE)
          ## Generating all possible capture histories.
          n.combins <- 2^n.traps
          combins <- matrix(NA, nrow = n.combins, ncol = n.traps)
          for (i in 1:n.traps){
            combins[, i] <- rep(rep(c(0, 1), each = 2^(n.traps - i)), times = 2^(i - 1))
          }
          full.bin.capt <- matrix(0, nrow = n.a, ncol = n.traps)
          for (i in 1:n.a){
            #if (sum(det.probs[i, ]) > 0){
            log.detprob.mat <- matrix(log.det.probs[i, ], nrow = n.combins, ncol = n.traps, byrow = TRUE)
            log.evadeprob.mat <- matrix(log.evade.probs[i, ], nrow = n.combins, ncol = n.traps, byrow = TRUE)
            log.prob.mat <- log.detprob.mat
            log.prob.mat[combins == 0] <- log.evadeprob.mat[combins == 0]
            ## Probabilities of each possible capture history.
            log.d.capt <- apply(log.prob.mat, 1, sum)
            d.capt <- exp(log.d.capt)
            ## Selecting a capture history.
            which.capt <- sample(2:n.combins, size = 1, prob = d.capt[2:n.combins])
            full.bin.capt[i, ] <- combins[which.capt, ]
            #}
          }
          full.ss.capt <- full.bin.capt
          full.ss.capt[full.ss.capt == 1] <- rtruncnorm(sum(full.bin.capt, na.rm = TRUE), a = cutoff,
                                                        mean = ss.mean[full.bin.capt == 1], sd = pars$sigma.ss)
        } else {
          ## If only first calls are required, simulate each call separately.
          ## Written in C++ as it was way too slow otherwise.
          full.ss.capt <- sim_ss(ss.mean, pars$sigma.ss, cutoff, freqs)
        }
      } else {
        ss.error <- rmvnorm(n.popn, sigma = sigma.mat)
        ## Filling ss.error for non-hetergeneity models for consistency with old versions.
        if (pars$sigma.b0.ss == 0){
          
          ss.error <- matrix(t(ss.error), nrow = n.popn, ncol = n.traps)
        }
        ## Creating SS capture history.
        full.ss.capt <- ss.mean + ss.error
      }
      captures <- which(apply(full.ss.capt, 1,
                              function(x, cutoff) any(x > cutoff),
                              cutoff = cutoff))
      full.bin.capt <- ifelse(full.ss.capt > cutoff, 1, 0)
      ss.capt <- full.ss.capt[captures, , drop=FALSE]
      if (length(captures) == 0){
        bin.capt <- ss.capt
      } else {
        bin.capt <- ifelse(ss.capt > cutoff, 1, 0)
      }
      ss.capt[ss.capt < cutoff] <- 0
      out <- list(bincapt = bin.capt, ss = ss.capt)
    }
  }
  ## Plot to test correct detection simulation.
  if (test.detfn){
    if (!is.null(het.source)){
      if (het.source){
        warning("Detection function testing for models with heterogeity in source strengths is not yet implemented.")
        test.detfn <- FALSE
      }
    }
  }
  if (test.detfn & n.popn != 0){
    capt.dists <- dists[full.bin.capt == 1]
    evade.dists <- dists[full.bin.capt == 0]
    all.dists <- c(capt.dists, evade.dists)
    capt.dummy <- c(rep(1, length(capt.dists)),
                    rep(0, length(evade.dists)))
    breaks <- seq(0, max(all.dists), length.out = 100)
    mids <- breaks[-length(breaks)] + 0.5*diff(breaks)
    breaks[1] <- 0
    split.dummy <- split(capt.dummy,
                         f = cut(all.dists, breaks = breaks))
    props <- sapply(split.dummy, mean)
    plot(mids, props, type = "l", xlim = c(0, max(all.dists)),
         ylim = c(0, 1))
    xx <- seq(0, max(all.dists), length.out = 100)
    lines(xx, ascr:::calc.detfn(xx, detfn, detpars, ss.link), col = "blue")
  }
  ## Total number of detections.
  n.dets <- sum(bin.capt)
  ## Keeping identities of captured individuals.
  capt.ids <- individual[captures]
  ## Locations of captured individuals.
  capt.popn <- popn[captures, ]
  ## IDs of captured individuals.
  ## Capture distances.
  capt.dists <- dists[captures, ]
  ## Simulating additional information.
  if (sim.bearings){
    bearings <- t(bearings(traps, as.matrix(capt.popn)))
    bearing.capt <- matrix(0, nrow = nrow(bin.capt),
                           ncol = ncol(bin.capt))
    bearing.capt[bin.capt == 1] <- (bearings[bin.capt == 1] +
                                      rvm(n.dets, mean = 0, k = pars$kappa)) %% (2*pi)
    out$bearing <- bearing.capt
  }
  if (sim.dists){
    dist.capt <- matrix(0, nrow = nrow(bin.capt),
                        ncol = ncol(bin.capt))
    betas <- pars$alpha/capt.dists[bin.capt == 1]
    dist.capt[bin.capt == 1] <- rgamma(n.dets, shape = pars$alpha,
                                       rate = betas)
    out$dist <- dist.capt
  }
  if (sim.toas){
    ## Time taken for sound to travel from source to detector.
    toa.capt <- capt.dists/sound.speed*bin.capt
    ## Adding in TOA error.
    toa.capt[bin.capt == 1] <-
      toa.capt[bin.capt == 1] + rnorm(n.dets, sd = pars$sigma.toa)
    out$toa <- toa.capt
  }
  if (sim.mrds){
    out$mrds <- capt.dists
  }
  if (keep.locs | keep.ids){
    out <- list(capt = out)
    if (keep.locs){
      out[["capt.locs"]] <- capt.popn
      out[["popn.locs"]] <- popn
    }
    if (keep.ids){
      out[["capt.ids"]] <- capt.ids
    }
  }
  out
}
