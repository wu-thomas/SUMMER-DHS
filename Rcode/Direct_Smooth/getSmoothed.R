### Code from Jessica

getSmoothed_sd <- function (inla_mod, nsim = 1000, weight.strata = NULL, weight.frame = NULL, 
                         verbose = FALSE, mc = 0, include_time_unstruct = FALSE, CI = 0.95, 
                         draws = NULL, save.draws = FALSE, include_subnational = TRUE, 
                         ...) 
{
  years <- region <- NA
  lowerCI <- (1 - CI)/2
  upperCI <- 1 - lowerCI
  save.draws.est <- save.draws
  if (!is.null(inla_mod$year_range)) {
    year_range <- inla_mod$year_range
  }
  else {
    warning("The fitted object was from an old version of SUMMER, please specify 'year_range' argument when calling getSmoothed()")
  }
  if (!is.null(inla_mod$year_label)) {
    year_label <- inla_mod$year_label
  }
  else {
    warning("The fitted object was from an old version of SUMMER, please specify 'year_label' argument when calling getSmoothed()")
  }
  if (!is.null(inla_mod$has.Amat)) {
    Amat <- inla_mod$Amat
  }
  else {
    warning("The fitted object was from an old version of SUMMER, please specify 'Amat' argument when calling getSmoothed()")
  }
  if (!is.null(inla_mod$family)) {
    if ("region.struct" %in% names(inla_mod$fit$summary.random) == 
        FALSE && !is.null(Amat)) {
      warning("No spatial random effects in the model. Set Amat to NULL", 
              immediate. = TRUE)
      Amat <- NULL
    }
    is.dynamic <- as.logical(inla_mod$strata.time.effect)
    if (length(is.dynamic) == 0) 
      is.dynamic = FALSE
    if (!is.dynamic) {
      stratalabels <- stratalabels.orig <- inla_mod$strata.base
      if (length(stratalabels) == 1 && stratalabels[1] == 
          "") {
        stratalabels <- "strata_all"
        stratalabels.orig <- "strata_all"
      }
      other <- rownames(inla_mod$fit$summary.fixed)
      other <- other[grep("strata", other)]
      other <- gsub("strata", "", other)
      stratalabels <- c(stratalabels, other)
      stratalabels.orig <- c(stratalabels.orig, other)
      frame.strata <- (inla_mod$fit$.args$data)[, c("age", 
                                                    "age.idx", "age.rep.idx")]
      frame.strata <- unique(frame.strata)
      frame.strata <- frame.strata[order(frame.strata$age.idx), 
                                   ]
      framelabels <- "frame_all"
      multi.frame <- FALSE
    }
    else {
      if ("frame" %in% colnames(inla_mod$fit$.args$data)) {
        frame.strata <- (inla_mod$fit$.args$data)[, c("age", 
                                                      "age.idx", "age.rep.idx", "strata", "frame")]
        frame.strata <- unique(frame.strata)
        frame.strata <- frame.strata[!is.na(frame.strata$frame), 
                                     ]
        frame.strata$strata.orig <- NA
        for (i in 1:dim(frame.strata)[1]) {
          frame.strata$strata.orig[i] <- gsub(paste0(frame.strata$frame[i], 
                                                     "-"), "", frame.strata$strata[i])
        }
      }
      else {
        frame.strata <- (inla_mod$fit$.args$data)[, c("age", 
                                                      "age.idx", "strata", "age.rep.idx")]
        frame.strata <- unique(frame.strata)
        frame.strata$strata.orig <- frame.strata$strata
      }
      frame.strata <- frame.strata[order(frame.strata$age.idx), 
                                   ]
      stratalabels <- as.character(unique(frame.strata$strata))
      stratalabels.orig <- as.character(unique(frame.strata$strata.orig))
      framelabels <- as.character(unique(frame.strata$frame))
      if (length(framelabels) == 0) {
        framelabels <- "frame_all"
        multi.frame <- FALSE
      }
      else if (length(framelabels) == 1) {
        multi.frame <- FALSE
      }
      else if (sum(frame.strata$frame != frame.strata$strata) == 
               0) {
        multi.frame <- FALSE
      }
      else {
        multi.frame <- TRUE
      }
    }
    if (inla_mod$is.yearly) 
      year_label <- c(year_range[1]:year_range[2], year_label)
    if (!inla_mod$is.yearly) 
      year_label <- year_label
    err <- NULL
    weight.strata.by <- NULL
    if (!is.null(weight.strata)) {
      if ((!"frame" %in% colnames(inla_mod$fit$.args$data)) && 
          "frame" %in% colnames(weight.strata)) {
        stop("frame variable is not in the fitted model, but exists in the weight.strata data frame. Please remove the column from weight.strata and use a single set of weights (that are not specific to sampling frames).")
      }
      if (!is.dynamic && sum(stratalabels %in% colnames(weight.strata)) != 
          length(stratalabels)) {
        stop(paste0("weight.strata argument not specified correctly. It requires the following columns: ", 
                    paste(stratalabels, collapse = ", ")))
      }
      if (is.dynamic && "frame" %in% colnames(frame.strata) && 
          !("frame" %in% colnames(weight.strata))) {
        stop("frame column was provided in the fitted object but not the weights.")
      }
      if (sum(c("region", "years") %in% colnames(weight.strata)) == 
          2) {
        for (i in year_label) {
          tmp <- colnames(Amat)[which(colnames(Amat) %in% 
                                        subset(weight.strata, years == i)$region == 
                                        FALSE)]
          if (length(tmp) > 0) 
            err <- c(err, paste(tmp, i))
        }
        if (!is.null(err)) {
          stop(paste0("The following region-year combinations are not present in the strata weights: ", 
                      paste(err, collapse = ", ")))
        }
        weight.strata.by <- c("region", "years")
      }
      else if ("region" %in% colnames(weight.strata) && 
               dim(weight.strata)[1] > 1) {
        warning("Time is not specified, assuming the strata weights are static.", 
                immediate. = TRUE)
        tmp <- colnames(Amat)[which(colnames(Amat) %in% 
                                      weight.strata$region == FALSE)]
        if (length(tmp) > 0) {
          stop(paste0("The following regions are not present in the strata weights: ", 
                      paste(tmp, collapse = ", ")))
        }
        weight.strata.by <- c("region")
      }
      else if ("years" %in% colnames(weight.strata)) {
        tmp <- year_label[which(year_label %in% weight.strata$years == 
                                  FALSE)]
        if (length(tmp) > 0) {
          stop(paste0("The following time periods are not present in the strata weights: ", 
                      paste(tmp, collapse = ", ")))
        }
        weight.strata.by <- c("years")
      }
      else {
        warning("no region or years in the strata proportion. Treat proportion as constant.", 
                immediate. = TRUE)
        weight.strata.by = "Constant"
      }
    }
    cs <- inla_mod$fit$misc$configs$contents$tag
    cs <- cs[cs != "Predictor"]
    cs <- cs[cs != "nugget.id"]
    select <- list()
    for (i in 1:length(cs)) {
      select[[i]] <- 0
      names(select)[i] <- cs[i]
    }
    if (is.null(draws)) {
      message("Starting posterior sampling...")
      sampAll <- INLA::inla.posterior.sample(n = nsim, 
                                             result = inla_mod$fit, intern = TRUE, selection = select, 
                                             verbose = verbose)
      message("Cleaning up results...")
    }
    else {
      message("Use posterior draws from input.")
      sampAll <- draws
      nsim <- length(draws)
    }
    fields <- rownames(sampAll[[1]]$latent)
    T <- max(inla_mod$fit$.args$data$time.struct)
    if (!"region.struct:1" %in% fields) 
      inla_mod$fit$.args$data$region.struct = 1
    rep.time <- length(unique(inla_mod$age.rw.group)) > 0
    cols <- c("time.struct", "time.unstruct", "region.struct", 
              "time.area", "strata", "age", "age.idx")
    if (rep.time) {
      cols <- c(cols, "age.rep.idx")
    }
    A <- unique(inla_mod$fit$.args$data[, cols])
    if (sum(stratalabels == "strata_all") == length(stratalabels)) {
      A$strata <- "strata_all"
    }
    if (inla_mod$strata.time.effect) {
      AA <- expand.grid(time.area = unique(A$time.area), 
                        age = unique(A$age))
      AA <- merge(AA, unique(A[, c("age", "strata")]), 
                  by = "age")
    }
    else {
      AA <- expand.grid(time.area = unique(A$time.area), 
                        strata = unique(A$strata), age = unique(A$age))
      if (sum(AA$strata != "") == 0) 
        AA$strata <- "strata_all"
    }
    AA <- merge(AA, unique(A[, c("time.struct", "time.unstruct", 
                                 "region.struct", "time.area")]), by = "time.area")
    if (rep.time) {
      AA <- merge(AA, unique(A[, c("age", "age.idx", "age.rep.idx")]), 
                  by = "age")
    }
    else {
      AA <- merge(AA, unique(A[, c("age", "age.idx")]), 
                  by = "age")
    }
    if (rep.time) {
      AA$time.struct <- AA$time.struct + (AA$age.rep.idx - 
                                            1) * T
    }
    AA$age <- paste0("age", AA$age, ":1")
    if (length(unique(AA$age)) == 1) 
      AA$age <- "(Intercept):1"
    AA.loc <- AA
    AA.loc$age <- match(AA.loc$age, fields)
    if (!is.dynamic) 
      AA.loc$strata <- paste0("strata", AA.loc$strata, 
                              ":1")
    AA.loc$strata <- match(AA.loc$strata, fields)
    AA.loc$time.area <- match(paste0("time.area:", AA.loc$time.area), 
                              fields)
    if ("region.int:1" %in% rownames(sampAll[[1]]$latent)) {
      AA.loc$time.area <- (AA.loc$time.unstruct - 1) * 
        dim(Amat)[1] + AA.loc$region.struct
      AA.loc$time.area <- match(paste0("region.int:", AA.loc$time.area), 
                                fields)
    }
    AA.loc$time.struct <- match(paste0("time.struct:", AA.loc$time.struct), 
                                fields)
    AA.loc$region.struct <- match(paste0("region.struct:", 
                                         AA.loc$region.struct), fields)
    AA.loc$time.unstruct <- match(paste0("time.unstruct:", 
                                         AA.loc$time.unstruct), fields)
    if (is.logical(include_time_unstruct)) {
      if (!include_time_unstruct) 
        AA.loc$time.unstruct <- NA
    }
    else {
      which.include <- which(year_label %in% as.character(include_time_unstruct))
      included <- which(AA$time.unstruct %in% which.include)
      not_included <- which(AA$time.unstruct %in% which.include == 
                              FALSE)
      AA.loc$time.unstruct[not_included] <- NA
      message(paste0("The IID temporal components are included in the following time periods: ", 
                     paste(year_label[which.include], collapse = ", ")))
    }
    slope <- grep("time.slope.group", fields)
    if (!is.null(slope)) {
      AA$tstar <- (AA$time.unstruct - (T + 1)/2)/(T + 1)
      AA$slope <- match(paste0("time.slope.group", AA$age.rep.idx, 
                               ":1"), fields)
    }
    else {
      AA$tstar <- AA$slope <- NA
    }
    st.slope <- grep("st.slope.id", fields)
    if (!is.null(st.slope)) {
      AA$ststar <- (AA$time.unstruct - (T + 1)/2)/(T + 
                                                     1)
      AA$st.slope <- match(paste0("st.slope.id:", AA$region.struct), 
                           fields)
    }
    else {
      AA$ststar <- AA$st.slope <- NA
    }
    AA.loc$age.idx <- AA.loc$age.rep.idx <- NA
    if (!include_subnational) {
      AA.loc$time.area <- NA
      AA.loc$region.struct <- NA
    }
    age <- match(paste0("age", frame.strata$age, ":1"), fields)
    age.nn <- inla_mod$age.n
    if (!is.null(Amat)) {
      N <- dim(Amat)[1]
    }
    else {
      N <- 1
    }
    if (length(age) == 0) {
      age.length <- 1
    }
    else {
      age.length <- length(age)
    }
    out1 <- expand.grid(strata = stratalabels, time = 1:T, 
                        area = 1:N)
    out2 <- expand.grid(frame = framelabels, time = 1:T, 
                        area = 1:N)
    out3 <- expand.grid(time = 1:T, area = 1:N)
    out1$lower <- out1$upper <- out1$mean <- out1$median <- out1$variance <- NA
    out2$lower <- out2$upper <- out2$mean <- out2$median <- out2$variance <- NA
    out3$lower <- out3$upper <- out3$mean <- out3$median <- out3$variance <- NA
    out1$years <- year_label[out1$time]
    out2$years <- year_label[out2$time]
    out3$years <- year_label[out3$time]
    if (N > 1) {
      out1$region <- colnames(Amat)[out1$area]
      out2$region <- colnames(Amat)[out2$area]
      out3$region <- colnames(Amat)[out3$area]
    }
    else {
      out1$region <- out2$region <- out3$region <- "All"
    }
    if (is.null(weight.strata)) {
      if (length(stratalabels) > 1) {
        message("No strata weights has been supplied. Set all weights to 0.")
      }
      else {
        message("No stratification in the model. Set all weights to 1.")
      }
      if (!is.null(Amat)) {
        weight.strata <- expand.grid(region = colnames(Amat), 
                                     frame = framelabels)
        weight.strata.by <- "region"
      }
      else {
        weight.strata <- data.frame(frame = framelabels)
        weight.strata.by <- "Constant"
      }
      for (tt in stratalabels.orig) {
        weight.strata[, tt] <- ifelse(length(stratalabels) > 
                                        1, 0, 1)
      }
    }
    if (weight.strata.by[1] == "Constant") {
      if ("frame" %in% colnames(weight.strata)) {
        out2 <- merge(out2, weight.strata, by = "frame")
      }
      else {
        out2 <- cbind(out2, data.frame(weight.strata))
      }
      strata.index <- match(stratalabels.orig, colnames(out2))
    }
    else {
      if (length(framelabels) == 1) {
        out2 <- merge(out2[, colnames(out2) != "frame"], 
                      weight.strata, by = weight.strata.by)
      }
      else {
        weight.strata.by <- c(weight.strata.by, "frame")
        out2 <- merge(out2, weight.strata, by = weight.strata.by)
      }
      strata.index <- match(stratalabels.orig, colnames(out2))
      out2 <- out2[with(out2, order(area, time)), ]
    }
    tau <- rep(NA, nsim)
    theta <- matrix(0, nsim, dim(AA)[1])
    for (i in 1:nsim) {
      draw <- sampAll[[i]]$latent
      theta[i, ] <- apply(AA.loc, 1, function(x, ff) {
        sum(ff[x], na.rm = TRUE)
      }, draw)
      add.slope <- draw[AA$slope] * AA$tstar
      add.slope[is.na(add.slope)] <- 0
      add.slope.st <- draw[AA$st.slope] * AA$ststar
      add.slope.st[is.na(add.slope.st)] <- 0
      theta[i, ] <- theta[i, ] + add.slope + add.slope.st
      if (inla_mod$family == "binomial") {
        tau[i] <- exp(sampAll[[i]]$hyperpar[["Log precision for nugget.id"]])
      }
    }
    if (inla_mod$family == "binomial" && (mc == 0)) {
      k <- 16 * sqrt(3)/15/base::pi
      theta <- theta/sqrt(1 + k^2/tau)
    }
    draw.temp <- draws.hazards <- NA
    draws.est <- NULL
    index.draws.est <- 1
    draws.est.overall <- NULL
    index.draws.est.overall <- 1
    index1 <- 1
    index2 <- 1
    if (N == 1 && AA$region.struct[1] == 0) 
      AA$region.struct <- 1
    for (j in 1:N) {
      if (inla_mod$family == "binomial" && mc > 0) {
        sd.temp <- matrix(1/sqrt(tau), nsim, mc)
        err.temp <- matrix(stats::rnorm(nsim * mc, mean = matrix(0, 
                                                                 nsim, mc), sd = sd.temp), nsim, mc)
      }
      for (i in 1:T) {
        sub <- which(AA$time.unstruct == i & AA$region.struct == 
                       j)
        AA.sub <- AA[sub, ]
        draws.sub <- theta[, sub, drop = FALSE]
        draws.sub.agg <- matrix(NA, nsim, length(stratalabels))
        for (k in 1:length(stratalabels)) {
          strata.sub <- which(AA.sub$strata == stratalabels[k])
          draws.hazards <- draws.sub[, strata.sub, drop = FALSE]
          if (inla_mod$family == "binomial" && mc > 0) {
            for (tt in 1:dim(draws.hazards)[2]) {
              draws.temp <- matrix(draws.hazards[, tt], 
                                   nsim, mc)
              draws.temp <- expit(draws.temp + err.temp)
              draws.hazards[, tt] <- apply(draws.temp, 
                                           1, mean)
            }
          }
          else {
            draws.hazards <- expit(draws.hazards)
          }
          draws.mort <- rep(1, dim(draws.hazards)[1])
          for (tt in 1:dim(draws.hazards)[2]) {
            draws.mort <- draws.mort * (1 - draws.hazards[, 
                                                          tt])^age.nn[AA.sub[strata.sub, "age.idx"][tt]]
          }
          draws.mort <- 1 - draws.mort
          draws.sub.agg[, k] <- draws.mort
          index1 <- which(out1$time == i & out1$area == 
                            j & out1$strata == stratalabels[k])
          if (save.draws.est) {
            draws.est[[index.draws.est]] <- list(years = year_label[i], 
                                                 region = colnames(Amat)[j], strata = stratalabels[k], 
                                                 draws = draws.mort)
            index.draws.est <- index.draws.est + 1
          }
          out1[index1, c("lower", "median", "upper")] <- quantile(draws.mort, 
                                                                  c(lowerCI, 0.5, upperCI))
          out1[index1, "mean"] <- mean(draws.mort)
          out1[index1, "variance"] <- var(draws.mort)
        }
        index2 <- which(out2$area == j & out2$time == 
                          i)
        draws.sub.agg.sum <- matrix(NA, nsim, length(index2))
        for (k in 1:length(index2)) {
          prop <- out2[index2[k], strata.index]
          if (!multi.frame) {
            cols <- match(stratalabels.orig, stratalabels)
          }
          else {
            cols <- match(paste(out2$frame[index2[k]], 
                                stratalabels.orig, sep = "-"), stratalabels)
          }
          draws.sub.agg.sum[, k] <- apply(draws.sub.agg[, 
                                                        cols, drop = FALSE], 1, function(x, p) {
                                                          sum(x * p)
                                                        }, prop)
          if (save.draws.est) {
            draws.est.overall[[index.draws.est.overall]] <- list(years = year_label[i], 
                                                                 region = colnames(Amat)[j], draws = draws.sub.agg.sum[, 
                                                                                                                       k])
            index.draws.est.overall <- index.draws.est.overall + 
              1
          }
          out2[index2[k], c("lower", "median", "upper")] <- quantile(draws.sub.agg.sum[, 
                                                                                       k], c(lowerCI, 0.5, upperCI))
          out2[index2[k], "mean"] <- mean(draws.sub.agg.sum[, 
                                                            k])
          out2[index2[k], "variance"] <- var(draws.sub.agg.sum[, 
                                                               k])
        }
        if (!is.null(weight.frame)) {
          index3 <- which(out3$area == j & out3$time == 
                            i)
          colnames(draws.sub.agg.sum) <- out2[index2, 
                                              "frame"]
          draws.sub.agg.sum2 <- rep(0, nsim)
          this.weight <- weight.frame
          if ("region" %in% colnames(weight.frame)) 
            this.weight <- subset(this.weight, region == 
                                    colnames(Amat)[j])
          if ("years" %in% colnames(weight.frame)) 
            this.weight <- subset(this.weight, years == 
                                    year_label[i])
          for (k in 1:dim(draws.sub.agg.sum)[2]) {
            draws.sub.agg.sum2 <- draws.sub.agg.sum2 + 
              logit(draws.sub.agg.sum[, k]) * as.numeric(this.weight[colnames(draws.sub.agg.sum)[k]])
          }
          draws.sub.agg.sum2 <- expit(draws.sub.agg.sum2)
          out3[index3, c("lower", "median", "upper")] <- quantile(draws.sub.agg.sum2, 
                                                                  c(lowerCI, 0.5, upperCI))
          out3[index3, "mean"] <- mean(draws.sub.agg.sum2)
          out3[index3, "variance"] <- var(draws.sub.agg.sum2)
        }
      }
    }
    out1$is.yearly <- !(out1$years %in% year_label)
    out1$years.num <- suppressWarnings(as.numeric(as.character(out1$years)))
    out2$is.yearly <- !(out2$years %in% year_label)
    out2$years.num <- suppressWarnings(as.numeric(as.character(out2$years)))
    out3$is.yearly <- !(out3$years %in% year_label)
    out3$years.num <- suppressWarnings(as.numeric(as.character(out3$years)))
    out1$years <- factor(out1$years, year_label)
    out2$years <- factor(out2$years, year_label)
    out3$years <- factor(out3$years, year_label)
    class(out1) <- c("SUMMERproj", "data.frame")
    class(out2) <- c("SUMMERproj", "data.frame")
    class(out3) <- c("SUMMERproj", "data.frame")
    out <- list(overall = out2, stratified = out1)
    if (!is.null(weight.frame)) {
      out$final = out3
    }
    if (save.draws) {
      out$draws = sampAll
    }
    if (save.draws.est) {
      out$draws.est <- draws.est
      out$draws.est.overall <- draws.est.overall
    }
    return(out)
  }
  else {
    if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
      stop("You need to install the packages 'INLA'. Please run in your R terminal:\n  install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
    }
    if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
      if (!is.element("INLA", (.packages()))) {
        attachNamespace("INLA")
      }
      if (is.null(Amat)) {
        region_names <- "All"
        region_nums <- 0
      }
      else {
        region_names <- colnames(Amat)
        region_nums <- 1:length(region_names)
      }
      is.yearly = inla_mod$is.yearly
      if (is.yearly) {
        timelabel.yearly <- c(year_range[1]:year_range[2], 
                              year_label)
      }
      else {
        timelabel.yearly <- year_label
      }
      results <- expand.grid(District = region_nums, Year = timelabel.yearly)
      results$median <- results$lower <- results$upper <- results$logit.median <- results$logit.lower <- results$logit.upper <- NA
      mod <- inla_mod$fit
      lincombs.info <- inla_mod$lincombs.info
      if(!is.null(draws)){
        warning("draws is an argument for inla_mod objects output from smoothCluster.")
      }

      if(save.draws.est){
        draws.est <- matrix(NA, 
                            nrow = length(region_names)*length(timelabel.yearly),
                            ncol = nsim)
      }
      for (i in 1:length(timelabel.yearly)) {
        for (j in 1:length(region_names)) {
          index <- lincombs.info$Index[lincombs.info$District == 
                                         region_nums[j] & lincombs.info$Year == i]
          tmp.logit <- INLA::inla.rmarginal(nsim, mod$marginals.lincomb.derived[[index]])
          tmp <- expit(tmp.logit)
          if(save.draws.est){
            draws.est[index,] <- tmp
          }
          results$median[results$District == region_nums[j] & 
                           results$Year == timelabel.yearly[i]] <- stats::median(tmp)
          results$upper[results$District == region_nums[j] & 
                          results$Year == timelabel.yearly[i]] <- stats::quantile(tmp, 
                                                                                  upperCI)
          results$lower[results$District == region_nums[j] & 
                          results$Year == timelabel.yearly[i]] <- stats::quantile(tmp, 
                                                                                  lowerCI)
          results$logit.median[results$District == region_nums[j] & 
                                 results$Year == timelabel.yearly[i]] <- stats::median(tmp.logit)
          results$logit.upper[results$District == region_nums[j] & 
                                results$Year == timelabel.yearly[i]] <- stats::quantile(tmp.logit, 
                                                                                        upperCI)
          results$logit.lower[results$District == region_nums[j] & 
                                results$Year == timelabel.yearly[i]] <- stats::quantile(tmp.logit, 
                                                                                        lowerCI)
        }
      }
      results$is.yearly <- !(results$Year %in% year_label)
      results$years.num <- suppressWarnings(as.numeric(as.character(results$Year)))
      if (region_names[1] != "All") {
        results$District <- region_names[results$District]
      }
      else {
        results$District <- "All"
      }
      colnames(results)[which(colnames(results) == "District")] <- "region"
      colnames(results)[which(colnames(results) == "Year")] <- "years"
      class(results) <- c("SUMMERproj", "data.frame")
      if(save.draws.est){
        out <- list(results = results,
                    draws.est = draws.est)
        return(out)
      }else{
        return(results)
      }
    }
  }
}
