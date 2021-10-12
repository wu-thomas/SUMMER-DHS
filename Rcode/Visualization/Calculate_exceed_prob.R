## By Richard

##########################################################
#  Function: to be included in SUMMER
##########################################################
#' @param est1 result from getSmoothed from cluster-level model with save.draws = TRUE.
#' @param est2 can be one of the following: 
#'  \itemize{ 
#' \item a numeric number
#' \item output from getDirect.
#' \item output from getSmoothed from cluster-level model with save.draws = TRUE.
#' }
#' For the latter two types of input, \code{grouping} needs to be specified.
#' @param grouping when comparing two set of estimates, grouping should contain the following columns: years.1, region.1, strata.1, years.2, region.2, strata.2. If no strata exist in the estimates, use 'strata_all' for all rows in the strata columns. Each row of this data frame correspond to one comparison of the following: Pr(estimate of years.1 in region.1 from the first model, i.e., est1, > estimate of years.2 in region.2 from the second model, i.e., est2)
#' 
exceedProb <- function(est1, est2, grouping = NULL,strat=TRUE){
  
    if(strat==TRUE){
    draw.list <- est1$draws.est.overall

        for(i in 1:length(draw.list)){
          draw.list[[i]]$strata <- 'strata_all'
        }
    
    
    }else{    draw.list <- est1$draws.est}
    
    if(is.null(draw.list)){
        stop("Posterior draws of the estimates are not saved. Please rerun getSmoothed function with 'save.draws = TRUE'")
    }
    out <- NULL
    # if another model is used as comparison
    if(!is.numeric(est2)){
        order.list <- NULL
        if(is.null(grouping)){
            stop("grouping needs to be specified comparing two set of estimates.")
        }
        if(sum(!c("years.1", "region.1", "years.2", "region.2") %in% colnames(grouping)) > 0){
            stop("grouping not correctly specified with the required column names. See ?exceedProb for details")
        }
        for(i in 1:length(draw.list)){
            years <- draw.list[[i]]$years
            region <- draw.list[[i]]$region
            strata <- draw.list[[i]]$strata
            which <- which(grouping$years.1 == years &
                           grouping$region.1 == region & 
                           grouping$strata.1 == strata)
            if(length(which) > 1){
                stop(paste0("There exist multiple rows corresponding to years = ", years, ", region = ", region, ", strata = ", strata))
            }
            if(length(which) == 0){
                which <- NA
            }
            order.list <- c(order.list, which)
        }
    }
    est2.tab <- NULL
    if(!is.data.frame(est2) && !is.numeric(est2)){
        for(i in 1:length(est2$draws.est)){
            est2.tab <- rbind(est2.tab, data.frame(
                years =  est2$draws.est[[i]]$years,
                region = est2$draws.est[[i]]$region,
                strata = est2$draws.est[[i]]$strata)) 
        }
    }else if(is.data.frame(est2)){
         est2.tab <- data.frame(years = est2$years, region = est2$region)
         if("strata" %in% colnames(est2)){
            est2.tab$strata <- est2$strata
         }else{
            est2.tab$strata <- "strata_all"
         }
    }
    for(i in 1:length(draw.list)){
        years <- draw.list[[i]]$years
        region <- draw.list[[i]]$region
        strata <- draw.list[[i]]$strata
        draws <- draw.list[[i]]$draws

        if(is.numeric(est2)){
            out <- rbind(out, data.frame(years = years, region = region, strata = strata, prob = sum(draws > est2) / length(draws)))
        }else{
            # find the correspondence in the list
            j <- which(
                est2.tab$years == grouping$years.2[order.list[i]] & 
                est2.tab$region == grouping$region.2[order.list[i]] &
                est2.tab$strata == grouping$strata.2[order.list[i]])
           if(length(j) > 1){
                stop(paste0("More than one rows of direct estimates corresponding to years = ", grouping$years.2[order.list[i]], ", region = ", grouping$region.2[order.list[i]], ", strata = ", grouping$strata.2[order.list[i]]))
                
            }
            # direct estimates
            if(is.data.frame(est2)){
                if(length(j) == 0){
                    draws2 <- NA
                }else{   
                    draws2 <- expit(rnorm(length(draws), mean = est2$logit.est[j], sd = sqrt(est2$var.est[j])))
                }
                out <- rbind(out, data.frame(years = years, region = region, strata = strata, prob = sum(draws > draws2) / length(draws)))

            # posterior draws
            }else{ 
                if(length(j) == 0){
                    draws2 <- NA
                }else{            
                    draws2 <- est2$draws.est[[j]]$draws
                }
                out <- rbind(out, data.frame(years = years, region = region, strata = strata, prob = sum(draws > draws2) / length(draws)))
            }
        }
    }
    return(out)
}
