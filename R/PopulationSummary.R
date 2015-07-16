
#' Take in the results from the population data and re-format it
#' @return list, data frame of logical (passed or not), data frame of numeric (all results)
#' @param results The list of results from running the package using BatchJobs
#' @keywords private
reformatData <- function(results) {
    all.samples.results <-
        as.data.frame(do.call("rbind", lapply(results, "[[", 1)))
    all.samples.names <-
        as.data.frame(do.call("rbind", lapply(results, "[[", 3)))
    rownames(all.samples.results) <- all.samples.names[,1]
    all.samples.pass <-
        as.data.frame(do.call("rbind", lapply(results, "[[", 2)))
    rownames(all.samples.pass) <- all.samples.names[,1]
    return(list("pass"=all.samples.pass, "results"=all.samples.results))
   }

