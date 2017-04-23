
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


#' Using a variable presumed to be associated with a batch effect, create a histogram to visualize the differences
#' @return ggplot2
#' @export
#' @param results a data frame with population level results
#' @param batcheffect meta data variable presumed to be associated with a batch effect
#' @param variable variable of interest from genotypeeval results
#' @param mycolors vector of colors for histogram
#' @examples
#' results <- read.table(system.file("ext-data", "example_results.txt", package="genotypeeval"), header=TRUE)
#' batcheffect <- "batch"
#' variable <- "medianReadDepth"
#' makeHistogram(results, batcheffect, variable)

makeHistogram <- function(results, batcheffect, variable, mycolors=c("red", "blue")) {
  dat.s <- data.frame("batch.var"=get(batcheffect, results), "var"=get(variable, results))
g <- ggplot(dat.s, aes(fill=as.factor(dat.s$batch.var), x=dat.s$var)) + geom_histogram() + 
  scale_x_continuous(variable)+ theme_minimal()+theme(legend.position="none", text=element_text(size=30))+ scale_y_continuous("Count")+ scale_fill_manual(values=mycolors) 
return(g)
}

#' Wrapper to create a tSNE plot from Rtsne.  See Rtsne for further documentation.
#' @return ggplot2
#' @export
#' @param dat a data frame with variables of interest for dimension reduction and meta variable for batch effect
#' @param batcheffect meta data variable presumed to be associated with a batch effect to color tsne plot
#' @param dims dims for tSNE plot, tuning parameter
#' @param initial_dims initial_dims for tSNE plot, tuning parameter
#' @param perplexity perplexity for tSNE plot, tuning parameter
#' @param max.iter max.iter for the tSNE plot, tuning parameter
#' @examples
#' results <- read.table(system.file("ext-data", "example_results.txt", package="genotypeeval"), header=TRUE)
#' results$sampleid <- NULL
#' g <- tsnePlot(results, "batch")
#' plot(g)

tsnePlot <- function(dat, batcheffect, dims=2, initial_dims=50, perplexity = 30, max.iter=5000) {
  batches <- get(batcheffect, dat)
  mynum <- which(names(dat) == batcheffect)
  X <- dat[,-mynum]
  tsne_out <- Rtsne(X, perplexity=perplexity, initial_dims=initial_dims)
  plotMe <- as.data.frame(tsne_out$Y)
  names(plotMe) <- c("tsne.1", "tsne.2")
  plotMe$batch <- batches
  g <- ggplot(data=plotMe, aes(x=plotMe$tsne.1, y=plotMe$tsne.2, color=as.factor(plotMe$batch))) + geom_point() + scale_color_manual(values=c("blue", "black")) + theme(text=element_text(size=15))
}






#' Wrapper to create a PCA plot, centers and scales variables by default
#' @return ggplot2
#' @export
#' @param dat a data frame with variables of interest for dimension reduction and meta variable for batch effect
#' @param batcheffect meta data variable presumed to be associated with a batch effect to color tsne plot
#' @examples
#' results <- read.table(system.file("ext-data", "example_results.txt", package="genotypeeval"), header=TRUE)
#' results$sampleid <- NULL
#' g <- pcaPlot(results, "batch")

pcaPlot <- function(dat, batcheffect) {
  batches <- as.factor(get(batcheffect, dat))
  mynum <- which(names(dat) == batcheffect)
  X <- dat[,-mynum]
  pcs <- prcomp(X, center=TRUE, scale=TRUE)
  plotMe <- pcs$x
  #plot(plotMe)
  pairs(plotMe, col=batches) 
}












