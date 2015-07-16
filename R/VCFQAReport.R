#Class to Read in VCF File

#' Declare class VCFQAReport which will evaluate a VCF stored as a ReadData object.
#' 
#' @param printnames List of tests applied to VCF
#' @param results Numeric vector of metrics calculated from VCF file
#' @param plots List of plots created from VCF File
#' @param tests TRUE (passed) or FALSE (failed) logical vector of whether VCF passed metrics using thresholds from VCFQAParam
#' @param fn Filename of VCF evaluated (for plot titles)
#############################################################################################
#############################################################################################
#Constructors
#############################################################################################

setClass(Class="VCFQAReport", representation=representation(printnames="character", results="numeric", plots="list", tests = "logical", fn="character")
)




#' Constructor for class.  Calls constructor for class.  Using the GENO fields present in the vcf header will evaluate the vcf file using metrics and generate plots.  Each metric will be tested against the params specified in the params class.  For example, if Read Depth is in the GENO header will calculate median read depth, percent in target (50 percent to 200 percent of the target specified in the params file) and generate a histogram of Read Depth.
#' @return Object of VCFQAReport.  
#' @export
#' @param myvcf Vcf file to evaluate
#' @param vcfparams object of VCFQAParam class.  Sets thresholds to evaluate the VCF File against.
#' @param gold.ref Object of class Gold that contains the 1000 Genomes reference
#' @param cds.ref Coding Region as GRanges 
#' @param masked.ref optional regions as GRanges to mask eg repeats, self chain, paralogs, etc.
#' @param admixture.ref VRanges with MAF for superpopulations (EAS, AFR, EUR)
#' @examples
#' vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
#' mydir <- paste(dirname(vcffn), "/", sep="")
#' myfile <-basename(vcffn)
#' svp <- ScanVcfParam(which=GRanges("22", IRanges(0,200e5)), geno="GT")
#' vcfparams <- VCFQAParam(count.limits=c(3014580000, Inf), readdepth.target = 30)
#' vcf <- ReadVCFData(mydir, myfile, "GRCh38")
#' ev <- VCFEvaluate(vcf, vcfparams)

VCFEvaluate <- function(myvcf, vcfparams, gold.ref=NA, cds.ref=NA, masked.ref =NA, admixture.ref=NA) {
    cat( " Evaluating File ...  \n")
    .Object <- new(Class="VCFQAReport", fn=myvcf@myfile)
              results <- NULL
    #do the admixture estimation prior to the GQ filter
             if (length(admixture.ref) > 1) {
                  results$admixture <- admixture(myvcf, admixture.ref)
              }
    #calc this prior to filters
                  if ("GT" %in% myvcf@genoString) {
                     results$numberOfHomRefs <- numberOfHomRefs(myvcf, vcfparams@homref.limits[1], vcfparams@homref.limits[2])
                      results$numberCalls <- numberCalls(myvcf, results$numberOfHomRefs[[3]], myvcf@n.dup, vcfparams@count.limits[1], vcfparams@count.limits[2])
                 }
              #apply the GQ filter
              if (vcfparams@gq.filter > 0) {
                 myvcf@vr <- subset(myvcf@vr, myvcf@vr$GQ > vcfparams@gq.filter)
             }
              #apply the DP filter
    if (vcfparams@dp.filter > -1) {
                  myvcf@vr <- subset(myvcf@vr, totalDepth(myvcf@vr) > vcfparams@dp.filter)
             }
              g <- list("chr" = callbyChrPlot(myvcf, myvcf@myfile))
              if ("DP" %in% myvcf@genoString) {
                  g$read_depth = readdepthPlot(myvcf, myvcf@myfile)
                  results$medianReadDepth <- readDepth(myvcf, vcfparams@readdepth.target)
                  results$percentInTarget <- percentInTarget(myvcf,vcfparams@readdepth.target, vcfparams@readdepth.percent.limits)
              }
              if ("GT" %in% myvcf@genoString) {

                  if (myvcf@chunked == FALSE) {   
                      g$variant_type = calltypePlot(myvcf, myvcf@myfile)
                      g$homref = homrefPlot(myvcf, myvcf@myfile)
                  }
                      results$percentHets <- percentHets(myvcf, results$numberOfHomRefs[[3]], vcfparams@percenthets.limits[1], vcfparams@percenthets.limits[2])
                  results$numberOfHets <- numberOfHets(myvcf, vcfparams@het.limits[1], vcfparams@het.limits[2])
                  results$numberOfHomVars <- numberOfHomVars(myvcf, vcfparams@homvar.limits[1], vcfparams@homvar.limits[2])


              if (length(masked.ref) > 1) {                     
                  results$masked <- hetsMasked(myvcf, vcfparams@masked.limits, vcfparams@non.masked.limits, masked.ref)
              }
              results$hetGap <- hetGap(myvcf, vcfparams@het.gap.limits)

              }
              if (class(gold.ref) != "logical" && gold.ref@rare == TRUE) {
                  results$rare <- rareCompare(myvcf, gold.ref@track.rare, gold.ref@goldparams@percent.het.rare.limits)
              }
              if ("GQ" %in% myvcf@genoString) {
                  g$gq = genotypeQualityPlot(myvcf, myvcf@myfile)
                  results$meanGQ <- meanGQ(myvcf, vcfparams@gq.limit)
              }
                  
    if (class(cds.ref) != "logical") {
        results$titv <- titv(myvcf, cds.ref, vcfparams@titv.noncoding.limits[1], vcfparams@titv.noncoding.limits[2], vcfparams@titv.coding.limits[1], vcfparams@titv.coding.limits[2])}
    if (class(gold.ref) != "logical") {
        results$gold.ref <- goldCompare(myvcf, gold.ref@track, cds.ref, gold.ref@goldparams@titv.confirmed.limits, gold.ref@goldparams@titv.unconfirmed.limits, gold.ref@goldparams@percent.confirmed.limits, "Gold comparator")
    }

              results$MULTI <- list("Multi Calls", TRUE, myvcf@n.dup)
              results.as.list <- do.call(Map, c(c,results))
              tests <- (unlist((results.as.list)[2]))
              printnames <- (unlist((results.as.list)[1]))
              results <- (unlist((results.as.list)[3]))
              .Object@tests <- tests
              .Object@printnames <- printnames
              .Object@results <- results
              .Object@plots <- g
              return(.Object)

}

#############################################################################################
#############################################################################################

#Methods.  All methods are private and called through the VCFEvaluate constructor
#############################################################################################
#' Total Calls.  This is the total number of calls in the file (including MULTIs so hom ref, hom var and hom alt might not add up).  These methods are private.  Users are not expected to provide the number of hom ref and number of duplicate calls - these functions are generally called through the function VCFEvaluate
#' @return list with name, logical (passed or not), numeric (total number of calls) 
#' @param x myvcf of interest
#' @param n.hom Number of homozygous reference.  Pass so don't need to recalculate
#' @param n.dup number of multi calls to add onto the total calls
#' @param ll lower limit to pass for hom ref
#' @param ul upper limit to pass for hom ref
#' @keywords internal
setGeneric(name="numberCalls", def=function(x, n.hom, n.dup, ll, ul) {standardGeneric("numberCalls")})

setMethod("numberCalls",
          signature=c(x="VCFData", n.hom="numeric", n.dup = "numeric", ll="numeric", ul="numeric"),
          definition=function(x,n.hom, n.dup, ll, ul) {
              pass = FALSE
              n = n.hom + n.dup+ length(x@vr[x@vr$variant %in% c(1,2),])
              if (ll < n && ul > n) {
                  pass=TRUE
              }
              return(list("Number of Calls", pass, n))
          }
          )

#' Count Number of Hom Ref
#' @return list with name, logical (passed or not), numeric (number of hom ref) 
#' @param x myvcf of interest
#' @param ll lower limit to pass for hom ref
#' @param ul upper limit to pass for hom ref
#' @keywords internal
setGeneric(name="numberOfHomRefs", def=function(x, ll, ul) {standardGeneric("numberOfHomRefs")})

setMethod("numberOfHomRefs",
          signature=c(x="VCFData", ll="numeric", ul="numeric"),
          definition=function(x,ll, ul) {
              pass = FALSE
              if (x@chunked == FALSE) {
              if ("END" %in% x@infoString) {
                  n = sum(as.numeric(width(x@vr[x@vr$variant %in% 0,])))

              }
              else {
                  n = length(x@vr[x@vr$variant %in% 0,])
              }
          }
              if (x@chunked == TRUE) {
                  n = x@n.homref
              }
              if (ll < n && ul > n) {
                  pass=TRUE
              }
              return(list("Number of Homozygous Reference", pass, n))
          }
          )

#' Count Number of Hets
#' @return list with name, logical (passed or not), numeric (number of hets) 
#' @param x myvcf of interest
#' @param ll lower limit to pass for hets
#' @param ul upper limit to pass for hets
#' @keywords internal
setGeneric(name="numberOfHets", def=function(x,ll, ul) {standardGeneric("numberOfHets")})

setMethod("numberOfHets",
          signature=c(x="VCFData",ll="numeric", ul="numeric"),
          definition=function(x,ll, ul) {
              pass = FALSE
              n <- length(x@vr[mcols(x@vr)$variant %in% 1,])
              if (ll < n && ul > n) {
                  pass=TRUE
              }
              return(list("Number of Hets", pass, n))
          }
          )

#' Count Number of Hom Vars
#' @return list with name, logical (passed or not), numeric (number of hom vars) 
#' @param x myvcf of interest
#' @param ll lower limit to pass for Hom Vars
#' @param ul upper limit to pass for Hom Vars
#' @keywords internal
setGeneric(name="numberOfHomVars", def=function(x, ll, ul) {standardGeneric("numberOfHomVars")})

setMethod("numberOfHomVars",
          signature=c(x="VCFData", ll="numeric", ul="numeric"),
          definition=function(x, ll, ul) {
              pass = FALSE
              vr.homalt <- x@vr[mcols(x@vr)$variant %in% 2,]
              n <- length(vr.homalt)
              if (ll < n && ul > n) {
                  pass=TRUE
              }
              return(list("Number of Homozygous Alternatives", pass, n))
          }
          )


#' Percent of Hets as Total number of variants
#' @return list with name, logical (passed or not), numeric (percent hets) 
#' @param x myvcf of interest
#' @param n.hom Pass the number of homozygous reference so don't need to recalculate
#' @param ll lower limit to pass for percent
#' @param ul upper limit to pass for percent
#' @keywords internal
setGeneric(name="percentHets", def=function(x, n.hom, ll, ul) {standardGeneric("percentHets")})

setMethod("percentHets",
          signature=c(x="VCFData",n.hom="numeric", ll="numeric", ul="numeric"),
          definition=function(x, n.hom, ll, ul) {
              pass = FALSE
              n.het <- length(x@vr[mcols(x@vr)$variant %in% 1,])
              n.nonref <- length(x@vr[x@vr$var.bin == 1,])
              het.p <- n.het/(n.nonref + n.hom)
              if (n.het > 0 && ll < het.p && ul > het.p) {
                  pass=TRUE
              }
              
              return(list("Percent of Hets out of Total Number of Variants", pass, het.p))
          }
          )



#' Transition transversion ratio in coding and non-coding
#' @return list with name, logical vector (passed or not), numeric vector of transition transversion ratio in coding and non coding regions
#' @param x myvcf of interest
#' @param coding GRanges of coding regions
#' @param ll.c lower limit to pass for transition transversion in coding
#' @param ul.c upper limit to pass for transition transversion in coding
#' @param ll.nc lower limit to pass for transition transversion in noncoding
#' @param ul.nc upper limit to pass for transition transversion in noncoding
#' @keywords internal
setGeneric(name="titv", def=function(x, coding, ll.nc, ul.nc, ll.c, ul.c) {standardGeneric("titv")})

setMethod("titv",
          signature=c(x="VCFData",coding="GRanges", ll.nc="numeric", ul.nc="numeric", ll.c="numeric", ul.c="numeric"),
          definition=function(x, coding, ll.nc, ul.nc, ll.c, ul.c) {
              pass = c(FALSE, FALSE)
              x@vr$coding <- x@vr %over% coding
              df <- transform(GenomicRanges::as.data.frame(x@vr[x@vr$var.bin == 1,]), change = paste(ref, alt, sep="/"))
              curtitv <- xtabs(~ change + coding, df)
              titv.ratios <- apply(curtitv, 2, computeTiTv)
              if (length(titv.ratios) < 2) {
                  titv.ratios <- c(titv.ratios[1], NA)
              }

              if (length(titv.ratios) < 1) {
                  titv.ratios <- c(NA, NA)
              }

              names(titv.ratios) <- c("non-coding", "coding")
              if (!is.na(titv.ratios[1])) {
              if (ll.nc < titv.ratios[1] && ul.nc > titv.ratios[1]) {
                  pass[1] = TRUE
              }}
              if (!is.na(titv.ratios[2])) {
                  if (ll.c < titv.ratios[2] && ul.c > titv.ratios[2]) {
                  pass[2] = TRUE
              }
              }
              names(pass) <- c("titv_non_coding", "titv_coding")
              return(list("Transition Transversion in Non-Coding and Coding", pass, titv.ratios))
          }
          )



#' Comparator to gold standard
#' @return list with name, logical vector (passed or not), numeric vector of transition transversion in coding and non coding regions in confirmed and non confirmed
#' @param x myvcf of interest
#' @param gr gold GRanges
#' @param coding Coding GRanges
#' @param limits.c all confirmed limits
#' @param limits.uc all unconfirmed limits
#' @param percent.c percent confirmed limit
#' @param name.of.gold user readable name
#' @keywords internal
setGeneric(name="goldCompare", def=function(x, gr, coding, limits.c, limits.uc, percent.c, name.of.gold) {standardGeneric("goldCompare")})

setMethod("goldCompare",
          signature=c(x="VCFData", gr="GRanges", coding="GRanges", limits.c = "vector", limits.uc="vector", name.of.gold="character"),
          definition=function(x, gr, coding, limits.c, limits.uc, percent.c, name.of.gold) {
              pass = rep(FALSE, 5)
              #seqlengths(gr) = seqlengths(vr)[names(seqlengths(gr))]
              x@vr$confirmed <- x@vr %over% gr
              x@vr$coding <- x@vr %over% coding
              counts_overlap <- table(x@vr[x@vr$var.bin == 1,]$confirmed)
              col.true = which(names(counts_overlap) == TRUE)
              percent.overlap <- counts_overlap[col.true]/sum(counts_overlap)
              #seqlengths(coding) = seqlengths(vr)[names(seqlengths(coding))]

              df <- transform(GenomicRanges::as.data.frame(x@vr[x@vr$var.bin == 1,]), change = paste(ref, alt, sep="/"))
              mutations <- xtabs(~ change + coding, subset(df,df$confirmed==1))
              titv.confirmed <- apply(mutations, 2, computeTiTv)
              if (length(titv.confirmed) < 2) {
                  titv.confirmed <- c(NA, NA)
              }
              names(titv.confirmed) <- c("non-coding", "coding")
              
              mutations <- xtabs(~ change + coding, subset(df,df$confirmed==0))
              titv.unconfirmed <- apply(mutations, 2, computeTiTv)
              if (length(titv.unconfirmed) < 2) {
                  titv.unconfirmed <- c(NA, NA)
              }


              names(titv.unconfirmed) <- c("non-coding", "coding")
              if (!is.na(titv.unconfirmed[1])) {
              if (limits.uc[3] < titv.unconfirmed[1] && limits.uc[4] > titv.unconfirmed[1]) {
                  pass[1] = TRUE
              }
              if (limits.uc[1] < titv.unconfirmed[2] && limits.uc[2] > titv.unconfirmed[2]) {
                  pass[2] = TRUE
              }}
              if (!is.na(titv.confirmed[1])) {
              if (limits.c[3] < titv.confirmed[1] && limits.c[4] > titv.confirmed[1]) {
                  pass[3] = TRUE
              }
        
              if (limits.c[1] < titv.confirmed[2] && limits.c[2] > titv.confirmed[2]) {
                  pass[4] = TRUE
              }
          }

              if (length(col.true) != 0) {
              if (percent.overlap > percent.c) {
                  pass[5] = TRUE
              }}
              if (length(col.true) == 0) {
                  percent.overlap = NA
              }
              mysum <- c(titv.unconfirmed, titv.confirmed, percent.overlap)
              names(pass) <- c("titv_noncoding_unconfirmed", "titv_coding_unconfirmed", "titv_noncoding_confirmed", "titv_coding_confirmed", "percent_confirmed")
              names(mysum) <- names(pass)
              return(list(paste("Comparison to ", name.of.gold, sep=""), pass, mysum))
          }
          )


#' Comparator to rare variants.  Rare is defined as 0.01 percent or less
#' @return list with name, logical (passed or not), numeric (proportion hets in rare) 
#' @param x myvcf of interest
#' @param gr GRanges of gold comparator
#' @param percent.h what percent higher hets than hom_alt?
#' @keywords internal
setGeneric(name="rareCompare", def=function(x, gr, percent.h) {standardGeneric("rareCompare")})

setMethod("rareCompare",
          signature=c(x="VCFData", gr="GRanges", percent.h="numeric"),
          definition=function(x, gr, percent.h) {
              pass = FALSE
              #seqlengths(gr) = seqlengths(vr)[names(seqlengths(gr))]
              vr.rare <- x@vr[x@vr %over% gr,]
              num.het <- length(vr.rare[mcols(vr.rare)$variant %in% 1,])
              num.homalt <- length(vr.rare[mcols(vr.rare)$variant %in% 2,])
              num.var <- num.het + num.homalt
              percent.het = num.het/num.var

              if ((percent.het > percent.h) & (num.het > 0)) {
                  pass = TRUE
              }
              return(list(paste("Rare Variants"), pass, percent.het))
          }
          )


#' Median read depth
#' @return list with name, logical (passed or not), numeric median read depth
#' @param x myvcf of interest
#' @param coverage median coverage cutoff
#' @keywords internal
setGeneric(name="readDepth", def=function(x, coverage) {standardGeneric("readDepth")})

setMethod("readDepth",
          signature=c(x="VCFData", coverage="numeric"),
          definition=function(x, coverage) {
              pass = FALSE
              rd <- as.data.frame(cbind(totalDepth(x@vr[x@vr$var.bin == 1,]), mcols(x@vr[x@vr$var.bin == 1,])$GT), stringsAsFactors=FALSE)
              rd[,1] <- as.numeric(as.character(rd[,1]))
              names(rd) <- c("DP", "GT")
              med.cov <- median(rd[,1], na.rm=TRUE)
              if (!is.na(med.cov)) {
              if (med.cov > coverage) {
                  pass = TRUE
              }}
              return(list(paste("Median Read Depth"), pass, med.cov))
          }
          )


#' Number hets in masked GRanges
#' @return list with name, logical vector (passed or not), numeric vector number hets in masked region and not in masked region
#' @param x myvcf of interest
#' @param masked.limits Limits for self chain region (lower and upper)
#' @param non.masked.limits Limits for non self chain region (lower and upper)
#' @param masked GRanges of difficult regions (repeat masked, paralogs, self chain, etc)
#' @keywords internal
setGeneric(name="hetsMasked", def=function(x, masked.limits, non.masked.limits, masked) {standardGeneric("hetsMasked")})

setMethod("hetsMasked",
          signature=c(x="VCFData", masked.limits="numeric", non.masked.limits="numeric", masked="GRanges"),
          definition=function(x, masked.limits, non.masked.limits, masked) {
              pass = c(FALSE, FALSE)
              x@vr$masked <- x@vr %over% masked
              n.het <- length(x@vr[mcols(x@vr)$variant %in% 1,])
              het.masked <- length(x@vr[mcols(x@vr)$masked == TRUE & mcols(x@vr)$variant %in% 1,])/n.het
              het.not.masked <-length(x@vr[mcols(x@vr)$masked == FALSE & mcols(x@vr)$variant %in% 1,])/n.het
              if (!is.na(het.masked)) {
              if (het.masked > masked.limits[1] && het.masked < masked.limits[2]) {
                  pass[1] = TRUE
              }}
              if (!is.na(het.not.masked)) {
              if (het.not.masked > non.masked.limits[1] && het.not.masked < non.masked.limits[2]) {
                  pass[2] = TRUE
              }}
              names(pass) <- c("perc.het.masked", "perc.het.not.masked")
              returnMe <- c(het.masked, het.not.masked)
              names(returnMe) <- c("het.masked", "het.not.masked")
              return(list(paste("Percent Masked Regions Het Count"), pass, returnMe))
          }
          )


#' percent in target range read depth
#' For 15 to 60 for 30x (50 percent to 200 percent)
#' @return list with name, logical (passed or not), numeric percent in target range
#' @param x myvcf of interest
#' @param target Coverage target
#' @param percent.target From params, the percent required to be in target
#' @keywords internal
setGeneric(name="percentInTarget", def=function(x, target, percent.target) {standardGeneric("percentInTarget")})

setMethod("percentInTarget",
          signature=c(x="VCFData", target="numeric", percent.target="numeric"),
          definition=function(x, target, percent.target) {
              pass = FALSE

              lower <- target*.5
              upper <- target*2
              ontarget <- subset(x@vr[x@vr$var.bin == 1,], totalDepth(x@vr[x@vr$var.bin == 1,]) > lower & totalDepth(x@vr[x@vr$var.bin == 1,]) < upper)
              sample.percent.target <- 100*(length(ontarget)/length(x@vr[x@vr$var.bin == 1,]))
              if (!is.na(sample.percent.target)) {
              if (sample.percent.target >= percent.target) {
                  pass=TRUE
              }}
              return(list(paste("Percent in Target Read Depth Range"), pass, sample.percent.target))
          }
          )

#' Gap between HETs by Chromosome
#' @return list with name, logical vector (passed or not), numeric vector maximum gap between two het calls by chromosome
#' @param x myvcf of interest
#' @param myLimits is largest allowable gap by chromosome
#' @keywords internal
setGeneric(name="hetGap", def=function(x, myLimits) {standardGeneric("hetGap")})

setMethod("hetGap",
          signature=c(x="VCFData", myLimits="numeric"),
          definition=function(x, myLimits="numeric") {
              myChrs <- c(seq(1:22), "X")
              gaps <- NULL
              for (i in 1:23) {
                  vr.sub <- subset(x@vr[x@vr$variant == 1,], seqnames(x@vr[x@vr$variant == 1,]) == myChrs[i])
                  curgap = NA
                  if (length(vr.sub) > 0) {
                      vr.sub$gap <- start(vr.sub)- c(start(vr.sub)[1], start(vr.sub)[1:length(vr.sub)-1]) + 1
                      curgap = max(mcols(vr.sub)$gap)
                  }
                  gaps <- c(gaps, curgap)
              }
              pass = rep(TRUE, 23)
              for (i in 1:23) {
                  if (!is.na(gaps[i]) && (gaps[i] > myLimits[i])) {
                      pass[i] = FALSE
                  }
              }
              chrnames <- paste("chr", myChrs, sep="")
              names(pass) <- chrnames
              names(gaps) <- chrnames
              return(list(paste("Largest Gap between Hets by Chr"), pass, gaps))
          }
          )


#' Mean Genotype Quality (GQ)
#' @return list with name, logical (passed or not), numeric (mean genotype quality)
#' @param x myvcf of interest 
#' @param myLimit From params, the bounds for Genotype Quality
#' @keywords internal
setGeneric(name="meanGQ", def=function(x, myLimit) {standardGeneric("meanGQ")})

setMethod("meanGQ",
          signature=c(x="VCFData", myLimit="numeric"),
          definition=function(x, myLimit) {
              pass = FALSE
              myMean = mean(mcols(x@vr[x@vr$var.bin == 1,])$GQ)
              if (!is.na(myMean)) {
              if (myMean > myLimit) {
                  pass=TRUE
              }}
              return(list(paste("Mean GQ"), pass, myMean))
          }
          )


#' Private function to calc likelihood for admixture
#' @return numeric, Estimates of log likelihood of ancestry coefficients
#' @param data current data frame
#' @param par parameters to estimate
#' @keywords internal
myf <- function(data, par) {
    g <- NULL
    f1 <- NULL
    f2 <- NULL
    f3 <- NULL
	with(data, -sum((g*(log(f1*(par[1]/sum(par)) + f2*(par[2]/sum(par)) + f3*(par[3]/sum(par)))) + (2-g)*log((1-f1)*(par[1]/sum(par)) + (1-f2)*(par[2]/sum(par)) + (1-f3)*(par[3]/sum(par)) ) )))

	}

#' Private function to calc coefficients for admixture
#' @return Ancestry coefficients
#' @param dat current data frame
#' @param seed seed to start constroptim
#' @keywords internal
getCoefs <- function(seed, dat) {
    amat <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)
    bvec <- c(0,0,0)
    #g is the minor allele count
    amat%*%seed - bvec
    result <- constrOptim(seed, data=dat, myf, NULL, ui=amat, ci=bvec)
    q1 <- result$par[1]/sum(result$par)
    q2 <- result$par[2]/sum(result$par)
    q3 <- result$par[3]/sum(result$par)
return(c(q1, q2, q3))
	}

#' admixture - estimate admixture components using supervised ADMIXTURE algorithm.
#' @return List, name of metric (admixture), list of true false values if passed, ancestry coefficient values
#' @param x myvcf of interest
#' @param vr.ref reference file for admixture estimation (uses subset of 1000 Genomes)
#' @keywords internal
setGeneric(name="admixture", def=function(x, vr.ref) {standardGeneric("admixture")})

setMethod("admixture",
          signature=c(x="VCFData", vr.ref="VRanges"),
          definition=function(x, vr.ref) {
              #overlaps.index <- findOverlaps(x@vr, vr.ref)
              #overlaps <- x@vr[queryHits(overlaps.index),]
              #overlaps.y <- vr.ref[subjectHits(overlaps.index),]
              ##overlaps$keysort <- paste(seqnames(overlaps), start(overlaps), sep=":")
              #overlaps.y$keysort <- paste(seqnames(overlaps.y), start(overlaps.y), sep=":")
              #overlaps <- overlaps[order(overlaps$keysort),]
              #overlaps.y <- overlaps.y[order(overlaps.y$keysort),]
              #overlaps$EAS_AF <- overlaps.y$EAS_AF
              #overlaps$AFR_AF <- overlaps.y$AFR_AF
              #overlaps$EUR_AF <- overlaps.y$EUR_AF
#dat <- data.frame(g=overlaps$variant, f1=overlaps$EAS_AF, f2=overlaps$AFR_AF, f3=overlaps$EUR_AF)

              overlaps <- x@vr[x@vr %over% vr.ref,]
overlaps$key <- paste(seqnames(overlaps), start(overlaps), end(overlaps), sep=":")
dat.o <- GenomicRanges::as.data.frame(overlaps)
vr.ref$key <- paste(seqnames(vr.ref), start(vr.ref), end(vr.ref), sep=":")
dat.ref <- GenomicRanges::as.data.frame(vr.ref)
dat.m <- merge(dat.o, dat.ref, by.x="key", by.y="key")
dat <- data.frame(g=dat.m$variant, f1=dat.m$EAS_AF, f2=dat.m$AFR_AF, f3=dat.m$EUR_AF)


seed <- c(.8, .1, .1)              
mycoefs <- getCoefs(seed, dat)
seed <- c(.1, .1, .8)
mycoefs <- rbind(mycoefs, getCoefs(seed, dat))
seed <- c(.1, .8, .1)
          mycoefs <- rbind(mycoefs, getCoefs(seed, dat))
returnMe <- apply(mycoefs, 2, mean)             
names(returnMe) <- c("EAS", "AFR", "EUR")
return(list("Admixture", rep(TRUE,6), returnMe))
        }
          )


#############################################################################################
#############################################################################################
#Plots, the plots are generated.  No cutoffs applied
##############################################################################################
#' Histogram of read depth by GT
#' @return ggplot of histogram of read depth
#' @param x myvcf of interest
#' @param fn The filename (for title of the plot)
#' @keywords internal
setGeneric(name="readdepthPlot", def=function(x, fn) {standardGeneric("readdepthPlot")})
setMethod("readdepthPlot",
          signature=c(x="VCFData", fn="character"),
          definition=function(x, fn) {
              DP <- NULL
              GT <- NULL
              g <- NULL
              rd <- as.data.frame(cbind(as.numeric(totalDepth(x@vr[x@vr$var.bin == 1,])), mcols(x@vr[x@vr$var.bin == 1,])$GT))
              rd[,1] <- as.numeric(as.character(rd[,1]))
              names(rd) <- c("DP", "GT")
              if (nrow(rd) > 0) {
                  g <- ggplot(rd, aes(x=DP)) + geom_histogram() + facet_wrap(~ GT) + ggtitle(paste("Read Depth by Genotype", fn))
              }
              return(g)

          }
          )

#' Bar plot of variants (counts)
#' @return ggplot bar plot of variant types
#' @param x myvcf of interest
#' @param fn filename for plot title
#' @keywords internal
setGeneric(name="calltypePlot", def=function(x, fn) {standardGeneric("calltypePlot")})
setMethod("calltypePlot",
          signature=c(x="VCFData", fn="character"),
          definition=function(x, fn) {
              variant_type <- NULL
              g <- NULL
              variants <- as.data.frame(x@vr$GT)
              names(variants) <- "variant_type"
              if (nrow(variants) > 0){
                  g <- ggplot(variants, aes(factor(variant_type))) + geom_bar() + ggtitle(paste("Counts of Variant Types", fn))
              }
              return(g)

          }
          )

#' Histogram of genotype qualities
#' @return genotype quality histogram as ggplot
#' @param x myvcf of interest
#' @param fn filename for plot title
#' @keywords internal
setGeneric(name="genotypeQualityPlot", def=function(x, fn) {standardGeneric("genotypeQualityPlot")})
setMethod("genotypeQualityPlot",
          signature=c(x="VCFData", fn="character"),
          definition=function(x, fn) {
              GQ <- NULL
              g <- NULL
              variants <- as.data.frame(x@vr[x@vr$var.bin == 1,]$GQ)
              names(variants) <- "GQ"
              if (nrow(variants) > 0) {
                  g <- ggplot(variants, aes(x=GQ)) + geom_histogram() + ggtitle(paste("Genotype Quality", fn))
              }
              return(g)

          }
          )

#' Dot plot of variant call counts (hom alt and het) by chromosome
#' @return ggplot of variant calls by chromosome
#' @param x myvcf of interest
#' @param fn filename for plot title
#' @keywords internal
setGeneric(name="callbyChrPlot", def=function(x, fn) {standardGeneric("callbyChrPlot")})
setMethod("callbyChrPlot",
          signature=c(x="VCFData", fn="character"),
          definition=function(x, fn) {
              dat <- table(seqnames(x@vr[x@vr$var.bin == 1,]))
              myNames <- names(dat)
              Chromosome <- NULL
              Count <- NULL
              plotMe <- as.data.frame(cbind(c(as.character(seq(1:22)), "X", "Y"), rep(0,24)))
              plotMe[,2] <- as.integer(plotMe[,2]) -1 
              for (i in 1:24) {
                  curChr = plotMe[i,1]
                  myMatch = match(as.character(curChr), myNames)
                  if (!is.na(myMatch)) {
                      plotMe[i, 2] <- dat[myMatch]
              }}
names(plotMe) <- c("Chromosome", "Count")
plotMe$Count <- as.numeric(as.character(plotMe$Count))
plotMe$Chromosome <- factor(plotMe$Chromosome, levels=plotMe$Chromosome)
g <- ggplot(plotMe, aes(x=Chromosome, y=Count)) + geom_point(color="blue") + ggtitle(fn)  + 
    geom_text(data=plotMe,aes(x=Chromosome,y=Count,label=Count),vjust=0, angle=90, color="black") + scale_y_continuous("Hom Alt + Het")
              return(g)
          }
          )


#' Dot plot of variant call counts (hom ref) by chromosome
#' @return ggplot of homozygous reference counts
#' @param x myvcf of interest
#' @param fn filename for plot title
#' @keywords internal
setGeneric(name="homrefPlot", def=function(x, fn) {standardGeneric("homrefPlot")})
setMethod("homrefPlot",
          signature=c(x="VCFData", fn="character"),
          definition=function(x, fn) {
              myNames <- table(seqnames(x@vr[x@vr$var.bin == 0,]))
              myNames <- names(myNames[myNames > 0])
              mycounts <- sapply(myNames, function(X) {sum(width(subset(x@vr, seqnames(x@vr) == X)))})
              #print(mycounts)
              Chromosome <- NULL
              Count <- NULL
              plotMe <- as.data.frame(cbind(c(as.character(seq(1:22)), "X", "Y"), rep(0,24)))
              plotMe[,2] <- as.integer(plotMe[,2]) -1 
              for (i in 1:24) {
                  curChr = plotMe[i,1]
                  myMatch = match(as.character(curChr), myNames)
                  if (!is.na(myMatch)) {
                      plotMe[i, 2] <- mycounts[myMatch]
              }}
names(plotMe) <- c("Chromosome", "Count")
plotMe$Count <- as.numeric(as.character(plotMe$Count))
plotMe$Chromosome <- factor(plotMe$Chromosome, levels=plotMe$Chromosome)
g <- ggplot(plotMe, aes(x=Chromosome, y=Count)) + geom_point(color="blue") + ggtitle(fn)  + 
    geom_text(data=plotMe,aes(x=Chromosome,y=Count,label=Count),vjust=0, angle=90, color="black") + scale_y_continuous("Hom Ref")
              return(g)
          }
          )
#############################################################################################
#############################################################################################

#Getters for VCFEvaluate class.   
#############################################################################################

#' Getter for VCFEvaluate class to check if Sample Passed.  Using thresholds from VCFQAParam object return a list.  First return whether each test was passed (TRUE) or failed (FALSE).  Then return an overall pass (TRUE) or fail (FALSE).
#' @return True or False if sample passed all thresholds
#' @export
#' @param Object an object of type VCFQAReport
#' @examples
#' vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
#' mydir <- paste(dirname(vcffn), "/", sep="")
#' myfile <-basename(vcffn)
#' svp <- ScanVcfParam(which=GRanges("22", IRanges(0,200e5)), geno="GT")
#' vcfparams <- VCFQAParam(count.limits=c(3014580000, Inf), readdepth.target = 30)
#' vcf <- ReadVCFData(mydir, myfile, "GRCh38")
#' ev <- VCFEvaluate(vcf, vcfparams)
#' didSamplePassOverall(ev)
didSamplePassOverall <- function(Object) {
              if (all(slot(Object, "tests"))) {
              cat("SAMPLE PASSED \n")
              return(TRUE)
          }
          else {cat("SAMPLE FAILED \n")
            return(FALSE)}
          }
          

#' Getter for VCFEvaluate class to check if Sample Passed.  Using thresholds from VCFQAParam object return a list.  First return whether each test was passed (TRUE) or failed (FALSE).  Then return an overall pass (TRUE) or fail (FALSE).
#' @export
#' @return Vector of True and False
#' @param Object an object of type VCFQAReport
#' @examples
#' vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
#' mydir <- paste(dirname(vcffn), "/", sep="")
#' myfile <-basename(vcffn)
#' svp <- ScanVcfParam(which=GRanges("22", IRanges(0,200e5)), geno="GT")
#' vcfparams <- VCFQAParam(count.limits=c(3014580000, Inf), readdepth.target = 30)
#' vcf <- ReadVCFData(mydir, myfile, "GRCh38")
#' ev <- VCFEvaluate(vcf, vcfparams)
#' didSamplePass(ev)
didSamplePass <- function(Object) {
          return(slot(Object, "tests"))      
                     }
 
#' Getter for VCFQAReport class to return results.  Return a list showing values that the sample was evaluated on.
#' @param Object an object of type VCFQAReport
#' @return numeric vector of results
#' @export
#' @examples
#' vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
#' mydir <- paste(dirname(vcffn), "/", sep="")
#' myfile <-basename(vcffn)
#' svp <- ScanVcfParam(which=GRanges("22", IRanges(0,200e5)), geno="GT")
#' vcfparams <- VCFQAParam(count.limits=c(3014580000, Inf), readdepth.target = 30)
#' vcf <- ReadVCFData(mydir, myfile, "GRCh38")
#' ev <- VCFEvaluate(vcf, vcfparams)
#' getResults(ev)
getResults <- function(Object) {
              return(Object@results)
          }

#' Getter for VCFQAReport class to return plots slot.
#' @param Object Object of Class VCFQAReport
#' @return List of named ggplots
#' @export
#' @examples
#' vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
#' mydir <- paste(dirname(vcffn), "/", sep="")
#' myfile <-basename(vcffn)
#' svp <- ScanVcfParam(which=GRanges("22", IRanges(0,200e5)), geno="GT")
#' vcfparams <- VCFQAParam(count.limits=c(3014580000, Inf), readdepth.target = 30)
#' vcf <- ReadVCFData(mydir, myfile, "GRCh38")
#' ev <- VCFEvaluate(vcf, vcfparams)
#' getPlots(ev)
getPlots <- function(Object) {
              return(Object@plots)
          }


#' Getter for VCFQAReport class to return filename slot
#' @param Object Object of class VCFQAReport
#' @return Name of file
#' @export
#' @examples
#' vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
#' mydir <- paste(dirname(vcffn), "/", sep="")
#' myfile <-basename(vcffn)
#' svp <- ScanVcfParam(which=GRanges("22", IRanges(0,200e5)), geno="GT")
#' vcfparams <- VCFQAParam(count.limits=c(3014580000, Inf), readdepth.target = 30)
#' vcf <- ReadVCFData(mydir, myfile, "GRCh38")
#' ev <- VCFEvaluate(vcf, vcfparams)
#' getName(ev)
getName <- function(Object) {
              return(Object@fn)
          }

#############################################################################################
#############################################################################################

#Helper functions
#############################################################################################

#' Private function to calc transition tranversion (titv) ratio
#' @return Transition Transversion Ratio
#' @param x Table of SNPs
#' @keywords internal
computeTiTv <- function(x) {
      ti <- c("A/G", "G/A", "T/C", "C/T")
      sum(x[ti]) / sum(x[!names(x) %in% ti])
  }


