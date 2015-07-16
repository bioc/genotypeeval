#Class to Read in Gold File 

#' Declare class Gold to store information from Gold" (1000 Genomes for example) along with the GoldDataParam
#' @return Object of class GoldData
#' @import VariantAnnotation
#' @import BiocGenerics
#' @import GenomicRanges
#' @import ggplot2
#' @import rtracklayer
#' @import GenomeInfoDb
#' @import IRanges
#' @import methods
#' @import BiocParallel
#' @param genome Genome build, GRCh37 or GRCh38
#' @param track Where the gold data is stored
#' @param goldparams The Param file with the limits to be applied
#' @param track.rare Stores the Gold data with MAF < 0.01 if MAF exists

setClass(Class="GoldData", representation=representation(genome="character", track="GRanges", goldparams="GoldDataParam", rare="logical", track.rare="GRanges"),
)


#' User Constructor for class
#' @return Object of class GoldData
#' @export
#' @param genome Genome build, GRCh37 or GRCh38
#' @param vcffilename path and filename of vcf file
#' @param goldparams GoldDataParam object setting thresholds for evaluation
#' @examples
#' gparam <-  GoldDataParam(percent.confirmed=0.792, percent.het.rare = 0.93)
#' g1000fn <- system.file("ext-data", "example_gold_file.vcf", package="genotypeeval")
#' g1000 <- ReadGoldData("GRCh38", g1000fn, gparam)
ReadGoldData <- function(genome, vcffilename, goldparams) {
    cat( " Gold Data Constructor ... Reading Gold Data in (can take awhile) \n")   
    .Object <- new(Class="GoldData", genome=genome, goldparams=goldparams)
    infoHeader = rownames(info(scanVcfHeader(vcffilename)))
    rare = FALSE
    if ("AF" %in% infoHeader) {
        rare = TRUE
    }
    .Object@rare = rare
    .Object@track <- readVcfGold(vcffilename, rare, genome)
    .Object@track.rare <- .Object@track[1]
    if (rare==TRUE){
        .Object@track.rare <- GenomicRanges::subset(.Object@track, mcols(.Object@track)$AF < 0.01)
    }
    return(.Object)
}

#' User Constructor for class.  Used to associate the gold params object with the gold granges and to check if MAF is present.
#' @return Object of class GoldData
#' @export
#' @param genome Genome build, GRCh37 or GRCh38
#' @param gold.granges Gold file as GRanges 
#' @param goldparams GoldDataParam object setting thresholds for evaluation
#' @examples
#' gparam <-  GoldDataParam(percent.confirmed=0.792, percent.het.rare = 0.93)
#' gr <- GRanges(seqnames="22", IRanges(1e7,5e7))
#' gold <- GoldDataFromGRanges("GRCh38", gr, gparam)

GoldDataFromGRanges <- function(genome, gold.granges, goldparams) {
    cat( " Gold Data Constructor ... \n")   
    .Object <- new(Class="GoldData", genome=genome, goldparams=goldparams)
    .Object@track <- gold.granges

    .Object@goldparams <- goldparams
    rare = FALSE

    if ("AF" %in% names(mcols(gold.granges))) {
        rare = TRUE
    }
    .Object@rare = rare
    .Object@track.rare <- .Object@track[1]
    if (rare==TRUE){
        .Object@track.rare <- GenomicRanges::subset(.Object@track, mcols(.Object@track)$AF < 0.01)
    }
    return(.Object)
}




#' Private method for class.  Read in Gold file - will read in the AF if it is detected in header
#' @return GRanges, read in from gold file
#' @param vcfFile Gold file to read in 
#' @param rare Whether AF is in the header
#' @param genome Build GRCh37 or GRCh38
#' @keywords internal
readVcfGold <- function(vcfFile,rare, genome) {
    if (rare == TRUE) {
        svp <- ScanVcfParam(geno=NA, info="AF")
    }
    else {
        svp <- ScanVcfParam(geno=NA, info=NA)
    }
    vr <- readVcfAsVRanges(vcfFile, genome=genome, param=svp)
    seqlevelsStyle(vr) = "NCBI"
    reg.chrs <- c(as.character(seq(1:22)), "X", "Y")
    vr <- keepSeqlevels(vr, reg.chrs)
    #keepStandardChromosomes(vr)
    genome(vr) <- genome
    return(vr)
}

