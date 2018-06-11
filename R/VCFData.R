#Class to Read in VCF File

#set the ggplot S3 class
setOldClass("gg")


#' Declare class
#' Reads in VCF using readVCFAsVRanges
#' @import VariantAnnotation
#' @import BiocGenerics
#' @import GenomicRanges
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom graphics pairs
#' @importFrom stats prcomp
#' @import rtracklayer
#' @import GenomeInfoDb
#' @import IRanges
#' @import methods
#' @import BiocParallel
#'
#' @return Object of class VCFData
#' @param mydir Directory of vcf file
#' @param myfile Filename of vcf file
#' @param vr.homref All SNPs from VCF with INDELs, MULTIs (seperately removed for variant and non variant), weird chromosomes removed
#' @param genoString A character vector of all genotype fields present (looks for AD, GQ, GT, DP)
#' @param infoString A character vector looking for "END" tag indicating file is a gVCF
#' @param genome Declare if the genome is GRCh37 or GRCh38
#' @param n.dup Counts the number of MULTIs removed
#' @param chunked Whether data was read in using ReadVCFDataChunk which means hom refs not in the admixture file were dropped

setClass(Class="VCFData", representation=representation(mydir="character", myfile="character", genome="character", vr = "VRanges",
                                 genoString = "character", infoString = "character", n.dup="numeric", chunked="logical", n.homref="numeric"),
)



#' User Constructor for class.  Calls VCFData constructor:
#' ReadVCFData is a wrapper for readVcfAsVRanges.  It removes indels, GL chromosomes, and MULTI calls.
#' It scans the header of the vcf file and adds in the following fields for analysis if present:  AD, GT, DP, GQ.
#' Looks for the "END" tag in the header and reads in file as gVCF if necessary.
#' @return Object of class VCFData
#' @export
#' @param mydir Directory of vcf file
#' @param myfile Filename of vcf file
#' @param genome GRCh37 or GRCh38
#' @examples
#' vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
#' mydir <- paste(dirname(vcffn), "/", sep="")
#' myfile <-basename(vcffn)
#' vcf <- ReadVCFData(mydir, myfile, "GRCh38")

ReadVCFData <- function(mydir, myfile, genome) {
    
    cat( "Reading VCF ... \n")
    .Object <- new(Class="VCFData", mydir=mydir, myfile=myfile, genome=genome)
    myfn <- paste(.Object@mydir, .Object@myfile, sep="")
    genoHeader = rownames(geno(scanVcfHeader(myfn)))
    genoString = NULL 
    if ("GT" %in% genoHeader) {
        genoString = c(genoString, "GT")}
    if ("DP" %in% genoHeader) {
        genoString = c(genoString, "DP")}
    if ("GQ" %in% genoHeader) {
        genoString = c(genoString, "GQ")}
    if ("AD" %in% genoHeader) {
        genoString = c(genoString, "AD")}
    .Object@genoString = genoString
    #validObject(.Object)

    infoString = NA_character_ 
    
    infoHeader = rownames(info(scanVcfHeader(myfn)))
    if ("END" %in% infoHeader)  {
        infoString =c("END")
    }
    .Object@infoString = infoString

    # svp <- ScanVcfParam(geno=genoString, info=infoString)
    #need to swap this single line if doing a pop vcf
    #svp <- ScanVcfParam(which=GRanges("1",IRanges(1,100e6)), geno=genoString, info=NA)
    #svp <- ScanVcfParam(which=GRanges("chr1",IRanges(1,100e5)), geno=genoString, info=infoString)
    seqprefix = ""
    if (length(grep("chr", rownames(meta(scanVcfHeader(myfn))$contig)[1])) > 0) {
        seqprefix = "chr"
    }
    vcf.chrs <- rownames(as.data.frame(seqinfo(scanVcfHeader(myfn))))
    #only want to bpl the chromosomes present
    reg.chrs <- NULL
    for (i in 1:22) {
        if ((paste("chr", as.character(i), sep="") %in% vcf.chrs) || (as.character(i) %in% vcf.chrs)) {
            reg.chrs <- c(reg.chrs, as.character(i))
        }}
    if (("chrX" %in% vcf.chrs) || ("X" %in% vcf.chrs)) {
    reg.chrs <- c(reg.chrs, "X")
    }

    if (("chrY" %in% vcf.chrs) || ("Y" %in% vcf.chrs)) {
    reg.chrs <- c(reg.chrs, "Y")
    }
    svp <- ScanVcfParam(geno=genoString, info=infoString)
 
    .Object@vr <- readVcfAsVRanges(myfn, genome=genome, param=svp)
    #remove any indels
    .Object@vr <- subset(.Object@vr, ref(.Object@vr) %in% c("A", "G", "T", "C"))
    .Object@vr <- subset(.Object@vr, alt(.Object@vr) %in% c("A", "G", "T", "C", "<NON_REF>"))
    #command no longer works
    #vr <- vr[isSNV(vr),]
    #relevel to UCSC and then remove chr so all consistent
    seqlevelsStyle(.Object@vr) = "NCBI"
    #get rid of all GLs
    #reg.chrs <- c(as.character(seq(1:22)))
    #reg.chrs <- c(as.character(seq(1:22)), "X", "Y")
    .Object@vr <- keepSeqlevels(.Object@vr, reg.chrs, pruning.mode="coarse")
    #keepStandardChromosomes(vr)
    hets <- c("0/1", "1|0", "0|1")
    hom_alts <- c("1/1", "1|1")
    hom_refs <- c("0/0", "0|0")

    
    #hets <- c("0/1", "1|0", "0|1", "0/2", "0/3", "0/4", "1/2", "1/4", "2/3", "2/4")
    #hom_alts <- c("1/1", "1|1", "2/2")
    #hom_refs <- c("0/0", "0|0")
    non_refs <- c(hets, hom_alts)

    .Object@vr <- .Object@vr[!(.Object@vr$GT %in% non_refs & (alt(.Object@vr) %in% "<NON_REF>")),]

    #keep highest quality MULTI
    .Object@vr$variant <- 99
    .Object@vr$var.bin <- 99

    if (length(.Object@vr[mcols(.Object@vr)$GT %in% hets,]) > 0) {
        .Object@vr[mcols(.Object@vr)$GT %in% hets,]$variant <- 1
           .Object@vr[mcols(.Object@vr)$GT %in% hets,]$var.bin <- 1
    }
    if (length(.Object@vr[mcols(.Object@vr)$GT %in% hom_alts,]) > 0) { 
        .Object@vr[mcols(.Object@vr)$GT %in% hom_alts,]$variant <- 2
           .Object@vr[mcols(.Object@vr)$GT %in% hom_alts,]$var.bin <- 1
    }
    if (length(.Object@vr[mcols(.Object@vr)$GT %in% hom_refs,]) > 0) {
        .Object@vr[mcols(.Object@vr)$GT %in% hom_refs,]$variant <- 0
           .Object@vr[mcols(.Object@vr)$GT %in% hom_refs,]$var.bin <- 0
    }
    .Object@vr$key <- paste(seqnames(.Object@vr), ":", start(.Object@vr), ":", end(.Object@vr), ":", .Object@vr$var.bin, sep="")
    dup.keys <- .Object@vr$key[duplicated(.Object@vr$key)]
    .Object@vr <- .Object@vr[!.Object@vr$key %in% dup.keys,]

    n.dup <- length(dup.keys)
    .Object@n.dup <- n.dup
    .Object@chunked = FALSE
    return(.Object)
}

#' chunkData is a private function to read in a chunk and process it.  This is a private function and is not meant to be called by the user.  An example is provided in line with bioconductor policies.
#' @return VRanges of single processed chromosome
#' @param myfn Filename/path of file
#' @param genome GRCh37 or GRCh38
#' @param svp Params specified for readVCFAsVRanges
#' @param admixture.ref VRanges with MAF for superpopulations (EAS, AFR, EUR)
#' @keywords internal 
chunkData <- function(myfn, genome, svp, admixture.ref) {
    vr <- readVcfAsVRanges(myfn, genome=genome, param=svp)
    vr.1 <- vr[1,]
    n.homref = 0
    n.dup = 0
    vr <- subset(vr, ref(vr) %in% c("A", "G", "T", "C"))
    vr <- subset(vr, alt(vr) %in% c("A", "G", "T", "C", "<NON_REF>"))
    seqlevelsStyle(vr) = "NCBI"
    seqlevelsStyle(vr.1) = "NCBI"

                                        #reg.chrs <- c(as.character(seq(1:22)))
    reg.chrs <- c(as.character(seq(1:22)), "X", "Y")
    #vr <- keepSeqlevels(vr, reg.chrs)
    #hets <- c("0/1", "1|0", "0|1", "0/2", "0/3", "0/4", "1/2", "1/4", "2/3", "2/4")   
    #hom_alts <- c("1/1", "1|1", "2/2")
    #hom_refs <- c("0/0", "0|0")

    hets <- c("0/1", "1|0", "0|1")
    hom_alts <- c("1/1", "1|1")
    hom_refs <- c("0/0", "0|0")
    non_refs <- c(hets, hom_alts)

    vr <- vr[!(vr$GT %in% non_refs & (alt(vr) %in% "<NON_REF>")),]

    if (length(vr) > 0) {
        

    vr$variant <- 99
    vr$var.bin <- 99
 

    if (length(vr[mcols(vr)$GT %in% hets,]) > 0) {
        vr[mcols(vr)$GT %in% hets,]$variant <- 1
           vr[mcols(vr)$GT %in% hets,]$var.bin <- 1
    }

    if (length(vr[mcols(vr)$GT %in% hom_alts,]) > 0) { 
        vr[mcols(vr)$GT %in% hom_alts,]$variant <- 2
           vr[mcols(vr)$GT %in% hom_alts,]$var.bin <- 1
    }
    if (length(vr[mcols(vr)$GT %in% hom_refs,]) > 0) {
        vr[mcols(vr)$GT %in% hom_refs,]$variant <- 0
           vr[mcols(vr)$GT %in% hom_refs,]$var.bin <- 0
    }
    n.dup <- 0
    if (length(vr) > 0) {        
    vr$key <- paste(seqnames(vr), ":", start(vr), ":", end(vr), ":", vr$var.bin, sep="")
    dup.keys <- vr$key[duplicated(vr$key)]
    vr <- vr[!vr$key %in% dup.keys,]
    n.dup <- length(dup.keys)
    }
   #count the number of hom refs and return it
   vr.1 <- vr[1,]
    if ("END" %in% vcfInfo(svp)) {
        n.homref = sum(as.numeric(width(vr[vr$variant %in% 0,])))
      }
     else {
        n.homref = length(vr[vr$variant %in% 0,])
    }

    vr$admix <- vr %over% admixture.ref
    vr.1 <- vr[1,]
    vr <- subset(vr, vr$admix == TRUE | vr$var.bin == 1)
  
}
    if (length(vr) == 0)  {
        vr = vr.1
        vr$variant <- 99
        vr$var.bin <- 99
        vr$key <- ""
        vr$admix <- FALSE
    }
    vr$n.homref = n.homref
    vr$n.dup = n.dup
    return(vr)
}


#' User Constructor for class.  Calls VCFData constructor:
#' ReadVCFDataChunk is a wrapper for readVcfAsVRanges.  It removes indels, GL chromosomes, and MULTI calls.
#' It scans the header of the vcf file and adds in the following fields for analysis if present:  AD, GT, DP, GQ.
#' Looks for the "END" tag in the header and reads in file as gVCF if necessary.  
#' This is a multi core version of readVCFData.  Note, input file must have been zipped and have a corresponding tabix file.  It will drop all hom ref sites not in the admixture file but retain the counts of homref and multi in the VCF file.  This means that a few of the metrics and the hom ref plot can no longer be calculated in VCFQAReport.  If the metrics can no longer be calculated, it will not be output.  Please note that if using a filter on the data (eg gq.filter) this will not be applied to the hom ref and total number of calls.  The filter is applied in the VCFQAReport step and the metrics number of hom ref and total number of calls is calculated while reading in the file.  When calling this function keep in mind the memory requirements.  For example, if numcores=6, then when submitting the job you may request 12 Gb each core (72 Gb total).  However the VCF in memory will need to fit back onto a single core or else R will not be able to allocate the memory.  The given example here does not make sense to run as it includes only chromosome 22.  
#' @return Object of type VCFData
#' @export
#' @param mydir Directory of vcf file
#' @param myfile Filename of vcf file (zipped)
#' @param genome GRCh37 or GRCh38
#' @param admixture.ref VRanges with MAF for superpopulations (EAS, AFR, EUR)
#' @param numcores Number of cores to read in VCF (passed to bplapply)
#' @examples
#' vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
#' mydir <- paste(dirname(vcffn), "/", sep="")
#' myfile <-basename(vcffn)
#' svp <- ScanVcfParam(which=GRanges("22", IRanges(0,1e5)), geno="GT")
#' vcf <- ReadVCFData(mydir, myfile, "GRCh38")
#' admix.var <- getVR(vcf)[getVR(vcf)$GT %in% c("0|1", "1|0", "1|1"),][,1:2]
#' admix.var$EAS_AF <- ifelse(admix.var$GT %in% c("1|1"), 1, .5)
#' admix.var$AFR_AF<- 0
#' admix.var$EUR_AF<- 0
#' admix.hom <- getVR(vcf)[getVR(vcf)$GT %in% c("0|0"),][,1:2]
#' admix.hom$EAS_AF<- 0
#' admix.hom$AFR_AF<- 1
#' admix.hom$EUR_AF<- 1
#' admix.ref <- c(admix.var, admix.hom)
#' ReadVCFDataChunk(mydir, myfile, "GRCh38", admix.ref, numcores=2)

ReadVCFDataChunk <- function(mydir, myfile, genome, admixture.ref, numcores) {
    
    cat( "Reading VCF ... \n")
    .Object <- new(Class="VCFData", mydir=mydir, myfile=myfile, genome=genome)
    myfn <- paste(.Object@mydir, .Object@myfile, sep="")
    genoHeader = rownames(geno(scanVcfHeader(myfn)))
    genoString = NULL 
    if ("GT" %in% genoHeader) {
        genoString = c(genoString, "GT")}
    if ("DP" %in% genoHeader) {
        genoString = c(genoString, "DP")}
    if ("GQ" %in% genoHeader) {
        genoString = c(genoString, "GQ")}
    if ("AD" %in% genoHeader) {
        genoString = c(genoString, "AD")}
    .Object@genoString = genoString
    .Object@genome <- genome

    infoString = NA_character_ 
    
    infoHeader = rownames(info(scanVcfHeader(myfn)))
    if ("END" %in% infoHeader)  {
        infoString =c("END")
    }
    #pull in  line of file
    myline = utils::read.table(myfn, nrows = 1, skip = 500, header=FALSE)
    seqprefix = ""
    if (!is.na(pmatch("chr",myline[,1]))) {
        seqprefix = "chr"
    }
    
    .Object@infoString = infoString
    vcf.chrs <- rownames(as.data.frame(seqinfo(scanVcfHeader(myfn))))
    #only want to bpl the chromosomes present
    reg.chrs <- NULL
    for (i in 1:22) {
        if ((paste("chr", as.character(i), sep="") %in% vcf.chrs) || (as.character(i) %in% vcf.chrs)) {
            reg.chrs <- c(reg.chrs, as.character(i))
        }}
    if (("chrX" %in% vcf.chrs) || ("X" %in% vcf.chrs)) {
    reg.chrs <- c(reg.chrs, "X")
    }

    if (("chrY" %in% vcf.chrs) || ("Y" %in% vcf.chrs)) {
    reg.chrs <- c(reg.chrs, "Y")
    }

    multicoreParam <- MulticoreParam(workers = numcores)
    .Object@vr <- do.call(c, (bplapply(reg.chrs, function(v) {
        svp <- ScanVcfParam(which = GRanges(paste(seqprefix, v, sep=""), IRanges(1,250e6)), geno=genoString, info=infoString)
        chunkData(myfn, genome, svp, admixture.ref) }
        , BPPARAM=multicoreParam)))
    .Object@n.dup <- sum(unique(data.frame(seqnames(.Object@vr), mcols(.Object@vr)$n.dup))[,2])
    .Object@n.homref <- sum(unique(data.frame(seqnames(.Object@vr), mcols(.Object@vr)$n.homref))[,2])
    .Object@chunked <- TRUE
    return(.Object)
}


#' getVr is a Getter. Returns vr slot.
#' @return VRanges
#' @export
#' @param x VCFData object
#' @examples
#' vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
#' mydir <- paste(dirname(vcffn), "/", sep="")
#' myfile <-basename(vcffn)
#' svp <- ScanVcfParam(which=GRanges("22", IRanges(0,1e5)), geno="GT")
#' vcf <- ReadVCFData(mydir, myfile, "GRCh38")
#' getVR(vcf)
getVR <- function(x) {
   return(x@vr)
}
