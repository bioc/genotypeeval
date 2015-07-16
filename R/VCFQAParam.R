#Class to set thresholds for VCF File evaluation

#' Declare class VCFQAParam which will store thresholds to apply to VCFEvaluate object.  This is intended for use in batch mode when a large number of vcf files needs to be screened and individual vcf files that fail flagged.
#' All limits follow the format lower limit than upper limit
#' @return VCFQAParam object
#' @param homref.limits lower limit, upper limit, number of homozygous reference
#' @param het.limits lower limit, upper limit, number of heterozygous calls
#' @param homvar.limits lower limit, upper limit, number of homozygous alternative
#' @param count.limits lower limit, upper limit, total number of counts
#' @param percenthets.limits lower limit, upper limit, Number of Heterozgyous / (Total Number of Counts) or percent het
#' @param titv.noncoding.limits lower limit, upper limit, Transition transversion ratio in noncoding regions
#' @param titv.coding.limits lower limit, upper limit, Transition transversion ratio in coding regions
#' @param readdepth.target The sequencing depth target (eg 30x)
#' @param readdepth.limits lower limit, upper limit, Mean read depth
#' @param readdepth.percent.limits lower limit, upper limit, Percent read depth in target (50 percent to 200 percent of target read depth)
#' @param gq.limit lower limit, Mean genotype quality (does not make sense to have an upper limit)
#' @param masked.limits lower limit, upper limit, (Number of heterozygous in masked regions)/(Total number of heterozygotes)
#' @param non.masked.limits lower limit, upper limit, (Number of heterozygous in non-self chained regions)/(Total number of heterozygotes)
#' @param het.gap.limits lower limit, upper limit, Largest gap within chromosome between two heterozygous calls
#' @param gq.filter filter for the VCF file on genotype quality (eg only GQ > 90)
#' @param dp.filter filter for the VCF file on read depth (eg only DP > 0)
setClass(Class="VCFQAParam", representation=representation(homref.limits="numeric", het.limits="numeric", homvar.limits="numeric", percenthets.limits = "numeric",
                                       titv.noncoding.limits="numeric", titv.coding.limits="numeric", readdepth.target="numeric", readdepth.limits="numeric", readdepth.percent.limits="numeric", gq.limit="numeric", masked.limits="numeric", non.masked.limits="numeric", het.gap.limits="numeric", count.limits = "numeric", gq.filter = "numeric", dp.filter="numeric"),
               )


#' User Constructor for class.  Call limits are set as default to pass.
#' @return Object of class VCFQAParam
#' @export
#' @param homref.limits lower limit, upper limit, number of homozygous reference
#' @param het.limits lower limit, upper limit, number of heterozygous calls
#' @param homvar.limits lower limit, upper limit, number of homozygous alternative
#' @param count.limits lower limit, upper limit, total number of counts
#' @param percenthets.limits lower limit, upper limit, Number of Heterozgyous / (Total Number of Counts) or percent het
#' @param titv.noncoding.limits lower limit, upper limit, Transition transversion ratio in noncoding regions
#' @param titv.coding.limits lower limit, upper limit, Transition transversion ratio in coding regions
#' @param readdepth.target The sequencing depth target (eg 30x)
#' @param readdepth.limits lower limit, upper limit, Mean read depth
#' @param readdepth.percent.limits lower limit, upper limit, Percent read depth in target (50 percent to 200 percent of target read depth)
#' @param gq.limit lower limit, Mean genotype quality (does not make sense to have an upper limit)
#' @param masked.limits lower limit, upper limit, (Number of heterozygous in self chained regions)/(Total number of heterozygotes)
#' @param non.masked.limits lower limit, upper limit, (Number of heterozygous in non-self chained regions)/(Total number of heterozygotes)
#' @param het.gap.limits lower limit, upper limit, Largest gap within chromosome between two heterozygous calls
#' @param gq.filter filter for the VCF file on genotype quality (eg only GQ > 90)
#' @param dp.filter filter for the VCF file on read depth (eg only DP > 0)
#' @examples
#' vcfparams <- VCFQAParam(count.limits=c(3014580000, Inf), readdepth.target = 30)
VCFQAParam <- function(homref.limits=c(-Inf, Inf), het.limits=c(-Inf, Inf), homvar.limits=c(-Inf, Inf), percenthets.limits=c(-Inf, Inf), titv.noncoding.limits=c(-Inf, Inf), titv.coding.limits=c(-Inf, Inf),
                             readdepth.target=-1, readdepth.limits=c(-Inf, Inf), readdepth.percent.limits=0, gq.limit=0, masked.limits=c(-Inf, Inf), non.masked.limits=c(-Inf, Inf), het.gap.limits=rep(Inf,24), count.limits=c(-Inf, Inf), gq.filter = 0, dp.filter = -1)
    {
        .Object <- new(Class="VCFQAParam", homref.limits=homref.limits, het.limits=het.limits, homvar.limits=homvar.limits, percenthets.limits=percenthets.limits, titv.coding.limits = titv.coding.limits, titv.noncoding.limits = titv.noncoding.limits, readdepth.target= readdepth.target, readdepth.limits=readdepth.limits, readdepth.percent.limits=readdepth.percent.limits, gq.limit=gq.limit, 
masked.limits=masked.limits, non.masked.limits=non.masked.limits, het.gap.limits=het.gap.limits, count.limits = count.limits, gq.filter = gq.filter, dp.filter = dp.filter)
    return(.Object)
}



