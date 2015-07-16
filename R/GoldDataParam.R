#' Declare class GoldDataParam which will store thresholds to apply to VCFEvaluate object.  This is intended for use in batch mode when a large number of vcf files needs to be screened and individual vcf files that fail flagged.
#' All limits follow the format lower limit than upper limit
#'
#' @return Object of type GoldDataParam
#' @param titv.confirmed.limits lower limit coding, upper limit coding, lower limit noncoding, upper limit noncoding, Transition transversion ratios for confirmed snps
#' @param titv.unconfirmed.limits lower limit coding, upper limit coding, lower limit noncoding, upper limit noncoding, Transition transversion ratios for unconfirmed snps
#' @param percent.confirmed.limits lower limit, upper limit, percent confirmed in Gold comparator
#' @param percent.het.rare.limits lower limit, upper limit, (Percent Het in Rare, MAF < 0.01 in Gold) / Total number of Heterozygotes
setClass(Class="GoldDataParam", representation=representation(titv.confirmed.limits="numeric", titv.unconfirmed.limits="numeric", percent.confirmed.limits="numeric", percent.het.rare.limits = "numeric"),
         ) 



#' User Constructor for class
#' @export
#' @return Object of type GoldDataParam
#' @param titv.coding.confirmed.l Lower limit of transition transversion ratio in coding confirmed
#' @param titv.coding.confirmed.u upper limit of transition transverion ratio coding confirmed
#' @param titv.noncoding.confirmed.l Lower limit of transition transversion ratio in noncoding confirmed
#' @param titv.noncoding.confirmed.u upper limit of transition transverion ratio noncoding confirmed
#' @param titv.coding.unconfirmed.l Lower limit of transition transversion ratio in coding unconfirmed
#' @param titv.coding.unconfirmed.u upper limit of transition transverion ratio coding unconfirmed
#' @param titv.noncoding.unconfirmed.l Lower limit of transition transversion ratio in noncoding unconfirmed
#' @param titv.noncoding.unconfirmed.u upper limit of transition transverion ratio noncoding unconfirmed
#' @param percent.confirmed.limits lower limit, upper limit, percent confirmed in Gold comparator
#' @param percent.het.rare.limits lower limit, upper limit, (Percent Het in Rare, MAF < 0.01 in Gold) / Total number of Heterozygotes
#' @examples
#' gparam <- GoldDataParam(percent.confirmed=0.792, percent.het.rare = 0.93)
GoldDataParam <- function(titv.coding.confirmed.l=0, titv.coding.confirmed.u=5, titv.noncoding.confirmed.l=0, titv.noncoding.confirmed.u=5,
                      titv.coding.unconfirmed.l=0, titv.coding.unconfirmed.u=5, titv.noncoding.unconfirmed.l=0, titv.noncoding.unconfirmed.u=5, percent.confirmed.limits=0,
                          percent.het.rare.limits=0) {
    titv.confirmed.limits <- c(titv.coding.confirmed.l,titv.coding.confirmed.u,titv.noncoding.confirmed.l,titv.noncoding.confirmed.u) 
    titv.unconfirmed.limits <- c(titv.coding.unconfirmed.l,titv.coding.unconfirmed.u,titv.noncoding.unconfirmed.l,titv.noncoding.unconfirmed.u)
   .Object = new(Class="GoldDataParam", titv.confirmed.limits=titv.confirmed.limits, titv.unconfirmed.limits=titv.unconfirmed.limits,  percent.het.rare.limits = percent.het.rare.limits, percent.confirmed.limits = percent.confirmed.limits)
   return(.Object)
}





