context("Unit Tests QA")
vcffn <- system.file("ext-data", "chr22.GRCh38.vcf.gz", package="genotypeeval")
vcfsplit <- strsplit(vcffn,"/")[[1]]
mydir <- paste(paste(vcfsplit[1:length(vcfsplit)-1], collapse="/"), "/", sep="")
myfile <- vcfsplit[length(vcfsplit)]
vcf <- ReadVCFData(mydir, myfile, "GRCh38")

test_that("Basic Counts of Variant Type Values", {
    expect_equal(numberOfHomRefs(vcf, 0, 3e11)[[3]], 348)
    expect_equal(numberOfHets(vcf, 0, 3e11)[[3]], 699)
    expect_equal(numberOfHomVars(vcf, 0, 3e11)[[3]], 196)
    expect_equal(numberCalls(vcf,348,0, 0,3e11)[[3]], 1243)
})

test_that("Basic Counts of Variant Type Limits", {
    expect_equal(numberOfHomRefs(vcf, 400, 3e11)[[2]], FALSE)
    expect_equal(numberOfHets(vcf, 0, 200)[[2]], FALSE)
})


test_that("Percent Hets", {
    expect_equal(percentHets(vcf, 348, 0, 1)[[3]], 0.5623492, tolerance=1e-6)
})

test_that("Transition Transversion", {
    myCoding <- GRanges(seqnames="22", ranges=IRanges(49999999, width=100000))
    myOutput <- c(2.360000, 2.716667)
    names(myOutput) <- c("non-coding", "coding")
    expect_equal(titv(vcf, myCoding, 0, 5, 0, 5)[[3]], myOutput, tolerance=1e-6)
})

test_that("Read Depth", {
    expect_equal(readDepth(vcf, 30)[[3]], 30)
})

test_that("Het Gap", {
    myOutput = 100925
    names(myOutput) = "chr22"
    expect_equal(hetGap(vcf, rep(1,22))[[3]][22], myOutput)
})

test_that("Gold Compare", {
    myCoding <- GRanges(seqnames="22", ranges=IRanges(49999999, width=100000))
    myGold <- GRanges(seqnames="22", ranges=IRanges(49932468, 49932468))
    myOutput <- 0.001117318
    names(myOutput) <- "percent_confirmed"
    expect_equal(goldCompare(vcf, myGold, myCoding, c(0,5,0,5), c(0,5,0,5), 0.8, "gold")[[3]][5], myOutput, tolerance=1e-6)
})

test_that("Admixture", {
    admix.var <- vcf@vr[vcf@vr$GT %in% c("0|1", "1|0", "0/1", "1|1"),]
    admix.var$EAS_AF<- rep(1, length(admix.var))
    admix.var[admix.var$GT %in% c("0|1", "1|0"),]$EAS_AF <- .5
    admix.var$AFR_AF<- rep(0, length(admix.var))
    admix.var$EUR_AF<- rep(0, length(admix.var))
    admix.hom <- vcf@vr[vcf@vr$GT %in% c("0/0", "0|0"),]
    admix.hom$EAS_AF<- rep(0, length(admix.hom))
    admix.hom$AFR_AF<- rep(1, length(admix.hom))
    admix.hom$EUR_AF<- rep(1, length(admix.hom))
 
    admix.ref <- c(admix.var, admix.hom)
    admix.ref$variant <- NULL
    admix.ref$var.bin <- NULL
    myOutput <- c(1, 0, 0)
    names(myOutput) <- c("EAS", "AFR", "EUR")
    expect_equal(admixture(vcf, admix.ref)[[3]], myOutput, tolerance=1e-6)
})
