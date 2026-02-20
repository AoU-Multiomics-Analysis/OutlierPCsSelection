library(tidyverse)
library(data.table)
library(magrittr)
library(optparse)
library(data.table)
library(rtracklayer)
library(RNOmni)
library(edgeR)

source('/opt/r/lib/ComputeBetasPCs.R')
source('/opt/r/lib/ExpressionPreprocessing.R')

Optlist <- function() {
option_list <- list(
  optparse::make_option(c("--GCTFile"), default=NULL,
    help="Minor allele frequency filter applied to genotype matrix", metavar="type",type="numeric"),
)
option_list
}


####### LOAD DATA, NORMALIZE AND COMPUTE PCs########

message('Loading and filtering count data ')
CountData <- LoadGCTFile(opt$GCtFile) %>% FilterCountData()
SampleByGenesCPM <- NormalizeCountsCPMs

message('Computing Principal components')
PCARes <- ComputePCs(SampleByGenesCPM)
Rotated <- PCARes$rotated %>% data.frame()


######## CALCUALTE Z SCORES ACROSS PCS #############
message('Computing Beta values for PCs')
BetaPCs <- ComputePCBetaValues(Rotated,SampleByGenesCPM)

PC_list <- seq(0,8900,by = 500)
message('Computing z scores across PCs')
ZscoreList <- ZScores_Kgrid(Rotated, 
                             BetaPCs,
                             SamplesByGenes, 
                             K_grid = PC_list, 
                             gene_block = 2000
                            )
ZscoreList %>% saveRDS('ZscoresAcrossPCs.rds')

