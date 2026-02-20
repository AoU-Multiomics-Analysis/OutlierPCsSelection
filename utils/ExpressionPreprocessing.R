LoadGCTFile <- function(PathGCT) {
message('Loading GCT File')
CountData <- fread(pathGCT,skip = 2,header = TRUE)
CountData
}


FilterCountData <- function(CountData,CountThresh = 6,PropSamples = 0.2) {
message('Transposing counts')
CountDataTransposed <- CountData %>%
    dplyr::select(-Description) %>% 
    column_to_rownames('Name') %>%
    t() %>% 
    data.frame()

message('Filtering Genes by read count')
CountDataFiltered <- CountDataTransposed %>%
        dplyr::select(where(~ mean(.x > 6) >= 0.2)) %>% 
        t() %>% 
        data.frame()
CountDataFiltered
}

NormalizeCountsCPMs <- function(CountData) {
message('Calculating normalization factors')
DataEdgeR <- edgeR::DGEList(CountData)
DataEdgeR <- edgeR::calcNormFactors(DataEdgeR)

message('Computing CPMs')
DataLogCPM <- edgeR::cpm(DataEdgeR, log=TRUE) %>% data.frame()
SamplesByGenes <- DataLogCPM %>% t() %>% data.frame()
rownames(SamplesByGenes) <- str_remove(rownames(SamplesByGenes),'^X')
SamplesByGenes
}
