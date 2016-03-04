#### OncoScore
####
#### Copyright (c) 2016, Luca De Sano, Carlo Gambacorti Passerini, 
#### Rocco Piazza, Daniele Ramazzotti, Roberta Spinelli
#### 
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.

#' Description
#'
#' @docType package
#' @name TRONCO
NULL


#' compute the OncoScore for a list of genes
#' 
#' @title compute.OncoScore
#'
#' @examples
#' 
#' @param data input data as result of the function perform.web.query
#' @param filter.threshold threshold to filter for a minimum number of citations for the genes
#' @param analysis.mode logaritmic scores to be computed, i.e., log10, log2, natural log or log5
#' @param cutoff.threshold threshold to be used to asses the oncogenes
#' @param make.report should I save the results to text files? 
#' @param results.file file where to save the results
#'
#' @return the computed OncoScores and the clusters for the genes
#' 
#' @export compute.OncoScore
#' 
compute.OncoScore <- function( data,
                               filter.threshold = NA,
                               analysis.mode = "Log2",
                               cutoff.threshold = 21.09,
                               make.report = TRUE,
                               results.file = "Results_OncoScore.txt" ) {
    
    # verify I was given a cutoff threshold
    if (is.na(cutoff.threshold)) {
        writeLines("Warning: the cutoff threshold to define the oncogenes was not provided.")
    }
    
    # make a matrix of cancer vs non cancer frequencies per gene
    DataGenesCancer = data$cancer
    colnames(DataGenesCancer) = c("GeneCounter", "GeneID", "CitationsGeneinCancer")
    DataGenes = data$all
    colnames(DataGenes) = c("GeneCounter", "GeneID", "CitationsGene")
    DataMatrixGenes = merge(DataGenesCancer,
                            DataGenes,by.x = "GeneID",
                            by.y = "GeneID",
                            all.x = TRUE,
                            all.y = TRUE)
    DataMatrixGenes$GeneCounter = as.numeric(DataMatrixGenes$GeneCounter.y)
    DataMatrixGenes$PercCit = DataMatrixGenes$GeneCounter
    data = DataMatrixGenes[, c(-2, -4, -6)]
    
    # set the percentage of citations
    for (i in 1:dim(data)[1]) {
        if (is.na(data[i, "CitationsGene"]) | data[i, "CitationsGene"]<=1) {
            data[i, "PercCit"] = 0
        }
        else if (is.na(data[i, "CitationsGeneinCancer"]) | data[i, "CitationsGeneinCancer"] == -1) {
            data[i,"PercCit"] = 0
        }
        else {
            data[i, "PercCit"] = data[i, "CitationsGeneinCancer"] * 100 / data[i, "CitationsGene"]
        }
    }
    
    # compute the scores based on the frequencies
    scores = compute.frequencies.scores(data, filter.threshold, analysis.mode)
    
    # estimate the oncogenes
    results = estimate.oncogenes(scores, cutoff.threshold)
    print(results)
    
    # write the results to file
    if (make.report) {
        write.table(results,file = results.file,
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)
    }
    
    return(results)
}


#' perform the OncoScore time series analysis for a list of genes and data times
#' 
#' @title perform.OncoScore.TimeSeries
#'
#' @examples
#' 
#' @param data input data as result of the function perform.time.series.query
#' @param filter.threshold threshold to filter for a minimum number of citations for the genes
#' @param analysis.mode logaritmic scores to be computed, i.e., log10, log2, natural log or log5
#' @param cutoff.threshold threshold to be used to asses the oncogenes
#' @param make.report should I save the results to text files? 
#' @param results.file file where to save the results
#'
#' @return the performed OncoScores time series analysis
#' 
#' @export perform.OncoScore.TimeSeries
#' 
perform.OncoScore.TimeSeries <- function( data,
                                          filter.threshold = NA,
                                          analysis.mode = "Log2",
                                          cutoff.threshold = 21.09,
                                          make.report = TRUE,
                                          results.file = "Results_TimeSeries_OncoScore.txt" ) {
    
    # structure where to save the results
    results = NULL
    
    # perform the analysis for each time point
    for (timepoint in names(data)) {
        writeLines(paste0("Computing oncoscore for timepoint ", timepoint))
        curr_time_result = compute.OncoScore(data[[timepoint]],
                                             filter.threshold = filter.threshold,
                                             analysis.mode = analysis.mode,
                                             cutoff.threshold = cutoff.threshold,
                                             make.report = make.report,
                                             results.file = paste0(gsub("/", "_", timepoint), "_", results.file))
        results[[timepoint]] = curr_time_result
    }
    
    return(results)
}


#' compute the logaritmic scores based on the frequencies of the genes
#' 
#' @title compute.frequencies.scores
#'
#' @examples
#' 
#' @param data input data as result of the function perform.web.query
#' @param filter.threshold threshold to filter for a minimum number of citations for the genes
#' @param analysis.mode logaritmic scores to be computed, i.e., log10, log2, natural log or log5
#'
#' @return the computed scores
#' 
#' @export compute.frequencies.scores
#' 
compute.frequencies.scores <- function( data,
                                        filter.threshold = NA, 
                                        analysis.mode = "Log2" ) {
    
    # filter for a minimum number of citations for the genes if requested
    if (!is.na(filter.threshold)) {
        data = data[data[, "CitationsGene"] >= filter.threshold, ]
    } else {
        filter.threshold = 0
    }
    
    # compute the scores
    if (any(analysis.mode %in% "Log2")) {
        # log base 2
        for (i in 1:dim(data)[1]) {
            if (data[i, "CitationsGene"] <= filter.threshold | is.na(data[i, "CitationsGene"])) {
                data[i, "Log2CitationsGene"]            = 0
                data[i, "1/Log2CitationsGene"]          = 0
                data[i, "PercCit*1/Log2Citass"]         = 0
                data[i, "PercCit-PercCit*1/Log2Citass"] = 0
            } else {
                data[i, "Log2CitationsGene"]            = log(data[i, "CitationsGene"], 2)
                data[i, "1/Log2CitationsGene"]          = 1 / log(data[i, "CitationsGene"], 2)
                data[i, "PercCit*1/Log2Citass"]         = data[i, "PercCit"] * data[i, "1/Log2CitationsGene"]
                data[i, "PercCit-PercCit*1/Log2Citass"] = data[i, "PercCit"] - data[i, "PercCit*1/Log2Citass"]
            }
        }
    } else if (any(analysis.mode %in% "Log")) {
        # natural log
        for (i in 1:dim(data)[1]) {
            if (data[i, "CitationsGene"] <= filter.threshold | is.na(data[i, "CitationsGene"])) {
                data[i, "LogCitationsGene"]            = 0
                data[i, "1/LogCitationsGene"]          = 0
                data[i, "PercCit*1/LogCitass"]         = 0
                data[i, "PercCit-PercCit*1/LogCitass"] = 0
            } else {
                data[i, "LogCitationsGene"]            = log(data[i, "CitationsGene"])
                data[i, "1/LogCitationsGene"]          = 1 / log(data[i, "CitationsGene"])
                data[i, "PercCit*1/LogCitass"]         = data[i, "PercCit"] * data[i, "1/LogCitationsGene"]
                data[i, "PercCit-PercCit*1/LogCitass"] = data[i, "PercCit"] - data[i, "PercCit*1/LogCitass"]
            }
        }
    } else if (any(analysis.mode %in% "Log5")) {
        # log base 5
        for (i in 1:dim(data)[1]) {
            if (data[i, "CitationsGene"] <= filter.threshold | is.na(data[i, "CitationsGene"])) {
                data[i, "Log5CitationsGene"]            = 0
                data[i, "1/Log5CitationsGene"]          = 0
                data[i, "PercCit*1/Log5Citass"]         = 0
                data[i, "PercCit-PercCit*1/Log5Citass"] = 0
            } else {
                data[i, "Log5CitationsGene"]            = log(data[i, "CitationsGene"], 5)
                data[i, "1/Log5CitationsGene"]          = 1 / log(data[i, "CitationsGene"], 5)
                data[i, "PercCit*1/Log5Citass"]         = data[i, "PercCit"] * data[i,"1/Log5CitationsGene"]
                data[i, "PercCit-PercCit*1/Log5Citass"] = data[i, "PercCit"] - data[i, "PercCit*1/Log5Citass"]
            }
        }
    }
    else if (any(analysis.mode %in% "Log10")) {
    # log base 10
        for (i in 1:dim(data)[1]) {
            if (data[i, "CitationsGene"] <= filter.threshold | is.na(data[i, "CitationsGene"])) {
                data[i, "Log10CitationsGene"]            = 0
                data[i, "1/Log10CitationsGene"]          = 0
                data[i, "PercCit*1/Log10Citass"]         = 0
                data[i, "PercCit-PercCit*1/Log10Citass"] = 0
            } else {
                data[i, "Log10CitationsGene"]            = log10(data[i, "CitationsGene"])
                data[i, "1/Log10CitationsGene"]          = 1 / log10(data[i, "CitationsGene"])
                data[i, "PercCit*1/Log10Citass"]         = data[i, "PercCit"] * data[i, "1/Log10CitationsGene"]
                data[i, "PercCit-PercCit*1/Log10Citass"] = data[i, "PercCit"] - data[i, "PercCit*1/Log10Citass"]
            }
        }
    }
    
    # structure where to save the results
    results = data
    colnames(results) = c("GeneID",
                          "CitationsGeneinCancer",
                          "CitationsGene",
                          "PercCit",
                          "alpha",
                          "1/alpha",
                          "PercCit*1/alpha",
                          "PercCitWeigth")
    return(results)
}


#' estimate the oncoscore for the genes
#' 
#' @title estimate.oncogenes
#'
#' @examples
#' 
#' @param data input data as result of the function compute.frequencies.scores
#' @param cutoff.threshold threshold to be used to asses the oncogenes
#'
#' @return the computed scores and oncogenes
#' 
#' @export estimate.oncogenes
#' 
estimate.oncogenes <- function( data,
                                cutoff.threshold = 21.09 ) {
    
    # assess the genes being oncogenes
    for (i in 1:dim(data)[1]) {
        if (!is.na(cutoff.threshold)) {
            if (data[i, "PercCitWeigth"] <= cutoff.threshold) {
                data[i, "Clustering"] = 0
            } else {
                data[i, "Clustering"] = 1
            }
        } else {
            data[i, "Clustering"] = NA
        }
    }
    
    # structure where to save the results
    results = data
    colnames(results) = c("GeneID",
                          "CitationsGeneinCancer",
                          "CitationsGene",
                          "PercCit",
                          "alpha",
                          "1/alpha",
                          "PercCit*1/alpha",
                          "OncoScore","Clustering")
    return(results)
}
