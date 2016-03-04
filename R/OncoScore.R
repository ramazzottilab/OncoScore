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
#' data(genes)
#' query = perform.web.query(genes[1:2])
#' compute.OncoScore(query)
#' 
#' @param data input data as result of the function perform.web.query
#' @param filter.threshold threshold to filter for a minimum number of citations for the genes
#' @param analysis.mode logaritmic scores to be computed, i.e., log10, log2, natural log or log5
#' @param cutoff.threshold threshold to be used to asses the oncogenes
#' @param file should I save the results to text files? 
#'
#' @return the computed OncoScores and the clusters for the genes
#' 
#' @export compute.OncoScore
#' 
compute.OncoScore <- function( data,
                               filter.threshold = NA,
                               analysis.mode = "Log2",
                               cutoff.threshold = 21.09,
                               file = NULL ) {
    
    # verify I was given a cutoff threshold
    if (is.na(cutoff.threshold)) {
        paste("Warning: the cutoff threshold to define the oncogenes was not provided.\n")
    }

    cat('### Processing data\n')
    
    data = cbind(data, 0)
    colnames(data)[3] = 'PercCit'
    
    # set the percentage of citations
    for (gene in rownames(data)) {
        if (is.na(data[gene, "CitationsGene"]) || data[gene, "CitationsGene"] <= 1) {
            data[gene, "PercCit"] = 0
        } else if (is.na(data[gene, "CitationsGeneInCancer"]) || data[gene, "CitationsGeneInCancer"] == -1) {
            data[gene,"PercCit"] = 0
        } else {
            data[gene, "PercCit"] = data[gene, "CitationsGeneInCancer"] * 100 / data[gene, "CitationsGene"]
        }
    }
  
    # compute the scores based on the frequencies
    scores = compute.frequencies.scores(data, filter.threshold, analysis.mode)
    
    # estimate the oncogenes
    results = estimate.oncogenes(scores, cutoff.threshold)
    
    # write the results to file
    if (!is.null(file) && is.character(file)) {
        write.table(results,file = file,
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)
    }
    
    cat('### Results:\n')
    for (gene in rownames(results)) {
        cat('\t', gene, '->', results[gene,'OncoScore'],'\n')
    }
    return(results)
}


#' perform the OncoScore time series analysis for a list of genes and data times
#' 
#' @title compute.OncoScore.TimeSeries
#'
#' @examples
#' data(genes)
#' data(timepoints)
#' query = perform.time.series.query(genes[1:2], timepoints[1:2])
#' compute.OncoScore.TimeSeries(query)
#' 
#' @param data input data as result of the function perform.time.series.query
#' @param filter.threshold threshold to filter for a minimum number of citations for the genes
#' @param analysis.mode logaritmic scores to be computed, i.e., log10, log2, natural log or log5
#' @param cutoff.threshold threshold to be used to asses the oncogenes
#' @param file should I save the results to text files? 
#'
#' @return the performed OncoScores time series analysis
#' 
#' @export compute.OncoScore.TimeSeries
#' 
compute.OncoScore.TimeSeries <- function( data,
                                          filter.threshold = NA,
                                          analysis.mode = "Log2",
                                          cutoff.threshold = 21.09,
                                          file = NULL ) {
    
    # structure where to save the results
    results = NULL
    
    # perform the analysis for each time point
    for (timepoint in names(data)) {
        cat("### Computing oncoscore for timepoint", timepoint, '\n')
        destination = NULL
        if (!is.null(file)) {
            destination = paste0(gsub("/", "_", timepoint), "_", file)
        }
        curr_time_result = compute.OncoScore(data[[timepoint]],
                                             filter.threshold = filter.threshold,
                                             analysis.mode = analysis.mode,
                                             cutoff.threshold = cutoff.threshold,
                                             file = destination)
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
    
    cat('### Computing frequencies scores \n')
    # filter for a minimum number of citations for the genes if requested
    if (!is.na(filter.threshold)) {
        data = data[data[, "CitationsGene"] >= filter.threshold, ]
    } else {
        filter.threshold = 0
    }

    data = cbind(data, matrix(0, ncol=4, nrow=nrow(data)))
    
    # structure where to save the results
    colnames(data)[4:7] = c("alpha",
                            "1/alpha",
                            "PercCit*1/alpha",
                            "OncoScore")

    
    # compute the scores
    if (analysis.mode == "Log2") {
        # log base 2
        for (gene in rownames(data)) {
            if (data[gene, "CitationsGene"] <= filter.threshold | is.na(data[gene, "CitationsGene"])) {
                data[gene, "alpha"]             = 0
                data[gene, "1/alpha"]           = 0
                data[gene, "PercCit*1/alpha"]   = 0
                data[gene, "OncoScore"]     = 0
            } else {
                data[gene, "alpha"]             = log(data[gene, "CitationsGene"], 2)
                data[gene, "1/alpha"]           = 1 / log(data[gene, "CitationsGene"], 2)
                data[gene, "PercCit*1/alpha"]   = data[gene, "PercCit"] * data[gene, "1/alpha"]
                data[gene, "OncoScore"]     = data[gene, "PercCit"] - data[gene, "PercCit*1/alpha"]
            }
        }
    } else if (analysis.mode == "Log") {
        # natural log
        for (gene in rownames(data)) {
            if (data[gene, "CitationsGene"] <= filter.threshold | is.na(data[gene, "CitationsGene"])) {
                data[gene, "alpha"]             = 0
                data[gene, "1/alpha"]           = 0
                data[gene, "PercCit*1/alpha"]   = 0
                data[gene, "OncoScore"]     = 0
            } else {
                data[gene, "alpha"]             = log(data[gene, "CitationsGene"])
                data[gene, "1/alpha"]           = 1 / log(data[gene, "CitationsGene"])
                data[gene, "PercCit*1/alpha"]   = data[gene, "PercCit"] * data[gene, "1/alpha"]
                data[gene, "OncoScore"]     = data[gene, "PercCit"] - data[gene, "PercCit*1/alpha"]
            }
        }
    } else if (analysis.mode == "Log5") {
        # log base 5
        for (gene in rownames(data)) {
            if (data[gene, "CitationsGene"] <= filter.threshold | is.na(data[gene, "CitationsGene"])) {
                data[gene, "alpha"]             = 0
                data[gene, "1/alpha"]           = 0
                data[gene, "PercCit*1/alpha"]   = 0
                data[gene, "OncoScore"]     = 0
            } else {
                data[gene, "alpha"]             = log(data[gene, "CitationsGene"], 5)
                data[gene, "1/alpha"]           = 1 / log(data[gene, "CitationsGene"], 5)
                data[gene, "PercCit*1/alpha"]   = data[gene, "PercCit"] * data[gene,"1/alpha"]
                data[gene, "OncoScore"]     = data[gene, "PercCit"] - data[gene, "PercCit*1/alpha"]
            }
        }
    }
    else if (analysis.mode == "Log10") {
    # log base 10
        for (gene in rownames(data)) {
            if (data[gene, "CitationsGene"] <= filter.threshold | is.na(data[gene, "CitationsGene"])) {
                data[gene, "alpha"]             = 0
                data[gene, "1/alpha"]           = 0
                data[gene, "PercCit*1/alpha"]   = 0
                data[gene, "OncoScore"]     = 0
            } else {
                data[gene, "alpha"]             = log10(data[gene, "CitationsGene"])
                data[gene, "1/alpha"]           = 1 / log10(data[gene, "CitationsGene"])
                data[gene, "PercCit*1/alpha"]   = data[gene, "PercCit"] * data[gene, "1/alpha"]
                data[gene, "OncoScore"]     = data[gene, "PercCit"] - data[gene, "PercCit*1/alpha"]
            }
        }
    }
    return(data)
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

    cat('### Estimating oncogenes\n')
    # structure where to save the results
    data = cbind(data, 0)
    colnames(data)[8] = 'Clustering'
    
    # assess the genes being oncogenes
    for (gene in rownames(data)) {
        if (!is.na(cutoff.threshold)) {
            if (data[gene, "OncoScore"] <= cutoff.threshold) {
                data[gene, "Clustering"] = 0
            } else {
                data[gene, "Clustering"] = 1
            }
        } else {
            data[gene, "Clustering"] = NA
        }
    }
    
    return(data)
}
