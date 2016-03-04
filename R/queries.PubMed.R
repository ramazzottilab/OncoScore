#### OncoScore
####
#### Copyright (c) 2016, Luca De Sano, Carlo Gambacorti Passerini, 
#### Rocco Piazza, Daniele Ramazzotti, Roberta Spinelli
#### 
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' perforn the query to PubMed
#' 
#' @title perform.web.query
#'
#' @examples
#' 
#' @param list.of.genes TODO
#' @param gene.num.limit TOOD
#'
#' @return TODO
#' 
#' @export perform.web.query
#' 
perform.web.query <- function( list.of.genes,
                               gene.num.limit = 100 ) {
    
    # perform the analysis
    cat("### Starting the queries for the selected genes.\n")
    
    # if the data are save in a file, read it
    if (is.character(list.of.genes) 
        && file.exists(list.of.genes)) {
        cat("### Reading the list of genes from file: ", paste(list.of.genes, collapse=' '), '\n')
        list.of.genes = read.table(file = list.of.genes,
                                   header = FALSE,
                                   row.names = NULL,
                                   check.names = FALSE,
                                   stringsAsFactors = FALSE)
    }

    if (!is.vector(list.of.genes)) {
        list.of.genes = unlist(list.of.genes)
        names(list.of.genes) = NULL
    }

    # set the genes number
    list.of.genes = unique(list.of.genes)
    genes.number = length(list.of.genes)
    
    if (genes.number > gene.num.limit) {
        stop("Too many genes to query! Please reduce the list and try again.")
    }
    
    # perform the query for the cancer topics
    search.fields = paste("[All Fields] AND ((lymphoma[MeSH Terms] OR lymphoma[All Fields])",
        "OR (lymphoma[MeSH Terms] OR lymphoma[All Fields] OR lymphomas[All Fields])",
        "OR (neoplasms[MeSH Terms] OR neoplasms[All Fields] OR cancer[All Fields])", 
        "OR (tumour[All Fields] OR neoplasms[MeSH Terms] OR neoplasms[All Fields]",
        "OR tumor[All Fields]) OR (neoplasms[MeSH Terms] OR neoplasms[All Fields]", 
        "OR neoplasm[All Fields]) OR (neoplasms[MeSH Terms] OR neoplasms[All Fields]", 
        "OR malignancy[All Fields]) OR (leukaemia[All Fields] OR leukemia[MeSH Terms]", 
        "OR leukemia[All Fields]) OR (neoplasms[MeSH Terms] OR neoplasms[All Fields]", 
        "OR cancers[All Fields]) OR (tumours[All Fields] OR neoplasms[MeSH Terms]", 
        "OR neoplasms[All Fields] OR tumors[All Fields]) OR (neoplasms[MeSH Terms]", 
        "OR neoplasms[All Fields] OR malignancies[All Fields]) OR (leukaemias[All Fields]", 
        "OR leukemia[MeSH Terms] OR leukemia[All Fields] OR leukemias[All Fields]))")
    

    cat("\n### Performing queries for cancer literature \n")
    ans = NULL

    for (gene in list.of.genes) {
        lc = GetPubMedDriverAnalysis(paste0(gene, search.fields),
                                     gene = gene)
        if (lc == -1) {
            warning(gene, ' not found on PubMed\n')
        }
        ans = rbind(ans, lc)
    }

    rownames(ans) = list.of.genes
    colnames(ans) = 'CitationsGeneInCancer'
    pubMedPanelGenes.Cancer = ans
    
    # perform the query for all the topics
    search.fields = "[All Fields]"
    
    cat("\n### Performing queries for all the literature \n")
    ans = NULL

    for (gene in list.of.genes) {
        lc = GetPubMedDriverAnalysis(paste0(gene, search.fields),
                                     gene = gene)
        if (lc == -1) {
            warning(gene, ' not found on PubMed\n')
        }
        ans = rbind(ans, lc)
    }
    
    rownames(ans) = list.of.genes
    colnames(ans) = 'CitationsGene'
    pubMedPanelGenes.All = ans
    
    pubMedPanelGenes = cbind(pubMedPanelGenes.All, pubMedPanelGenes.Cancer)
    return(pubMedPanelGenes)
}


#' perforn the query to PubMed for the time series analysis
#' 
#' @title perform.time.series.query
#'
#' @examples
#' 
#' @param list.of.genes TODO
#' @param list.of.datatimes TODO 
#' @param gene.num.limit TODO
#' @param timepoints.limit TODO
#'
#' @return TODO
#' 
#' @export perform.time.series.query
#'
perform.time.series.query <- function( list.of.genes,
                                       list.of.datatimes,
                                       gene.num.limit = 100,
                                       timepoints.limit = 10 ) {
    
    # perform the analysis
    writeLines("### Starting the queries for the selected genes and timepoints.")
    
    # if the data are save in a file, read it
    if (!(is.data.frame(list.of.genes) || is.matrix(list.of.genes)) && is.character(list.of.genes)) {
        writeLines(paste0("### Reading the list of genes from file: ",list.of.genes))
        list.of.genes = read.table(file = list.of.genes,
                                   header = FALSE,
                                   row.names = NULL,
                                   check.names = FALSE,
                                   stringsAsFactors = FALSE)
    }
    if (!(is.data.frame(list.of.datatimes) || is.matrix(list.of.datatimes)) && is.character(list.of.datatimes)) {
        writeLines(paste0("### Reading the list of timepoints from file: ",list.of.datatimes))
        list.of.datatimes = read.table(file = list.of.datatimes,
                                       header = FALSE,
                                       row.names = NULL,
                                       check.names = FALSE,
                                       stringsAsFactors = FALSE)
    }
    
    # set the genes number
    list.of.genes = unique(list.of.genes)
    genes.number = dim(list.of.genes)[1]
    list.of.datatimes = unique(list.of.datatimes)
    all.times = dim(list.of.datatimes)[1]
    
    if (genes.number > gene.num.limit) {
        stop("Too many genes to query! Please reduce the list and try again.")
    }
    
    if (all.times > timepoints.limit) {
        stop("Too timepoints to be considered! Please reduce the list and try again.")
    }
    
    # repeat the query for all the time points
    pubMedPanelGenes = NULL
    for(l in 1: all.times) {
        writeLines(paste0("### Quering PubMed for timepoint ",
                          list.of.datatimes[l, 1]))
    
        # perform the query for the cancer topics
        search.fields = paste0("[All Fields] AND ((lymphoma[MeSH Terms] OR lymphoma[All Fields]) 
            OR (lymphoma[MeSH Terms] OR lymphoma[All Fields] OR lymphomas[All Fields]) 
            OR (neoplasms[MeSH Terms] OR neoplasms[All Fields] OR cancer[All Fields]) 
            OR (tumour[All Fields] OR neoplasms[MeSH Terms] OR neoplasms[All Fields] 
            OR tumor[All Fields]) OR (neoplasms[MeSH Terms] OR neoplasms[All Fields] 
            OR neoplasm[All Fields]) OR (neoplasms[MeSH Terms] OR neoplasms[All Fields] 
            OR malignancy[All Fields]) OR (leukaemia[All Fields] OR leukemia[MeSH Terms] 
            OR leukemia[All Fields]) OR (neoplasms[MeSH Terms] OR neoplasms[All Fields] 
            OR cancers[All Fields]) OR (tumours[All Fields] OR neoplasms[MeSH Terms] 
            OR neoplasms[All Fields] OR tumors[All Fields]) OR (neoplasms[MeSH Terms] 
            OR neoplasms[All Fields] OR malignancies[All Fields]) OR (leukaemias[All Fields] 
            OR leukemia[MeSH Terms] OR leukemia[All Fields] OR leukemias[All Fields])) AND ",
            paste0("(0001/01/01[PDAT] : ",list.of.datatimes[l, 1], "[PDAT])"))
        
        ans = NULL
        for (gene.counter in 1:genes.number) {
            lc = GetPubMedDriverAnalysis(paste(list.of.genes[gene.counter, 1],
                                               search.fields,
                                               sep = ' '),
                                         gene = list.of.genes[gene.counter, 1])
            rowans = data.frame(gene.counter,
                                list.of.genes[gene.counter, 1],
                                lc,
                                stringsAsFactors = FALSE)
            ans = rbind(ans, rowans)
        }
        
        pubMedPanelGenes.Cancer = ans
        
        # perform the query for all the topics
        search.fields = paste0("[All Fields]",
                               paste0("(0001/01/01[PDAT] : ",
                               list.of.datatimes[l, 1],
                               "[PDAT])"))
        
        ans = NULL
        for(gene.counter in 1:genes.number) {
                lc = GetPubMedDriverAnalysis(paste(list.of.genes[gene.counter, 1],
                                                   search.fields,
                                                   sep = ' '),
                                             gene = list.of.genes[gene.counter, 1])
                rowans = data.frame(gene.counter,
                                    list.of.genes[gene.counter, 1],
                                    lc,
                                    stringsAsFactors = FALSE)
                ans = rbind(ans,rowans)
            }
        
        pubMedPanelGenes.All = ans
        pubMedPanelGenes[[list.of.datatimes[l, 1]]] = list(cancer = pubMedPanelGenes.Cancer,
                                                           all = pubMedPanelGenes.All)
        
    }
    return(pubMedPanelGenes)
}
