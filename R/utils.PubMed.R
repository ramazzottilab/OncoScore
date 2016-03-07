#### OncoScore
####
#### Copyright (c) 2016, Luca De Sano, Carlo Gambacorti Passerini, 
#### Rocco Piazza, Daniele Ramazzotti, Roberta Spinelli
#### 
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' query PubMed for a list of genes
#' 
#' @title get.pubmed.driver.analysis
#' 
#' @param keywords TODO
#' @param gene TODO
#'
#' @return TODO
#'
get.pubmed.driver.analysis <- function(keywords,
                                       gene) {

    # setup the parameters for the queries to PubMed
    options("scipen" = 10)

    # setup the query to PubMed
    keywords = gsub(" ", "+", keywords)
    getURL = paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                    "esearch.fcgi?db=pubmed&term=",
                    keywords,
                    "&retmax=3000000",
                    "&rettype=FASTA",
                    "&tool=OncoScore")

    # perform the query to PubMed
    webget = try.scan(getURL)
    lc = get.list.from.xml(webget)

    cat("\tNumber of papers found in PubMed for",
        gene,
        "was:",
        lc,
        "\n")
    return(lc)
}

#' try to query the given URL
#' 
#' @title try.scan
#' 
#' @param getURL The given URL
#' 
#' @return TODO
#'
try.scan <- function(getURL) {


    # try the query
    thispage = try(scan(getURL,
                        what = "character",
                        sep = "\n",
                        quiet = TRUE),
                        silent = TRUE)

    retrycount = 0
    while (class(thispage) == "try-error" && retrycount < 10) {
        Sys.sleep(0.5)
        retrycount = retrycount + 1
        thispage = try(scan(getURL,
                            what = "character",
                            sep = "\n",
                            quiet = TRUE),
                       silent = TRUE)
    }
    if (class(thispage) == "try-error") {
        stop("Error: The PubMed website could not be reached.")
    }
    return(thispage)
}


#' process the result of the query
#' 
#' @title get.list.from.xml
#' 
#' @param webget TODO
#' 
#' @return TODO
#'
get.list.from.xml <- function(webget) {
    
    # get the 0 length returns
    get.zero <- function( webget ) {
        z = length(grep("<Count>0</Count>", webget))
        return(z)
    }
    
    # get the count of list items
    get.count.of.list <- function( webget ) {
        a = as.numeric(gsub("^[[:print:]]*<Count>([[:digit:]]+)[[:print:]]*$",
                            "\\1",
                            webget[grep("<Count>", webget)][1]))
        return(a)
    }
    
    # get the phrases not found in the text
    get.phrase.not.found <- function( webget ) {
        b = as.numeric(grep("PhraseNotFound", webget))
        return(b)
    }
    
    # process the result of the query
    webget = gsub("\t", "", webget)
    webget = gsub("<IdList>", "\t<IdList>", webget)
    webget = gsub("</IdList>", "\t</IdList>", webget)
    webget = unlist(strsplit(webget, "\t"))
    s1 = grep("^<IdList>$", webget)[1]
    s2 = grep("^</IdList>", webget)[1]
    
    if ((is.na(s1) || is.na(s2)) 
        && (get.count.of.list(webget) == 0) 
        && (length(get.phrase.not.found(webget)) != 0)) {
        lc = -1
    } else if ((is.na(s1) || is.na(s2)) 
        && (get.zero(webget) != 0) 
        && (length(get.phrase.not.found(webget)) != 0)) {
        lc = -1
    } else if ((is.na(s1) || is.na(s2)) 
        && (get.zero(webget) != 0)) {
        lc = 0
    } else if ((is.na(s1) || is.na(s2)) 
        && (get.count.of.list(webget) == 0)) {
        lc = get.count.of.list(webget)
    } else if (is.na(s1) || is.na(s2)) {
        stop("error: parsing XML problem.")
    } else {
        ListItems = gsub("[[:print:]]*<Id>([[:print:]]*)</Id>$",
                         "\\1",
                         webget[(s1 + 1):(s2 - 1)])
        ListItems = ListItems[ListItems != ""]
        lc = get.count.of.list(webget)
        if (lc > 3e+06) {
            cat("Warning: Number of items was greater than expected. PARTIAL RESULTS USED [MaxRet needs to be increased].\n")
            cat("Actual Length:", lc, "\n")
            cat("Returned Length:", MaxRet, "\n")
        }

    }
    
    if (!is.na(lc) 
        && lc != 0 
        && lc != -1
        && lc != length(ListItems)) {
            stop("Warning: Number of items was different than the number of references.")
    }
    return(lc)
}
