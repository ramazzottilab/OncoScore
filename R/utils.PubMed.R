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
#' @title GetPubMedDriverAnalysis
#' 
#' @param keywords TODO
#' @param gene TODO 
#' @param showurl TODO
#' @param MaxRet TODO
#' @param sme TODO
#' @param smt TODO
#'
#' @return TODO
#'
GetPubMedDriverAnalysis <- function(keywords,
                                    gene,
                                    showurl = FALSE,
                                    MaxRet = 3e+06,
                                    sme = FALSE,
                                    smt = FALSE) {

    # setup the parameters for the queries to PubMed
    options("scipen" = 10)
    URLdef = ncbi2r.driver.options()

    # setup the query to PubMed
    keywords = gsub(" ", "+", keywords)
    getURL = paste(URLdef$front,
                   "esearch.fcgi?db=pubmed&term=",
                   keywords,
                   "&retmax=",
                   MaxRet,
                   "&rettype=FASTA",
                   URLdef$back,
                   sep = "")

    # perform the query to PubMed
    webget = get.driver.file(getURL,
                             showurl = showurl,
                             clean = FALSE)
    lc = getListFromXMLdriver(webget,
                              sme = sme,
                              smt = smt)

    print(paste("Number of papers found in PubMed for",
                gene,
                "was:",
                lc,
                sep = " "))
    return(lc)
}


#' setup the parameters for the queries to PubMed to a global variable
#' 
#' @title ncbi2r.driver.options
#' 
#' @return TODO
#'
ncbi2r.driver.options <- function() {

    if (!exists("ncbi2r.options")) {
        ncbi2r.options <<- ncbi2r.driver.options.default()
    }
    return(ncbi2r.options)
}


#' setup the default parameters for the queries to PubMed
#' 
#' @title ncbi2r.driver.options.default
#' 
#' @return TODO
#'
ncbi2r.driver.options.default <- function() {

    baseurl     = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    tool        = "NCBI2R"
    email       = "ncbi2r@gmail.com"
    tailurl     = paste("&tool=", tool, "&email=", email, sep = "")
    TimeStampA  = Sys.time() - 3
    
    return(list(front = baseurl,
                back = tailurl,
                TimeStampA = TimeStampA))
}


#' get the URL webpage
#' 
#' @title get.driver.file
#' 
#' @param getURL TODO
#' @param showurl TODO
#' @param sep TODO
#' @param quiet TODO
#' @param clean TODO
#' @param verbose TODO
#' 
#' @return TODO
#'
get.driver.file <- function(getURL,
                            showurl = FALSE,
                            sep = "\n",
                            quiet = TRUE,
                            clean = FALSE,
                            verbose = FALSE) {

    # write the URL to output if requered
    if (showurl) {
        writeLines(getURL)
    }

    # try the query to get the URL webpage
    webget = suppressWarnings(tryScan(getURL,
                                      sep = sep,
                                      quiet = quiet,
                                      verbose = verbose))

    if (clean) {
        webget = clean.xml(webget)
    }

    return(webget)
}


#' try to query the given URL
#' 
#' @title tryScan
#' 
#' @param getURL TODO
#' @param sep TODO
#' @param quiet TODO
#' @param retry TODO
#' @param error TODO
#' @param verbose TODO
#' 
#' @return TODO
#'
tryScan <- function(getURL,
                    sep = "\n",
                    quiet = TRUE,
                    retry = 20,
                    error = TRUE,
                    verbose = FALSE) {

    # write the URL to output if requered
    if (verbose) {
        writeLines(getURL)
    }

    # try the query
    thispage = try(scan(getURL,
                        what = "character",
                        sep = sep,
                        quiet = quiet),
                   silent = TRUE)

    retrycount = 0
    while (class(thispage) == "try-error" & retrycount < retry) {
        DelayDriver(0.5)
        retrycount = retrycount + 1
        if (verbose) {
            writeLines(as.character(retrycount))
        }
        thispage = try(scan(getURL,
                            what = "character",
                            sep = sep,
                            quiet = quiet),
                       silent = TRUE)
    }
    if (class(thispage) == "try-error" && error) {
        stop("NCBI2R error: The PubMed website could not be reached.")
    }
    if (class(thispage) == "try-error" && !error) {
        thispage <- "Page not found"
    }

    return(thispage)
}

#' repeat the query after a given time
#' 
#' @title DelayDriver
#' 
#' @param seconds TODO
#' 
#' @return TODO
#'
DelayDriver <- function(seconds) {
    
    ncbi2r.driver.options()
    ncbi2r.options$TimeStampA <<- Sys.time()
    while (Sys.time() < ncbi2r.options$TimeStampA + seconds) {
        OnlyForDelay = 1
    }
    
}


#' clean the webpage returned by the query
#' 
#' @title DelayDriver
#' 
#' @param webget TODO
#' 
#' @return TODO
#'
clean.xml <- function(webget) {
    
    webget = gsub("&lt;", "<", webget)
    webget = gsub("&gt;", ">", webget)
    webget = gsub("                <", "<", webget)
    webget = gsub("                <", "<", webget)
    webget = gsub("              <", "<", webget)
    webget = gsub("            <", "<", webget)
    webget = gsub("          <", "<", webget)
    webget = gsub("        <", "<", webget)
    webget = gsub("      <", "<", webget)
    webget = gsub("    <", "<", webget)
    webget = gsub("  <", "<", webget)
    webget = gsub(" <", "<", webget)
    
    return(webget)
}


#' process the result of the query
#' 
#' @title DelayDriver
#' 
#' @param webget TODO
#' @param smt TODO
#' @param sme TODO
#' @param MaxRet TODO
#' @param return.data TODO
#' 
#' @return TODO
#'
getListFromXMLdriver <- function(webget,
                                 smt = FALSE,
                                 sme = FALSE, 
                                 MaxRet = 3e+06,
                                 return.data = FALSE) {
    
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
        stop("NCBI2R error: parsing XML problem.")
    } else {
        ListItems = gsub("[[:print:]]*<Id>([[:print:]]*)</Id>$",
                         "\\1",
                         webget[(s1 + 1):(s2 - 1)])
        ListItems = ListItems[ListItems != ""]
        dummy = showMessages(webget, smt = smt, sme = sme)
        if (class(dummy$errors) == "data.frame") {
            writeLines("Error using search query.")
            print(dummy$errors)
            if (class(dummy$warnings) == "data.frame") {
                writeLines("Additional warnings were also found.")
                print(dummy$warnings)
            }
            lc = NA
        } else {
            lc = get.count.of.list(webget)
            if (lc > MaxRet) {
                writeLines("Warning: Number of items was greater than expected. PARTIAL RESULTS USED [MaxRet needs to be increased].")
                writeLines(paste("Actual Length:", lc))
                writeLines(paste("Returned Length:", MaxRet))
            }
        }
    }
    
    if (!is.na(lc) && lc != 0 && lc != -1) {
        if (lc == length(ListItems)) {
            if (return.data == TRUE) {
                return(list(ListItems = ListItems,
                            errors = dummy$errors,
                            warnings = dummy$warnings))
            } else {
                return(lc)
            }
        } else {
            writeLines("Warning: Number of items was different than the number of references.")
            stop("")
        }
    } else {
        return(lc)
    }
}


#' show error and warnings messages to standard output
#' 
#' @title showMessages
#' 
#' @param webget TODO
#' @param smt TODO
#' @param sme TODO
#' 
#' @return TODO
#'
showMessages <- function(webget,
                         smt = FALSE,
                         sme = FALSE) {
    
    cleanedset = parse.items(webget)

    if (length(grep("<QueryTranslation>", webget)) == 0) {
        answer = "******No Query Translation available******"
        rem = webget
    } else {
        answer = gsub("^[[:print:]]*<QueryTranslation>([[:print:]]*)</QueryTranslation>[[:print:]]*",
                      "\\1",
                      webget[grep("<QueryTranslation>", webget)])
        rem = gsub("^[[:print:]]*<QueryTranslation>([[:print:]]*)</QueryTranslation>([[:print:]]*)",
                   "\\2",
                   webget[grep("<QueryTranslation>", webget)])
    }
    s1 = "No errors"
    s2 = "No warnings"
    if (smt) {
        print(paste("QueryTranslation was:", answer))
    }
    if (length(grep("<ErrorList>", rem)) == 1 || length(grep("<WarningList>", rem)) == 1) {
        Y = unlist(strsplit(rem, "<"))
        if (length(grep("<ErrorList>", rem)) == 1) {
            subtxt = Y[(grep("ErrorList>", Y)[1] + 1):(grep("ErrorList>", Y)[2] - 1)]
            s1 = as.data.frame(matrix(subtxt, ncol = 2, byrow = TRUE), stringsAsFactors = FALSE)
            s1$tag = gsub("/|>", "", s1$V2)
            s1$text = gsub("^[[:print:]]*>([[:print:]]*)$", "\\1", s1$V1)
            s1 = s1[, c("tag", "text")]
            if (sme) { 
                writeLines(paste("Error On Query (", s1$tag, "): ", s1$text, sep = ""))
            }
        }
        
        if (length(grep("<WarningList>", rem)) == 1) {
            subtxt = Y[(grep("WarningList>", Y)[1] + 1):(grep("WarningList>", Y)[2] - 1)]
            s2 = as.data.frame(matrix(subtxt, ncol = 2, byrow = TRUE), stringsAsFactors = FALSE)
            s2$tag = gsub("/|>", "", s2$V2)
            s2$text = gsub("^[[:print:]]*>([[:print:]]*)$", "\\1", s2$V1)
            s2 = s2[, c("tag", "text")]
            if (sme) { 
                writeLines(paste("Warning On Query (", s2$tag, "): ", s2$text, sep = ""))
            }
        }
    }
    return(list(errors = s1, warnings = s2))
}


#' utility function to parse items
#' 
#' @title parse.items
#' 
#' @param dataset TODO
#' 
#' @return TODO
#'
parse.items <- function( dataset ) {
    
    while (length(splitfirst(dataset[length(dataset)], "  ")) == 2) {
        a = splitfirst(dataset[length(dataset)], "  ")
        counter = length(dataset)
        while (substr(a[2], 1, 1) == " ") {
            a[2] = substr(a[2], 2, nchar(a[2]))
        }
        dataset[counter + 1] = a[2]
        dataset[counter] = a[1]
        counter = counter + 1
    }
    
    return(dataset)
}

#' utility function to split items
#' 
#' @title splitfirst
#' 
#' @param lineoftext TODO
#' @param searchfor TODO
#' @param characterposition TODO
#' 
#' @return TODO
#'
splitfirst <- function( lineoftext,
                        searchfor = " ",
                        characterposition = 1 ) {
    
    splitnow = FALSE
    content = rep("", 1000)
    word = 1
    while (characterposition <= nchar(lineoftext) && !splitnow) {
        if (substr(lineoftext, characterposition, characterposition + nchar(searchfor) - 1) == searchfor) {
            splitnow = TRUE
            word = word + 1
            content[word - 1] = substr(lineoftext, 1, characterposition - 1)
            content[word] = substr(lineoftext, characterposition + nchar(searchfor), nchar(lineoftext))
        }
        characterposition = characterposition + 1
    }
    
    return(content[content!=""])
}
