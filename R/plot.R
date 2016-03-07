#### OncoScore
####
#### Copyright (c) 2016, Luca De Sano, Carlo Gambacorti Passerini, 
#### Rocco Piazza, Daniele Ramazzotti, Roberta Spinelli
#### 
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' plot the OncoScore for a list of genes
#' 
#' @title plot.oncoscore
#'
#' @examples
#' data(query)
#' result = compute.OncoScore(query)
#' plot.oncoscore(result)
#' 
#' @param x input data as result of the function compute.OncoScore
#' @param gene.number number of genes to print
#' @param main the title
#' @param xlab description of x asix (defaul score)
#' @param ylab description of y asix (defaul genes)
#' @param ... additional parameter to pass to the barplot function
#' 
#' @export plot.oncoscore
#' 
plot.oncoscore <- function(x,
                           gene.number = 5,
                           main = 'OncoScore',
                           xlab = 'score',
                           ylab = 'genes',
                           ...) {

    x = x[, "OncoScore"]
    x = sort(x, decreasing = TRUE)
    if (length(x) > gene.number) {
        x = x[1:gene.number]
    }
    barplot(sort(x), 
            horiz = TRUE,
            main = main,
            xlab = xlab,
            ylab = ylab,
            ...)
}


#' plot the OncoScore for a list of genes
#' 
#' @title plot.oncoscore.timeseries
#'
#' @examples
#' data(query.timepoints)
#' result = compute.OncoScore.TimeSeries(query.timepoints)
#' plot.oncoscore.timeseries(result)
#' 
#' @param x input data as result of the function compute.OncoScore
#' @param gene.number number of genes to print
#' @param main the title
#' @param xlab description of x asix (defaul score)
#' @param ylab description of y asix (defaul genes)
#' @param ... additional parameter to pass to the lines function
#' 
#' @export plot.oncoscore.timeseries
#' 
plot.oncoscore.timeseries <- function(x,
                                      gene.number = 5,
                                      main = 'OncoScore',
                                      xlab = 'score',
                                      ylab = 'genes',
                                      ...) {

    print(x)

    sorted = sort(x[[length(x)]][, "OncoScore"], decreasing = TRUE)
    if (length(sorted) > gene.number) {
        sorted = sorted[1:gene.number]
    }
    sorted = names(sorted)

    minval = min(sapply(x, function(v){min(v[sorted,'OncoScore'])})) - 3
    maxval = max(sapply(x, function(v){max(v[sorted,'OncoScore'])})) + 3

    plot(c(1, length(names(x))),
         c(minval, maxval),
         type="n", 
         main = main,
         xlab=xlab,
         ylab=ylab,
         axes = FALSE) 

    axis(1, at=1:length(names(x)), lab=names(x))
    axis(2, at=seq(round(minval), round(maxval), by=5), las=1)

    v = NULL


    for (timepoint in names(x)) {
        data = x[[timepoint]][sorted, "OncoScore", drop=FALSE]
        v = cbind(v, data)
        print(data)
    }

    colnames(v) = names(x)
    rownames(x) = names(sorted)

    print(v)

    color = rainbow(nrow(v))
    names(color) = rownames(v)

    for (gene in rownames(v)) {
        lines(v[gene,], col=color[gene], lwd=2)
    }

    legend(1, maxval, rownames(v), col=color,  lty=1, lwd=2, cex=0.8, horiz = TRUE)

}