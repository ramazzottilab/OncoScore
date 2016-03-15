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
#' result = compute.oncoscore(query)
#' plot.oncoscore(result)
#' 
#' @param x input data as result of the function compute.OncoScore
#' @param gene.number number of genes to print
#' @param main the title
#' @param xlab description of x asix (defaul score)
#' @param ylab description of y asix (defaul genes)
#' @param file where to save the plot
#' @param ... additional parameter to pass to the barplot function
#' 
#' @return A plot
#' 
#' @export plot.oncoscore
#' @importFrom graphics barplot
#' @importFrom grDevices dev.copy2pdf
#' 
plot.oncoscore <- function(x,
                           gene.number = 5,
                           main = 'OncoScore',
                           xlab = 'score',
                           ylab = 'genes',
                           file = NA,
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

    if (!is.na(file)) {
        dev.copy2pdf(file = file)
    }
}


#' plot the OncoScore for a list of genes
#' 
#' @title plot.oncoscore.timeseries
#'
#' @examples
#' data(query.timepoints)
#' result = compute.oncoscore.timeseries(query.timepoints)
#' plot.oncoscore.timeseries(result)
#' 
#' @param x input data as result of the function compute.OncoScore
#' @param gene.number number of genes to print
#' @param incremental display the OncoScore increment
#' @param relative dispaly the incrementa as relative value
#' @param main the title
#' @param xlab description of x asix (defaul score)
#' @param ylab description of y asix (defaul genes)
#' @param legend.pos Position of the legend
#' @param file where to save the plot
#' @param ... additional parameter to pass to the lines function
#' 
#' @return A plot
#' 
#' @export plot.oncoscore.timeseries
#' @importFrom grDevices dev.copy2pdf rainbow
#' @importFrom graphics axis legend lines plot
#' 
plot.oncoscore.timeseries <- function(x,
                                      gene.number = 5,
                                      incremental = FALSE,
                                      relative = FALSE,
                                      main = 'OncoScore',
                                      xlab = 'timepoints',
                                      ylab = 'score',
                                      legend.pos = 'top',
                                      file = NA,
                                      ...) {


    sorted = sort(x[[length(x)]][, "OncoScore"], decreasing = TRUE)
    if (length(sorted) > gene.number) {
        sorted = sorted[1:gene.number]
    }
    sorted = names(sorted)

    v = NULL

    for (timepoint in names(x)) {
        data = x[[timepoint]][sorted, "OncoScore", drop=FALSE]
        v = cbind(v, data)
    }

    colnames(v) = names(x)
    rownames(v) = sorted

    if (incremental) {
        values = v
        v = matrix(0, nrow = nrow(values))
        rownames(v) = sorted
        for(i in 2:ncol(values)) {
            if (relative) {
              v = cbind(v, values[, i, drop = FALSE] / values[, i - 1, drop = FALSE])
            } else {
              v = cbind(v, values[, i, drop = FALSE] - values[, i - 1, drop = FALSE])
            }
        }
    }

    minval = min(v) - 3
    maxval = max(v) + 3

    plot(c(1, length(names(x))),
         c(minval, maxval),
         type="n", 
         main = main,
         xlab=xlab,
         ylab=ylab,
         axes = FALSE) 

    axis(1, at=1:length(names(x)), labels=names(x))
    axis(2, at=seq(round(minval), round(maxval), by=5), las=1)


    color = rainbow(nrow(v))
    names(color) = rownames(v)

    for (gene in rownames(v)) {
        lines(v[gene,], col=color[gene], lwd=2)
    }

    legend(legend.pos, rownames(v), col=color,  lty=1, lwd=2, cex=0.8, horiz = TRUE)

    if (!is.na(file)) {
        dev.copy2pdf(file = file)
    }
}
