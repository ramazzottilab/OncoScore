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
#' @param ... additional parameter to pass to the barplot function
#' 
#' @export plot.oncoscore
#' 
plot.oncoscore <- function(x,
                           gene.number = 5,
                           ...) {

    x = x[, "OncoScore"]
    x = sort(x, decreasing = TRUE)
    if (length(x) > gene.number) {
        x = x[1:gene.number]
    }
    barplot(sort(x), 
            horiz = TRUE,
            xlab = 'OncoScore',
            ylab = 'Genes',
            ...)
}