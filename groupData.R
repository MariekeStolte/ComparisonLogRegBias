# groupData - Helper function that groups data
# Input: data - (data.frame) Ungrouped data (All columns for covariates
#                 should be categorical)
#        y.name - (character) Name of the binary target variable (default: y)
#        cases - (numeric or character) value that cases have in target variable
#                   (default: 1)
# Output: (data.frame) Grouped data set with 1st column Cases, 2nd column Controls
#
groupData <- function(data, y.name = "y", cases = 1) {
  ind <- which(colnames(data) == y.name)
  combs <- apply(data[, -ind, drop = FALSE], 1, function(x) paste(x, collapse = ""))
  tab <- table(combs, factor(data[, ind], levels = 0:1))
  x <- do.call(rbind, lapply(strsplit(rownames(tab), ""), as.numeric))
  grouped <- data.frame(Cases = tab[, 2], Controls = tab[, 1], x)
  colnames(grouped)[-(1:2)] <- colnames(data)[-ind]
  rownames(grouped) <- NULL
  return(grouped)
}


