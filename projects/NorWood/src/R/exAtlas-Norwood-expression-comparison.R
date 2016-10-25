#' ---
#' title: "exAtlas and NorWood expressin profiles"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---

#' # Library and helper files
library(RColorBrewer)
source("~/Git/UPSCb/src/R/plot.multidensity.R")
pal <- brewer.pal(8,"Dark2")

#' # Data
exAtlas <- read.csv("/mnt/picea/projects/spruce/genome/HTSeq/Pabies-01_variance-stabilized-normalized-count-table.csv.gz",
                    row.names=1)
norWood <- read.delim("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/expression/vst.raw.expression_2016_07_06.txt",
                      header=FALSE)

#' # Plot
plot.multidensity(lapply(1:ncol(exAtlas),
                         function(i,m){m[,i]},exAtlas),col=sample(pal,ncol(exAtlas),replace = TRUE),
                  xlab="vst expression",main="exAtlas")

plot.multidensity(lapply(1:ncol(norWood),
                         function(i,m){m[,i]},norWood),col=sample(pal,ncol(norWood),replace = TRUE),
                  xlab="vst expression",main="norWood")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
#' 
