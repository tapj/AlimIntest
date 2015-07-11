## ----library_loading, include=FALSE--------------------------------------
user=Sys.info()[["user"]]
require(rmarkdown)
require(markdown)
library(ade4)
library(RColorBrewer)
library(BiotypeR)
library(gridExtra)
library(reshape2)
library(vegan)
library(knitr)
library(ggplot2)
library(scales)
library(dplyr)
#source("src/mclapply.hack.R")

## ----knit_options, include=FALSE-----------------------------------------
opts_chunk$set(fig.width=12, fig.height=8, warning=FALSE, 
message=FALSE, echo=FALSE, dev=c('png', 'pdf', 'win.metafile', 'tiff'), 
fig.cap=NULL, cache=TRUE, dpi=300)

