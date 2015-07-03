## ----load_library_data---------------------------------------------------

library(ggplot2)
library(scales)
library(ade4)
library(reshape2)
library(grid)
library(knitr)
data(metatrans_alimintest)


## ----knit_options, include=FALSE-----------------------------------------
opts_chunk$set(fig.width=12, fig.height=8, warning=FALSE, 
message=FALSE, echo=FALSE, dev=c('png', 'pdf', 'win.metafile', 'tiff'), 
fig.cap=NULL, cache=TRUE, dpi=300)

