#!/usr/bin/env Rscript

library(tidyverse)


calcXS <- function(df) {
    df %>% group_by(run, event)  %>% filter(row_number()==1) %>% select(weight) -> wgts
    sum(wgts$weight) / max(wgts$run)
}

processFile <- function(fname) {
    read.table(fname, 
               col.names=c("run", "event", "id", "charge", "weight", "posx", "posy", "posz",
                           "E", "px", "py", "pz", "history", "prodid", "Enu")
               ) %>% filter(weight>0) %>% calcXS
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

processFile(args[1])
