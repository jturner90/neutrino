#!/usr/bin/env Rscript

library(tidyverse)



translate <- function(id, q, na.rm=TRUE) {
    if      (id==903)          {r <-   15}
    else if (id==901)          {r <-   11}
    else if (id==902)          {r <-   13}
    else if (id==1   & q==0)   {r <-   2112}
    else if (id==1   & q==1)   {r <-   2212}
    else if (id==-1  & q==-1)  {r <-  -2212}
    else if (id==-1  & q==0)   {r <-  -2112}

    else if (id==101 & q==-1)  {r <-  -211}
    else if (id==101 & q==0)   {r <-   111}
    else if (id==101 & q==1)   {r <-   211}
    
    else if (id==114 & q==1)   {r <-   411}
    else if (id==114 & q==0)   {r <-   421}
    
    else if (id==115 & q==-1)  {r <-  -411}
    else if (id==115 & q==0)   {r <-  -421}
    
    else if (id==118 & q==1)   {r <-   431}
    
    else if (id==119 & q==-1)   {r <-  -431}

    else if (id==2   & q==2)   {r <- 2224}
    else if (id==2   & q==1)   {r <- 2214}
    else if (id==2   & q==0)   {r <- 2114}
    else if (id==2   & q==-1)  {r <- 1114}


    else if (id==110  & q==0)   {r <- 311}
    else if (id==110  & q==1)   {r <- 321}
   
    else if (id==32   & q==0)   {r <- 3122}

    else if (id==111  & q==-1)  {r <- -321}
    else if (id==111  & q==0)   {r <- -311}
    
    else if (id==53  & q==0)    {r <-  3322}
    else if (id==53  & q==-1)   {r <-  3312}
    
    else if (id==56  & q==1)   {r <- 4122}

    else if (id==57  & q==2)   {r <-  4222}
    else if (id==57  & q==1)   {r <-  4212}
    else if (id==57  & q==0)   {r <-  4112}
    
    else if (id==59  & q==0)   {r <-  4132}
    else if (id==59  & q==1)   {r <-  4232}
    #else if (id==59  & q==-1)  {r <- -4232}
    
    else if (id==911)            {r<-   12}
    else if (id==-911)           {r<-  -12}
    else if (id==912)            {r<-   14}
    else if (id==-912)           {r<-  -14}
    else if (id==913)            {r<-   16}
    else if (id==-913)           {r<-  -16}
    
    else if (id==33   & q==1)   {r <-  3222}
    else if (id==33   & q==0)   {r <-  3212}
    else if (id==33   & q==-1)  {r <-  3112}
    else if (id==-33  & q==-1)  {r <- -3222}
    else if (id==-33  & q==0)   {r <- -3212}
    else if (id==-33  & q==1)   {r <- -3112}
    
    else if (id==-32   & q==0)  {r <- -3122}
    
    else if (id==-53  & q==0)   {r <- -3322}
    else if (id==-53  & q==1)   {r <- -3312}
    
    else if (id==-56  & q==-1)   {r <- -4122}
    else if (id==55  & q==-1)   {r <- 3334}
    else {
        print(id)
        print(q)
        stop("oh no")
    }
    return(r)
}

radius <- function(x,y,z,na.rm=TRUE) {
    sqrt(x*x + y*y + z*z)
}

calcXS <- function(df) {
    df %>% group_by(run, event)  %>% filter(row_number()==1) %>% select(weight) -> wgts
    sum(wgts$weight) / max(wgts$run)
}

mkDataFrame <- function(fname) {
    read.table(fname, col.names=c("run", "event", "id", "charge", "weight", "posx", "posy", "posz", "E", "px", "py", "pz", "history", "prodid", "Enu")) -> df
    df %>% filter(weight>0) -> df
    xs <- calcXS(df)
    df %>% mutate(xsec=xs) %>% unite(run_event, c(run, event)) -> df
    df %>% mutate(radius=radius(posx, posy, posz)) -> df
    df[, "pdgid"] <- mapply(translate, df[,"id"], df[,"charge"])
    df
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

mkDataFrame(args[1]) -> df
fout_taus      <- paste(args[2],"_tauleptons", ".csv", sep="")
fout_particles <- paste(args[2],"_particles",  ".csv", sep="")
df %>% filter (pdgid==15) %>% select(run_event, Enu, E,px,py,pz) %>% write.table(fout_taus,      row.names=FALSE, col.names=TRUE, sep=",")
df %>% filter (pdgid!=15) %>% select(run_event, pdgid,E,px,py,pz,Enu,weight,xsec) %>% write.table(fout_particles, row.names=FALSE, col.names=TRUE, sep=",")
