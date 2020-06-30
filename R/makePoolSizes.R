
### This file contains the functions for determining the phase two pool sizes
### There are two ways to determine pool sizes for a given k and max.nrs:
### Fast way: does not consider family memberships
### Slow way: considers family memberships,
### and increases the pool sizes accordingly.

### 1. getMaxNforK(k,max.num.rs) - not exported
### 2. extractValidIM.OC(rm,fm,im)
### 3. extractValidIM(rm,fm,im,n.splits)
### 4. getOneExhaustiveIM(rm,k)
### 5. getPSforK(rm.ordered,k,max.nrs.ee,max.compute)
### 6. getPoolSizes(rm.ordered,k.max,max.nrs.ee,max.compute){


#getPoolSizeForK(rm.ordered,max.nrs,k,fastWay)
#(rm.ordered,k.max,max.num.rs.ee)


### 1.
### This function takes as input a k value and an upper bound, max.num.rs,
### and it findes the largest N such that (N choose K) < max.num.rs
getMaxNforK <- function(k,max.num.rs){
  if(missing(max.num.rs)) max.num.rs <- 20000
  go <- 1
  N <- k+1
  while(go){
    temp <- choose(N,k)
    if(temp > max.num.rs) go <- 0
    N <- N + 1
  }
  max.N <- N-2 ### N-2 was last one under max.num.rs
  return(max.N)
}

### 2.
extractValidIM.OC <- function(rm,fm,im){
  fam.rules <- which(rowSums(fm)>1)
  fam.rules <- intersect(fam.rules,im) ### ignore rules that don't appear in im
  if(length(fam.rules)==0) return(im)
  for(rule in fam.rules){
    rule.fams <- setdiff(which(fm[rule,]==1),rule)
    rule.fams <- setdiff(rule.fams,c(1:rule))
    if(length(rule.fams)>0){
      rule.idc <- which(apply(im,1,function(x)rule %in% x)==TRUE)
      im.rule <- im[rule.idc,]
      if(length(rule.idc)==1) im.rule <- t(as.matrix(im.rule))
      if(length(rule.idc)==0) next
      
      im.nonrule <- im[setdiff(1:nrow(im),rule.idc),]

      bad.idc <- which(apply(im.rule,1,function(x)length(intersect(rule.fams,x)))>0)
      new.im.rule <- im.rule[setdiff(1:nrow(im.rule),bad.idc),]
      new.im <- rbind(im.nonrule,new.im.rule)
      im <- new.im
      #print(dim(im))
    }
  }
  return(im)
}

### 3.

#' @importFrom foreach getDoParWorkers foreach %dopar%
extractValidIM <- function(rm,fm,im,n.splits){
  #im.cur <- NULL
  n.cores <- getDoParWorkers()
  if(missing(n.splits)) n.splits <- min(n.cores,nrow(im))
  im.list <- splitMatIntoList(im,n.splits)
  new.im <- foreach(im.cur = im.list,.combine = 'rbind',.multicombine = TRUE,.export=ls(envir=globalenv())) %dopar% {
    extractValidIM.OC(rm,fm,im.cur)
  }
  return(new.im)
}

### 4.
makeOneExhaustiveIM <- function(rm,k){
  im.raw <- t(combn(nrow(rm),k))
  fm <- getFamMat(rm)
  #im <- extractValidIM(rm,fm,im.raw) -- need to figure out whats wrong
  im <- extractValidIM.OC(rm,fm,im.raw)

  return(im)
}

### 5.
getPSforK <- function(rm.ordered,k,max.nrs.ee,max.compute){
  #if(missing(max.compute)) max.compute <- 10^6
  #if(missing(max.nrs.ee)) max.nrs.ee <- 10^5
  this.max.compute <- min(max.compute,20*max.nrs.ee)
  if(k <= 6) this.max.compute <- min(max.compute,5*max.nrs.ee)

  beg <- Sys.time()
  max.n <- getMaxNforK(k,this.max.compute)
  max.n <- min(max.n,nrow(rm.ordered))
  rm <- rm.ordered[1:max.n,]
  im.raw <- t(combn(nrow(rm),k))
  fm <- getFamMat(rm)
  ### Easier way, but maybe slower for large k
  temp <- apply(im.raw,1,function(x)isRSvalid(rm,x,fm))
  im <- im.raw[which(temp==TRUE),]
  #if(sum(temp)==0) return(0)


  #im <- makeOneExhaustiveIM(rm,k)

  if(nrow(im)==0){
    pool.size <- 0
    nrspk <- 0
  }
  if(nrow(im)>0){
    all.maxs <- apply(im,1,max)
    num.rs.below <- rep(NA,length=max.n)
    for(j in 1:max.n) num.rs.below[j] <- length(which(all.maxs <= j))
    pool.size <- max(which(num.rs.below <= max.nrs.ee))
    im <- im[which(all.maxs<= pool.size),]
    nrspk <- nrow(im)
  }
  #print(difftime(Sys.time(),beg,units = "min"))
  #print(c(pool.size,nrspk))
  return(pool.size)
}

# 6.
#' @title Get pool sizes for phase 2
#'
#' @param rm.ordered binary rule matrix ordered from phase 1
#' @param k.max maximum rule set size
#' @param max.nrs.ee max number of rule sets per k
#' @param max.compute maximum raw rule sets considered per k
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' rm.ordered <- rm.full # Skip phase one in this example
#' getPoolSizes(rm.ordered,k.max = 7,max.nrs.ee = 10000)
#' # [1] 60  60  40  23  18  16  15
#' @export
getPoolSizes <- function(rm.ordered,k.max,max.nrs.ee,max.compute){
  if(missing(k.max)) k.max <- 12
  if(missing(max.nrs.ee)) max.nrs.ee <- 10^5
  if(missing(max.compute)) max.compute <- 10*max.nrs.ee
  #beg <- Sys.time()
  pool.sizes <- rep(NA,length=k.max)
  names(pool.sizes) <- paste0("K.",1:k.max)
  pool.sizes[1] <- nrow(rm.ordered)
  for(k in 2:k.max) pool.sizes[k] <- getPSforK(rm.ordered,k,max.nrs.ee,max.compute)
  #print(paste0("Time to make pool sizes: " ,signif(difftime(Sys.time(),beg,units="min"),4)," min"))
  return(pool.sizes)
}
