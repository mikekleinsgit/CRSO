
# Exported:
### 1. makePhaseThreeImList(D,Q,rm.ordered,til.ee,pool.sizes,max.stored,max.nrs.p3,msa,shouldPrint)

# Not exported
### 2. makeFullD1.IM(im.d1.core,borrow.pool,K)
### 3. makeFullD2.IM(im.d2.core,borrow.pool,K)
### 4. makeFullD3.IM(im.d3.core,borrow.pool,K)

### 5 updateTopIM.K3(D,Q.mat,rm,pool.sizes,max.nrs,til,max.stored,msa)
### 6. updateTopIM(D,Q.mat,rm,fm,pool.sizes,max.nrs,til,max.stored,K,msa)



#' @title Make phase 3 im list from phase 2 im list
#'
#' @param D binary matrix of events by samples
#' @param Q penalty matrix of events by samples
#' @param rm.ordered matrix of rules ordered by phase one
#' @param fm.ordered family matrix of rm.ordered
#' @param til.ee list of rule set matrices (im list) from phase two
#' @param pool.sizes pool sizes for phase two
#' @param max.stored max number of rule sets saved
#' @param max.nrs.p3 max number of new rule sets per k, default is 10^5
#' @param msa minimum samples assigned per rule in rule set
#' @param shouldPrint Print progress updates? Default is TRUE
#'
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.full,pool.sizes=c(60,10,10),
#'           max.stored=10,msa=15,shouldPrint = FALSE)
#' fm.ordered <- getFamMat(rm.full)
#' til.p3 <- makePhaseThreeImList(D,Q,rm.ordered = rm.full,fm.ordered,
#' til.ee = til.p2, pool.sizes=c(60,20,20),
#'          max.stored=100,max.nrs.p3=100,msa=15,shouldPrint = TRUE)
#' @export
#' @return phase 3 top im list
makePhaseThreeImList <- function(D,Q,rm.ordered,fm.ordered,til.ee,pool.sizes,max.stored,max.nrs.p3,msa,shouldPrint){
  if(missing(shouldPrint)) shouldPrint <- TRUE
  if(missing(max.nrs.p3)) max.nrs.p3 <- 10^5
  if(shouldPrint) print("Starting Phase 3: neighbor expansion")
  ### Step 3: make updated top im lists
  grand.beg <- Sys.time()
  til.final <- til.ee ### til = top im list
  for(K in 3:length(til.final)){
    beg <- Sys.time()
    til.final[[K]] <- updateTopIM(D,Q,rm.ordered,fm.ordered,pool.sizes,max.nrs.p3,til.final,max.stored,K,msa)
    if(shouldPrint) print(paste0("Updated top im for K = ",K,", time = ",signif(difftime(Sys.time(),beg,units="min"),4)," min"))
  }
  if(shouldPrint) print(paste0("Total Phase 3 Time: " ,signif(difftime(Sys.time(),grand.beg,units="min"),4)," min"))
  return(til.final)
}

### 2. This function takes in a matrix of unique cores of length K-1, 
### and make IM of adding each borrow rule to every core
### All rs are valid
makeFullD1.IM <- function(rm,fm,im.d1.core,full.borrow.pool,K){
  full.d1.im <- matrix(NA,nrow=0,ncol=K)
  for(j in 1:nrow(im.d1.core)){
    cur.core <- im.d1.core[j,]
    fam.idc <- which(colSums(fm[cur.core,])>0)
    borrow.pool <- setdiff(full.borrow.pool,fam.idc)
    core.im <- matrix(NA,nrow=length(borrow.pool),ncol=K)
    for(z in 1:nrow(core.im)) core.im[z,] <- c(cur.core,borrow.pool[z])
    full.d1.im <- rbind(full.d1.im,core.im)
  }
  full.d1.im <- t(apply(full.d1.im,1,sort))
  return(full.d1.im)
}

### 3. Takes a core matrix of K-2 cores, and a borrow pool and makes all distance two rule sets
### Checks that all rs are valid
#' @importFrom utils combn
makeFullD2.IM <- function(rm,fm,im.d2.core,full.borrow.pool,max.nrs.d2,K){
  full.d2.im <- matrix(NA,nrow=0,ncol=K)
  if(max.nrs.d2 <= 2) return(full.d2.im)
  
  #partial.borrow.im <- t(combn(borrow.pool,2))
  for(j in 1:nrow(im.d2.core)){
    cur.core <- im.d2.core[j,]
    fam.idc <- which(colSums(fm[cur.core,])>0)
    borrow.pool <- setdiff(full.borrow.pool,fam.idc)
    bp.length <- min(getMaxNforK(2,max.nrs.d2),length(borrow.pool))
    borrow.pool <- borrow.pool[1:bp.length]
    core.im <- t(combn(borrow.pool,2))        
    ### Check here that core.im of new rules is valid
    is.valid.vec <- apply(core.im,1,function(x)isRSvalid(rm,x,fm))
    num.valid <- sum(is.valid.vec)
    if(num.valid==0) next
    if(num.valid==1){
      core.im <- core.im[which(is.valid.vec==TRUE),]
      core.im <- c(cur.core,core.im)
      full.d2.im <- rbind(full.d2.im,core.im)
      next
    }
    core.im <- core.im[which(is.valid.vec==TRUE),]
    core.im <- t(apply(core.im,1,function(x)c(cur.core,x)))
    #for(z in 1:nrow(core.im)) core.im[z,] <- c(cur.core,partial.borrow.im[z,])
    full.d2.im <- rbind(full.d2.im,core.im)
  }
  full.d2.im <- t(apply(full.d2.im,1,sort))
  return(full.d2.im)
}

### 4. Takes a core matrix of K-3 cores, and a borrow pool and makes all distance 3 rule sets
### all distance 3 rule sets are valid
#' @importFrom utils combn
makeFullD3.IM <- function(rm,fm,im.d3.core,full.borrow.pool,max.nrs.d3,K){
  full.d3.im <- matrix(NA,nrow=0,ncol=K)
  if(max.nrs.d3 <= 3) return(full.d3.im)
  for(j in 1:nrow(im.d3.core)){
    
    cur.core <- im.d3.core[j,]
    if(K>4)fam.idc <- which(colSums(fm[cur.core,])>0)
    if(K==4)fam.idc <- which(fm[cur.core,]>0)
    
    borrow.pool <- setdiff(full.borrow.pool,fam.idc)
    bp.length <- min(getMaxNforK(3,max.nrs.d3),length(borrow.pool))
    borrow.pool <- borrow.pool[1:bp.length]
    core.im <- t(combn(borrow.pool,3))
    ### Check here that core.im of new rules is valid
    is.valid.vec <- apply(core.im,1,function(x)isRSvalid(rm,x,fm))
    num.valid <- sum(is.valid.vec)
    if(num.valid==0) next
    if(num.valid==1){
      core.im <- core.im[which(is.valid.vec==TRUE),]
      core.im <- c(cur.core,core.im)
      full.d3.im <- rbind(full.d3.im,core.im)
      next
    }
    core.im <- core.im[which(is.valid.vec==TRUE),]
    core.im <- t(apply(core.im,1,function(x)c(cur.core,x)))
    full.d3.im <- rbind(full.d3.im,core.im)
  }
  full.d3.im <- t(apply(full.d3.im,1,sort))
  return(full.d3.im)
}

### 5. Special handling of K = 3

#' @importFrom utils combn
updateTopIM.K3 <- function(D,Q,rm,fm,pool.sizes,max.nrs,til,max.stored,msa){
  K <- 3
  top.im <- til[[K]]
  prev.im <- til[[K-1]]
  if(nrow(rm) <= pool.sizes[K]) return(top.im) ### if all rules have been exhaustively analyzed, should alwasy be = not less
  rm.cm <- makeRSCoverageMat(D,rm)
  first.new.rule <- pool.sizes[K] + 1
  cur.top.rs <- top.im[1,] ### Length K
  prev.top.rs <- prev.im[1,] ### Length K-1

  ### Special handling, only do this for K = 3:
  all.best.rs <- union(cur.top.rs,prev.top.rs)
  full.borrow.pool <- c(first.new.rule:nrow(rm)) ### always can do full for dist1

  ### Make D1 core:
  im.d1.core <- t(combn(cur.top.rs,K-1))
  im.d1.core <- rbind(im.d1.core,prev.top.rs)
  im.d1.core <- im.d1.core[!duplicated(im.d1.core),]
  ### Get dist 1 top im:
  full.d1.im <- makeFullD1.IM(rm,fm,im.d1.core,full.borrow.pool,K)
  full.d1.im <- t(apply(full.d1.im,1,sort))
  top.im.d1 <- getTopIm(D,rm,rm.cm,full.d1.im,Q,max.stored,msa,should.check = FALSE)
  
  ### Get dist 2 top im:
  im.d2.core <- t(combn(cur.top.rs,K-2))   ### Let's get all k-1 cores
  im.d2.core <- rbind(im.d2.core,t(combn(prev.top.rs,K-2)))
  im.d2.core <- im.d2.core[!duplicated(im.d2.core),] ### just a vector of single rules
  max.nrs.d2 <- floor(max.nrs/length(im.d2.core)) ### Adjust max nrs for number of cores
  ### Make full im.d2 in this cas
  full.d2.im <- matrix(NA,nrow=0,ncol=K)
  for(j in 1:length(im.d2.core)){
    cur.core <- im.d2.core[j]
    fam.idc <- which(fm[cur.core,]>0)
    borrow.pool <- setdiff(full.borrow.pool,fam.idc)
    bp.length <- min(getMaxNforK(2,max.nrs.d2),length(borrow.pool))
    borrow.pool <- borrow.pool[1:bp.length]
    core.im <- t(combn(borrow.pool,2))
    core.im <- t(apply(core.im,1,function(x)c(cur.core,x)))
    full.d2.im <- rbind(full.d2.im,core.im)
  }
  full.d2.im <- t(apply(full.d2.im,1,sort))
  ### Check validity here
  is.valid.vec <- apply(full.d2.im,1,function(x)isRSvalid(rm,x,fm))
  full.d2.im <- full.d2.im[which(is.valid.vec==TRUE),]
  top.im.d2 <- getTopIm(D,rm,rm.cm,full.d2.im,Q,max.stored,msa,should.check = FALSE)
  
  if(max.stored < nrow(top.im)) top.im <- top.im[1:max.stored,]
  new.top.im <- rbind(top.im,top.im.d1,top.im.d2)
  new.top.im <- new.top.im[!duplicated(new.top.im),]
  new.top.im <- getTopIm(D,rm,rm.cm,new.top.im,Q,max.stored,msa,should.check = FALSE)
  return(new.top.im)
}


### 6. For all K, special handling of K = 3
updateTopIM <- function(D,Q,rm,fm,pool.sizes,max.nrs,til,max.stored,K,msa){
  if(K==3) return(updateTopIM.K3(D,Q,rm,fm,pool.sizes,max.nrs,til,max.stored,msa))
  top.im <- til[[K]]
  prev.im <- til[[K-1]]
  if(nrow(rm) <= pool.sizes[K]) return(top.im) ### if all rules have been exhaustively analyzed, should alwasy be = not less
  rm.cm <- makeRSCoverageMat(D,rm)
  cur.top.rs <- top.im[1,] ### Length K
  prev.top.rs <- prev.im[1,] ### Length K-1
  full.borrow.pool <- c((pool.sizes[K] + 1):nrow(rm)) ### always can do full for dist1

  ### Make D1 core:
  im.d1.core <- t(combn(cur.top.rs,K-1))
  im.d1.core <- rbind(im.d1.core,prev.top.rs)
  im.d1.core <- im.d1.core[!duplicated(im.d1.core),]
  ### Get dist 1 top im:
  full.d1.im <- makeFullD1.IM(rm,fm,im.d1.core,full.borrow.pool,K)
  top.im.d1 <- getTopIm(D,rm,rm.cm,full.d1.im,Q,max.stored,msa,should.check = FALSE)

  ### Get dist 2 top im:
  im.d2.core <- t(combn(cur.top.rs,K-2))   ### Let's get all k-1 cores
  im.d2.core <- rbind(im.d2.core,t(combn(prev.top.rs,K-2)))
  im.d2.core <- as.matrix(im.d2.core[!duplicated(im.d2.core),])
  max.nrs.d2 <- floor(max.nrs/nrow(im.d2.core)) ### Adjust max nrs for number of cores
  ### Get borrow pool based on max nrs
  full.d2.im <- makeFullD2.IM(rm,fm,im.d2.core,full.borrow.pool,max.nrs.d2,K)
  top.im.d2 <- getTopIm(D,rm,rm.cm,full.d2.im,Q,max.stored,msa,should.check = FALSE)


  ### Get dist 3 top im:
  im.d3.core <- t(combn(cur.top.rs,K-3))   ### Let's get all k-1 cores
  im.d3.core <- rbind(im.d3.core,t(combn(prev.top.rs,K-3)))
  im.d3.core <- im.d3.core[!duplicated(im.d3.core),]
  if(K==4) im.d3.core <- as.matrix(im.d3.core) ### because duplicated makes it vector
  max.nrs.d3 <- max.nrs/nrow(im.d3.core) ### Adjust max nrs for number of cores
  ### Get borrow pool based on max nrs
  full.d3.im <- makeFullD3.IM(rm,fm,im.d3.core,full.borrow.pool,max.nrs.d3,K)
  
  top.im.d3 <- getTopIm(D,rm,rm.cm,full.d3.im,Q,max.stored,msa,should.check = FALSE)

  if(max.stored < nrow(top.im)) top.im <- top.im[1:max.stored,]
  new.top.im <- rbind(top.im,top.im.d1,top.im.d2,top.im.d3)
  new.top.im <- new.top.im[!duplicated(new.top.im),]
  new.top.im <- getTopIm(D,rm,rm.cm,new.top.im,Q,max.stored,msa,should.check = FALSE)
  return(new.top.im)
}


