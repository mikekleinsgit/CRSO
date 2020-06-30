

### 1. evaluateIM(D,rm,rm.cm,im,Q,msa,n.splits)
### 2. getTopIm (D,rm,rm.cm,full.im,Q,max.stored,msa,should.check)
### 3. makePhaseTwoImList(D,Q,rm.ordered,k.max,pool.sizes,max.stored,msa,shouldPrint)
### 4. evaluateListOfIMs(D,Q,rm,im.list,msa)







### 1. Get performance of rule sets in IM
### IM is normally a matrix, but this can handle if it is a vector

#' @param D D
#'
#' @param rm rm
#' @param rm.cm rm.cm
#' @param im index matrix
#' @param Q Q 
#' @param msa msa
#' @param n.splits number of splits for parallelization
#'
#' @importFrom foreach foreach getDoParWorkers %dopar%
evaluateIM <- function(D,rm,rm.cm,im,Q,msa,n.splits){
  n.cores <- foreach::getDoParWorkers()
  if(missing(n.splits)) n.splits <- n.cores
  W.0 <- sum(Q)

  if(ncol(as.matrix(im))==1){
    rm.temp <- rm[im,]
    perfs <- getSingleRuleWs(D,rm.temp,Q,msa)
    return(perfs)
  }

  perfs <- rep(NA,length=nrow(im))
  if(nrow(im)<2*n.cores){
    for(j in 1:length(perfs)) perfs[j] <- getWofRS(D,rm,rm.cm,rs.idc=im[j,],Q,msa,W.0)
  }
  if(nrow(im)>=2*n.cores){
    im.list <- splitMatIntoList(im,n.cores)
    im.cur <- 0
    perfs <- foreach::foreach(im.cur = im.list,.combine = 'c',.export=ls(envir=globalenv())) %dopar% {
      mini.perfs <- rep(NA,length=nrow(im.cur))
      for(j in 1:length(mini.perfs)) mini.perfs[j] <- getWofRS(D,rm,rm.cm,rs.idc=im.cur[j,],Q,msa,W.0)
      mini.perfs
    }
  }
  return(perfs)
}

### 2. evaluate the given im and store the top "max.stored" rule sets.
### returns the top im ordered by performance. Also has option to pre-screen for invalid rules
#' Title
#'
#' @param D D
#' @param rm rm
#' @param rm.cm rm.cm
#' @param full.im starting index matrix
#' @param Q Q 
#' @param max.stored max top rules stored
#' @param msa msa
#' @param should.check Check if valid rule sets? Default is false.
#'
#' @return
#' @export
#'
#' @examples
getTopIm <- function(D,rm,rm.cm,full.im,Q,max.stored,msa,should.check){
  if(missing(should.check)) should.check <- FALSE ### should check means do we need to screen im to remove family rules
  if(should.check){
    fm <- getFamMat(rm)
    temp <- apply(full.im,1,function(x)isRSvalid(rm,x,fm))
    full.im <- full.im[which(temp==TRUE),]
  }
  full.perfs <- evaluateIM(D,rm,rm.cm,full.im,Q,msa)
  # idc <- which(full.perfs>0)
  # full.im <- full.im[idc,]
  # full.perfs <- full.perfs[idc]
  # if(nrow(full.im)==0) return(rep(1,ncol(full.im)))
  full.im <- full.im[order(full.perfs,decreasing = TRUE),]
  max.row.idc <- min(max.stored,nrow(full.im))
  top.im <- full.im[1:max.row.idc,]
  return(top.im)
}



### 3. This important function outputs a list of top ims for each k in 1 to k.max subject to pool sizes
#' @title Output list of top rule sets for each k in 1:k.max
#'
#' @param D binary matrix of events by samples
#' @param Q penalty matrix of events by samples
#' @param rm.ordered matrix of rules ordered by phase one
#' @param pool.sizes vector of the number of top rules evaluated for each k
#' @param max.stored max number of rule sets saved
#' @param msa minimum samples assigned
#' @param shouldPrint Print progress updates? Default is TRUE
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.start <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.start,
#'          pool.sizes=c(60,20,20),max.stored=10,msa=15,shouldPrint = TRUE)
#' @export
#' @return largest n such that n choose k < max.num.rs
makePhaseTwoImList <- function(D,Q,rm.ordered,pool.sizes,max.stored,msa,shouldPrint){
  if(missing(shouldPrint)) shouldPrint <- TRUE
  grand.beg <- Sys.time()
  list.of.ims <- vector("list",length=length(pool.sizes))
  names(list.of.ims) <- paste0("K.",1:length(list.of.ims))
  list.of.perfs <- list.of.ims
  ### Do k = 1: assume fewer than max.num.stored
  im <- c(1:nrow(rm.ordered))
  perfs <- getSingleRuleWs(D,rm.ordered,Q,msa)
  list.of.ims[[1]] <- im[order(perfs,decreasing = TRUE)]
  rm.cm.ordered <- makeRSCoverageMat(D,rm.ordered)
  
  ### k = 2 through k max
  if(shouldPrint) print("Starting Phase 2: Exhaustive Evaluation")
  for(k in 2:length(pool.sizes)){
    pool.size <- pool.sizes[k]
    rm <- rm.ordered[1:pool.size,]
    full.im <- makeOneExhaustiveIM(rm,k) ### full im contains only valid rule sets
    if(nrow(full.im)==0){
      print("Reached max K with valid rule set, consider raising pool sizes.")
      return(list.of.ims[1:(k-1)])
    }
    rm.cm <- rm.cm.ordered[1:pool.size,]

    beg <- Sys.time()
    top.im <- getTopIm(D,rm,rm.cm,full.im,Q,max.stored, msa = msa)
    if(shouldPrint) print(paste0("Evaluation time for k = ",k,": ",signif(difftime(Sys.time(),beg,units="min"),4)," min"))
    list.of.ims[[k]] <- top.im
  }
  if(shouldPrint) print(paste0("Total Phase 2 Time: " ,signif(difftime(Sys.time(),grand.beg,units="min"),4)," min"))

  return(list.of.ims)
}


### 4. Evaluate perfs of im list
#' @title Evaluate list of rule set matrices
#'
#' @param D binary matrix of events by samples
#' @param Q penalty matrix of events by samples
#' @param rm matrix of rules ordered by phase one
#' @param im.list list of rule set matrices
#' @param msa minimum samples assigned
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.start <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.start,
#'          pool.sizes=c(60,20,20),max.stored=10,msa=15,shouldPrint = TRUE)
# p2.performance.list <- evaluateListOfIMs(D,Q,rm.full,p2.im.list,msa)
#' @export
#' @return list of Js for each rule set matrix
evaluateListOfIMs <- function(D,Q,rm,im.list,msa){
  rm.cm <- makeRSCoverageMat(D,rm)
  perfs.list <- vector("list",length=length(im.list))
  names(perfs.list) <- names(im.list)
  perfs.list[[1]] <- evaluateIM(D,rm,rm.cm,im.list[[1]],Q,msa)
  for(k in 2:length(im.list)){
    im <- im.list[[k]]
    ### Case where im is 1 row only:
    if(ncol(as.matrix(im))==1) {perfs.list[[k]] <- getWofRS(D,rm,rm.cm,im,Q,msa)
    }else{perfs.list[[k]] <- evaluateIM(D,rm,rm.cm,im,Q,msa)}
  }
  return(perfs.list)
}

#########################################################################################
