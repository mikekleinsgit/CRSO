
### 10/15/2018
###################################################################################
### Exact version is now the default version, see old libraries for other versions
###################################################################################

### Not Exported:
### 1. getRCs(D,rm,rm.cm,ss,spr,rule.ids,Q.mat,W.0)
### 2. getRCs.Parallel(D,rm,rm.cm,ss,spr,rule.ids,Q.mat,W.0,n.splits)
### 3. makeScaledRSM(D,rm,spr,Q.mat,n.splits)

### Exported:
### 4. runPhaseOne(D,Q,rm.start,cut.size,spr,trn,n.splits,shouldPrint)



# 1.
### Get rule scores for ss
### rule ids is the rules we want to get scores for
getRCs <- function(D,rm,rm.cm,ss,spr,rule.ids,Q,W.0){
  if(missing(W.0)) W.0 <- sum(Q)
  p1.msa <- 0 ### For phase one we need msa to be 0 
  mat <- matrix(NA,nrow=length(rule.ids),ncol=spr)
  for(r in 1:nrow(mat)){
    rule.id <- rule.ids[r]
    for(i in 1:ncol(mat)){
      rs.idc <- c(rule.id,sample(setdiff(1:nrow(rm),rule.id),ss-1))
      full.perf <- getWofRS(D,rm,rm.cm,rs.idc,Q,p1.msa,W.0)
      perf.without <- getWofRS(D,rm,rm.cm,setdiff(rs.idc,rule.id),Q,p1.msa,W.0)
      score <- 100*(full.perf-perf.without)/full.perf
      mat[r,i] <- score
    }
  }
  scores <- rowMeans(mat)
  names(scores) <- rownames(rm[rule.ids,])
  return(scores)
}


# 2.

#' @importFrom foreach foreach %dopar% getDoParWorkers
getRCs.Parallel <- function(D,rm,rm.cm,ss,spr,rule.ids,Q,W.0,n.splits){
  vec <- NULL
  if(missing(n.splits)) n.splits <- getDoParWorkers()
  if(missing(W.0)) W.0 <- sum(Q)
  #n.splits <- min(n.cores,length(rule.ids))
  vec.list <- splitVectorIntoList(rule.ids,n.splits)
  scores <- foreach::foreach(vec=vec.list,.combine = c,.export=ls(envir=globalenv())) %dopar% {
    getRCs(D,rm,rm.cm,ss,spr,vec,Q,W.0)
  }
  return(scores)
}


# 3)
makeScaledRSM <- function(D,rm,spr,Q,n.splits){
  n.cores <- foreach::getDoParWorkers()
  if(missing(n.splits)) n.splits <- n.cores
  W.0 <- sum(Q)

  ### ss.vec is set here now
  nr <- nrow(rm)
  if(nr > 150)   ss.vec <- c(8,12,16)#ceiling(.01*c(4,7,10)*nr)
  if(nr > 300)   ss.vec <- c(15,25,35)
  if(nr <= 150) ss.vec <- c(6,8,10,12)

  #print(ss.vec)
  n.splits <- min(n.splits,floor((1/2)*nr)) ### The two is arbitrary
  #print(paste0("n.splits = ", n.splits))
  rm.cm <- makeRSCoverageMat(D,rm)
  rsm <- matrix(NA,nrow=nrow(rm),ncol=length(ss.vec))
  rownames(rsm) <- rownames(rm); colnames(rsm) <- paste0("ss.",ss.vec)
  for(j in 1:length(ss.vec)){
    scores <- getRCs.Parallel(D,rm,rm.cm,ss.vec[j],spr,c(1:nrow(rm)),Q,W.0,n.splits)
    rsm[,j] <- scale(scores)
  }
  return(rsm)
}

#' @title Order rules according to phase one importance ranking
#'
#' @param D Binary matrix of N events and M samples
#' @param Q Penalty matrix, negative log of passenger probability matrix.
#' @param rm.start Starting binary rule matrix (i.e., rule library)
#' @param cut.size Fraction of rules eliminated per phase 1 iteration
#' @param spr Random rule sets per rule in each phase one iteration. Default is 40.
#' @param trn Target rule number for stopping iterating. Default is 16.
#' @param n.splits number of splits for parallelization. Default is all available cpus.
#' @param shouldPrint Print progress updates? Default is TRUE
#' @examples
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.06) # Rule library matrix, dimension: 36s x 71
#' rm.ordered <- runPhaseOne(D,Q,rm.full,cut.size=1,spr = 1,trn = 30,shouldPrint = TRUE)
#' # note, for real applications, spr should be at least 10.
#' @export
#' @return binary rule matrix ordered by phase one importance ranking
runPhaseOne <- function(D,Q,rm.start,cut.size,spr,trn,n.splits,shouldPrint){
  if(missing(shouldPrint)) shouldPrint <- TRUE
  if(missing(spr)) spr <- 10
  if(missing(trn)) trn <- 20 ### Target rule number
  if(missing(cut.size)) cut.size <- 0.5
  if(missing(n.splits)) n.splits <- foreach::getDoParWorkers()
  if(nrow(rm.start)<=trn){
    if(shouldPrint) print("RM fewer than TRN. Proceed to phase two with rm.start.")
    return(NULL)
  }
  beg <- Sys.time()
  rm <- rm.start
  order.vec <- c()
  go <- 1
  if(shouldPrint) print(paste0("Starting rules = ",nrow(rm)))

  while(go){
    rsm <- makeScaledRSM(D,rm,spr,Q,n.splits)
    scores <- rowMeans(rsm)

    ntbr <- min(floor(cut.size*nrow(rm)),nrow(rm)-trn) # number of rules to be removed
    bad.thresh <- sort(scores,decreasing=FALSE)[ntbr]
    order.vec <- c(order.vec,names(sort(scores,decreasing=FALSE)[1:ntbr]))
    rm <- rm[which(scores>bad.thresh),]
    if(shouldPrint) print(paste0("Remaining rules = ",nrow(rm)))
    if(nrow(rm)<=trn) go <- 0 ### Should always be equal not less
  }

  rsm <- makeScaledRSM(D,rm,spr,Q,n.splits)
  scores <- rowMeans(rsm)
  order.vec <- c(order.vec,names(sort(scores,decreasing = FALSE)))
  ### Reverse order of order vec
  order.vec <- order.vec[length(order.vec):1]
  rm.final <- rm.start[order.vec,]

  if(shouldPrint) print(difftime(Sys.time(),beg,units="min"))
  return(rm.final)
}


