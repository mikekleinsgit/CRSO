
### 1. addOneKtoTIL(D,Q,rm,rm.cm,til,max.stored, max.nrs,K,msa)
### 2. makePhaseFourImList(D, Q, rm.ordered, p3.iml,k.max.4,max.stored.p4,max.nrs.p4,msa,shouldPrint)

### Phase 4: add next K to til
addOneKtoTIL <- function(D,Q,rm,rm.cm,til,max.stored, max.nrs,K,msa){
  if(K != length(til)+1) print("Something is wrong")
  prev.im <- til[[K-1]]
  prev.top.rs <- prev.im[1,] ### Length K-1
  #rm.cm <- makeRSCoverageMat(D,rm)
  full.borrow.pool <- c(1:nrow(rm))
  im.d1.core <- prev.im[1:10,]
  fm <- getFamMat(rm)
  
  ### Get dist 1 top im:
  full.d1.im <- makeFullD1.IM(rm,fm,im.d1.core,full.borrow.pool,K)
  full.d1.im <- full.d1.im[!duplicated(full.d1.im),]
  top.im.d1 <- getTopIm(D,rm,rm.cm,full.d1.im,Q,max.stored,msa)
  
  ### Get dist 2 top im:
  cur.top.rs <- top.im.d1[1,]
  im.d2.core <- t(combn(cur.top.rs,K-2))   ### Let's get all k-1 cores
  im.d2.core <- rbind(im.d2.core,t(combn(prev.top.rs,K-2)))
  im.d2.core <- as.matrix(im.d2.core[!duplicated(im.d2.core),])
  max.nrs.d2 <- floor(max.nrs/nrow(im.d2.core)) ### Adjust max nrs for number of cores
  ### Get borrow pool based on max nrs
  full.d2.im <- makeFullD2.IM(rm,fm,im.d2.core,full.borrow.pool,max.nrs.d2,K)
  full.im <- rbind(top.im.d1,full.d2.im)
  full.im <- full.im[!duplicated(full.im),]
  top.im <- getTopIm(D,rm,rm.cm,full.im,Q,max.stored,msa)
  
  ### Get dist 3 top im:
  cur.top.rs <- top.im[1,]
  im.d3.core <- t(combn(cur.top.rs,K-3))   ### Let's get all k-1 cores
  im.d3.core <- rbind(im.d3.core,t(combn(prev.top.rs,K-3)))
  im.d3.core <- im.d3.core[!duplicated(im.d3.core),]
  #if(K==4) im.d3.core <- as.matrix(im.d3.core) ### because duplicated makes it vector
  max.nrs.d3 <- floor(max.nrs/nrow(im.d3.core)) ### Adjust max nrs for number of cores
  ### Get borrow pool based on max nrs
  full.d3.im <- makeFullD3.IM(rm,fm,im.d3.core,full.borrow.pool,max.nrs.d3,K)
  full.im <- rbind(top.im,full.d3.im)
  full.im <- full.im[!duplicated(full.im),]
  final.top.im <- getTopIm(D,rm,rm.cm,full.im,Q,max.stored,msa)
  
  if(max.stored < nrow(final.top.im)) final.top.im <- final.top.im[1:max.stored,]
  til[[K]] <- final.top.im
  names(til)[K] <- paste0("K.",K)
  return(til)
}

#' Title
#'
#' @param D D
#' @param Q Q
#' @param rm.ordered rm.ordered
#' @param p3.iml p3 index matrix list
#' @param k.max.4 k max 4
#' @param max.stored.p4 max rules stored per K
#' @param max.nrs.p4 max rules evaluated
#' @param msa minimum samples assigned
#' @param shouldPrint Print progress updates? Default is TRUE
#'
#' @return
#' @export
#'
#' @examples
makePhaseFourImList <- function(D, Q, rm.ordered, p3.iml,k.max.4,max.stored.p4,max.nrs.p4,msa,shouldPrint){
  if(missing(shouldPrint)) shouldPrint <- TRUE
  til <- p3.iml
  rm.cm.ordered <- makeRSCoverageMat(D,rm.ordered)
  range <- c((length(p3.iml)+1):k.max.4)
  grand.beg <- Sys.time()
  stopflag <- FALSE
  for(K in range){
    beg <- Sys.time()
    til <- addOneKtoTIL(D,Q,rm.ordered,rm.cm.ordered,til,max.stored.p4, max.nrs.p4,K,msa)
    rs.idc <- til[[K]][1,]
    if(shouldPrint) print(paste0("Time for adding K = ",K,": ",difftime(Sys.time(),beg,units = "min")))
    score <- getWofRS(D,rm.ordered,rm.cm.ordered,rs.idc,Q,msa)
    if(score==0){
      til <- til[1:(K-1)]
      print(paste0("no rs satisfies K = ",K))
      stopflag <- TRUE
      break
    }
    if(stopflag)break 
  }
  if(shouldPrint) print(paste0("Total Phase 4 Time = ",difftime(Sys.time(),grand.beg,units = "min")))
  p4.iml <- til
  return(p4.iml)
}


