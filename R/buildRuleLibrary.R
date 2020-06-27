
### This file contains the functions needed for building the rule library from D.

### 1. number2binary(number, noBits) - not exported
### 2. getSamplesCovered(D,rule) - not exported
### 3. makeExhaustiveRM(D,n) - not exported
### 4. getNextRM(rm,n) - not exported
### 5. makeRSCoverageMat(D,rs.mat) - not exported
### 6. orderAndTrim(D, freq.thresh) - not exported

### Exported functions:
### 7. buildRuleLibrary(D,rule.thresh,min.epr)
### 8. getRulesAsStrings(rm)




# 1) represent number as binary vector of fixed length
number2binary = function(number, num.bits) {
  if(num.bits > 32){
    stop('num.bits cannot be greater than 32')
  }
  if(number == 0) return(rep(0,num.bits))
  binary_vector = rev(as.numeric(intToBits(number))) ### binary vector of length 32
  min.bits.required <- max(which(binary_vector[32:1]==1))
  if(num.bits < min.bits.required){
    stop('insufficient number of bits')
  }
  return(binary_vector[-(1:(length(binary_vector) - num.bits))])
}

# 2) Returns a binary vector indicating if each sample is covered by the rule.
getSamplesCovered <- function(D,rule){
  ### If it is null rule then all samples are covered
  if(sum(rule)==0) return(rep(1,ncol(D)))
  ### If not null rule, only one gene in rule:
  if(length(which(rule==1))==1) return(D[which(rule==1),])
  covered <- apply(D[which(rule==1),],2,prod)
  return(covered)
}


# 3) Make rule mat for D that consists of all rules
#     containing first n events of D (2^n - 1)
#     output is rm, which is a binary rule matrix
makeExhaustiveRM <- function(D,n){
  N <- nrow(D)
  rm <- matrix(NA,nrow=2^n,ncol=n)
  for(j in 1:nrow(rm))   rm[j,] <- number2binary(j-1,n)
  rm <- cbind(rm,matrix(0,nrow=nrow(rm),ncol=N-n))
  colnames(rm) <- rownames(D)
  rm <- rm[2:nrow(rm),] ### the zero rule is excluded
  rownames(rm) <- paste0("r",1:nrow(rm))
  return(rm)
}

# 4) add n+1 event to rm
getNextRM <- function(rm,n){
  n.plus.1 <- c(rep(0,nrow(rm)),rep(1,nrow(rm)))
  next.rm <- rbind(rm,rm)
  next.rm[,n+1] <- n.plus.1
  ### add rule that is only new event
  rule <- rep(0,ncol(rm)); rule[n+1] <- 1
  rm <- rbind(next.rm,rule)
  return(rm)
}


# 5) This function returns a k x N coverage matrix indicating which rules each sample
# is covered by. A rule set is represented as a k x N matrix, so that the full rm is not necessary.
makeRSCoverageMat <- function(D,rs.mat){
  ### The following if statement checks if the rule set is of size 1
  if(ncol(as.matrix(rs.mat))==1) return(getSamplesCovered(D,rs.mat))
  cov.mat <- t(apply(rs.mat,1,function(x)getSamplesCovered(D,x)))
  return(cov.mat)
}

# 6) order D by frequency and exclude events below frequency threshold
orderAndTrim <- function(D,freq.thresh){
  if(missing(freq.thresh)) freq.thresh <- 0
  freqs <- apply(D,1,sum)/ncol(D)
  idc <- which(freqs>=freq.thresh)
  if(length(idc)<=1) stop('Error 1 or 0 events above thresh')
  D <- D[which(freqs>=freq.thresh),]
  freqs <- apply(D,1,sum)/ncol(D)
  D <- D[order(freqs,decreasing=TRUE),]
  return(D)
}

# 7)
#' @title Make full rule library of all rules that satisfy minimum coverage threshold.
#'
#' @param D Binary matrix of N events and M samples
#' @param rule.thresh Minimum fraction of rules covered. Default is .03
#' @param min.epr minimum events per rule. Default is 2.
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # build rule library
#' dim(rm.full) # Should be matrix with dimension 60 x 71
#' @export
buildRuleLibrary <- function(D,rule.thresh,min.epr){
  ### rule.thresh is fraction of samples covered, e.g., .05
  if(missing(min.epr)) min.epr <- 2
  if(missing(rule.thresh)) rule.thresh <- .03
  ### First thing we need to make sure every event satisfies rule.thresh
  ### If an event does not satisfy rule thresh it cannot be in any rule that does
  original.D <- D

  D <- orderAndTrim(original.D,rule.thresh)
  if(identical(D,-1)) stop('no eligible rules')
  if(nrow(D)<=2) stop('no eligible rules')
  D <- D[nrow(D):1,]  ### For computational efficiency we consider events according to decreasing frequency
  max.n <- nrow(D)

  ### Make rm of first 2 events
  rm <- makeExhaustiveRM(D,2)
  all.cov.mat <- makeRSCoverageMat(D,rm)
  valid.idc <- which(rowMeans(all.cov.mat)>=rule.thresh) ### we know at least 2 rules satisfy this (the ones that contain 1 event each)
  n.rules <- length(valid.idc)
  rm <- rm[valid.idc,]
  prev.cov.mat <- all.cov.mat[valid.idc,]
  ### n.rules and prev.cov.mat are used to make this faster

  for(n in 3:(max.n)){
    rm <- getNextRM(rm,n-1)
    new.rm <- rm[c((n.rules+1):nrow(rm)),]
    new.cov.mat <- makeRSCoverageMat(D,new.rm)
    all.cov.mat <- rbind(prev.cov.mat,new.cov.mat)
    #all.cov.mat <- makeRSCoverageMat(D,rm)
    valid.idc <- which(rowMeans(all.cov.mat)>=rule.thresh)
    n.rules <- length(valid.idc)
    rm <- rm[valid.idc,]
    prev.cov.mat <- all.cov.mat[valid.idc,]
  }

  n.events <- nrow(original.D)
  rm.add.on <- matrix(0,nrow=nrow(rm),ncol=n.events-ncol(rm))
  add.on.names <- setdiff(rownames(original.D),rownames(D))
  rm.full <- cbind(rm,rm.add.on)
  colnames(rm.full) <- c(colnames(rm),add.on.names)

  ### Need the ordering of rm.full columns to match the order of original.D rows
  rm.full <- rm.full[,rownames(original.D)]
  if (!identical(rownames(original.D),colnames(rm.full))) stop('Names mismatch')
  ### Need rules with between min.epr and max.epr events
  events.per.rule <- rowSums(rm.full)
  idc <- which(events.per.rule >= min.epr)
  if(length(idc)<=1)stop('no eligible rules')
  rm.full <- rm.full[idc,]
  ### Order rm according to coverage

  rm.full <- rm.full[order(rowSums(makeRSCoverageMat(original.D,rm.full)),decreasing = T),]
  rownames(rm.full) <- paste0("r",1:nrow(rm.full))
  return(rm.full)
}

# 8)
#' @title Represent binary rule matrix as strings
#'
#' @param rm binary rule matrix
#' @return vector or rules represented as strings
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.1) # Small rule library matrix, dimension: 5 x 71
#' getRulesAsStrings(rm.full)
#' # output should be: "BRAF-M.CDKN2A-MD"   "CDKN2A-MD.NRAS-M"
#' # "BRAF-M.PTEN-MD"    "ADAM18-M.BRAF-M" "ADAM18-M.CDKN2A-MD"
#' @export
getRulesAsStrings <- function(rm){
  if(nrow(t(rm))==1){
    return(paste0(sort(names(which(rm==1))),collapse = "."))
  }
  rm.strings <- rep(NA,nrow(rm))
  names(rm.strings) <- rownames(rm)
  for(j in 1:length(rm.strings)){
    rm.strings[j] <- paste0(sort(names(which(rm[j,]==1))),collapse = ".")
  }
  return(rm.strings)
}



