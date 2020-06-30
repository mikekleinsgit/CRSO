### This file contains the functions for evaluating the performance of a rule set.
### W refers to the objective function score.

### None of the functions in this file are exported.

### 1) makeOptimalAssignment(rm,rm.cm,rs.idc,Q)
### 2) makePassengerD(D,rm,rs.idc,assign.vec)
### 3) getWofAssignment(D,rm,rs.idc,assign.vec,Q,W.0)
### 4) getWofRS(D,rm,rm.cm,rs.idc,Q,msa,W.0)
### 5) getSingleRuleWs(D,rm,Q,msa)




# 1
#' Title
#'
#' @param rm rm
#' @param rm.cm rm.cm
#' @param rs.idc rs.idc
#' @param Q penalty matrix 
#'
#' @return
#' @export
#'
#' @examples
makeOptimalAssignment <- function(rm,rm.cm,rs.idc,Q){
  if(length(rs.idc)==0) return(rep(0,ncol(rm.cm)))
  if(length(rs.idc)==1) {
    assign.vec <- as.numeric(rm.cm[rs.idc,])
    assign.vec[which(assign.vec==1)] <- rs.idc
    return(assign.vec)
  }
  assign.vec <- rep(0,length=ncol(rm.cm))
  rs.cm <- rm.cm[rs.idc,]
  for(j in 1:ncol(rs.cm)){
    candidate.rule.idc <- rs.idc[which(rs.cm[,j]==1)]
    if(length(candidate.rule.idc)==1) assign.vec[j] <-  candidate.rule.idc
    if(length(candidate.rule.idc)>1){
      ### here we score each rule
      rule.scores <- rep(NA,length=length(candidate.rule.idc))
      for(y in 1:length(rule.scores)){
        rule.id <- candidate.rule.idc[y]
        rule.events <- which(rm[rule.id,]==1)
        rule.scores[y] <- sum(Q[rule.events,j])
      }
      assign.vec[j] <- candidate.rule.idc[which.min(rule.scores)]
    }
  }
  return(assign.vec)
}

# 2.
#' Title
#'
#' @param D D
#' @param rm rm
#' @param rs.idc rs.idc
#' @param assign.vec assignment vector
#'
#' @return
#' @export
#'
#' @examples
makePassengerD <- function(D,rm,rs.idc,assign.vec){
  #rs.mat <- rm[rs.idc,]
  passenger.D <- D
  ### if rs in no rules
  if(length(rs.idc)==0)return(D)
  ### if rs is only one rule
  if(length(rs.idc)==1){
    rule.events.idc <- which(rm[rs.idc,]==1)
    passenger.D[rule.events.idc,which(assign.vec==rs.idc)] <- 0
    return(passenger.D)
  }
  for(j in 1:ncol(D)){
    if(assign.vec[j]>0){
      rule <- assign.vec[j]
      rule.events.idc <- which(rm[rule,]==1)
      passenger.D[rule.events.idc,j] <- 0
    }
  }
  return(passenger.D)
}


# 3.
#' Title
#'
#' @param D binary input D
#' @param rm rule library matrix
#' @param rs.idc rule set indices
#' @param assign.vec assignment vector
#' @param Q Q 
#' @param W.0 null penalty
#'
#' @return
#' @export
#'
#' @examples
getWofAssignment <- function(D,rm,rs.idc,assign.vec,Q,W.0){
  if(missing(W.0)) W.0 <- sum(Q)
  ### We pass along Q.0 to avoid
  D.pass <- makePassengerD(D,rm,rs.idc,assign.vec)
  Q.assigned <- Q
  Q.assigned[which(D.pass==0)] <- 0
  W <- sum(Q.assigned) - W.0
  return(W)
}


# 4.
#' Title
#'
#' @param D D
#' @param rm rule library matrix
#' @param rm.cm rule coverage matrix
#' @param rs.idc rule set indices
#' @param Q Q
#' @param msa msa
#' @param W.0 null penalty
#'
#' @return
#' @export
#'
#' @examples
getWofRS <- function(D,rm,rm.cm,rs.idc,Q,msa,W.0){
  if(missing(W.0)) W.0 <- sum(Q)
  #if(missing(msa)) msa <- 3
  rs.idc <- unique(rs.idc)
  assign.vec <- makeOptimalAssignment(rm,rm.cm,rs.idc,Q)
  temp <- assign.vec[which(assign.vec!=0)]
  if(msa>0){
    if(length(unique(temp))<length(rs.idc)) return(0) ### Rules have 0 assignments
    if(min(table(temp))<msa) return(0)
  }

  #assign.vec <- makeHierarchyAssignment(rm,rm.cm,rs.idc)
  W <- getWofAssignment(D,rm,rs.idc,assign.vec,Q,W.0)
  return(W)
}

# 5.
#' Title
#'
#' @param D D
#' @param rm rm
#' @param Q Q
#' @param msa msa 
#'
#' @return
#' @export
#'
#' @examples
getSingleRuleWs <- function(D,rm,Q,msa){
  #if(missing(msa)) msa <- 3
  W.0 <- sum(Q)
  perfs.vec <- rep(NA,nrow(rm))
  names(perfs.vec) <- rownames(rm)
  rm.cm <- makeRSCoverageMat(D,rm)
  for(rule.id in 1:nrow(rm))perfs.vec[rule.id] <- getWofRS(D,rm,rm.cm,rs.idc=rule.id,Q,msa,W.0)
  return(perfs.vec)
}



