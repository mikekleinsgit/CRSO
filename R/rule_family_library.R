
### This file contains functions for making family matrix (fm) associated with rm.

### None of the functions are exported
### 1. getFamMat(rm)
### 2. isRSvalid(rm,rs,fam.mat)
### 3. getDaughtersOfRule(rule,rm)
### 4. getParentsOfRule(rule,rm)

# 1) This function takes in rule space and makes matrix of family members.
# Rules are considered to be family members with themselves.

#' @title Get family relationship matrix
#'
#' @param rm rule matrix
#'
#' @return
#' @export
#'
#' @examples
getFamMat <- function(rm){
  fam.mat <- matrix(0,nrow=nrow(rm),ncol=nrow(rm))
  rownames(fam.mat) <- rownames(rm)
  colnames(fam.mat) <- rownames(rm)

  for(j in 1:nrow(fam.mat)){
    daughters <- getDaughtersOfRule(rm[j,],rm)
    parents <- getParentsOfRule(rm[j,],rm)
    fams <- setdiff(union(daughters,parents),NA)
    fams <- c(j,fams)
    fam.mat[j,fams] <- 1
    #if(j%%50==0)print(j)
  }
  return(fam.mat)
}

# 2) Check if rule set is valid, i.e., there are no two family members in the rule set
# rs is represented as a vector of indices of rm
isRSvalid <- function(rm,rs,fm){
  ### Assume rs in indices of rules in rm
  if(length(rs)==1) return(TRUE)
  for(j in 1:length(rs)){
    cur.rule <- rs[j]
    fams <- setdiff(which(fm[cur.rule,]==1),cur.rule)
    if(length(intersect(fams,rs))>0) return(FALSE)
  }
  return(TRUE)
}



# 3)
getDaughtersOfRule <- function(rule,rm){
  off.idc <- sort(which(rule==0))
  if(length(off.idc)==0) return(c(1:(nrow(rm)-1))) ### If all are on then everything is daughter
  daughter.idc <- c()

  for(j in 1:nrow(rm)){
    temp <- intersect(which(rm[j,]==0),off.idc)
    if(length(temp)==length(off.idc))   {
      if(!identical(rm[j,],rule)) daughter.idc <- c(daughter.idc,j)
    }
  }
  if(length(daughter.idc)==0)return(NA)
  return(daughter.idc)
}

# 4)
getParentsOfRule <- function(rule,rm){
  on.idc <- sort(which(rule==1))
  if(length(on.idc)==length(rule)) return(NA) ### If all are on then everything is daughter
  parent.idc <- c()
  for(j in 1:nrow(rm)){
    temp <- intersect(which(rm[j,]==1),on.idc)
    if(length(temp)==length(on.idc))   {
      if(!identical(rm[j,],rule)) parent.idc <- c(parent.idc,j)
    }
  }
  if(length(parent.idc)==0) return(NA)
  return(parent.idc)
}
