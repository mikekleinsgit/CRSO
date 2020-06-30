

### Parallel library: these functions help parallelization by splitting
### matrices and vectors into evenly sized lists.

### 1. splitMatIntoList(mat,n.splits)
### 2. splitVectorIntoList(vec,n.splits)




########################################################################
### Parallel friendly functions
########################################################################

# 1)
splitMatIntoList <- function(mat,n.splits){
  ### Takes a matrix and splits into n.splits evenly sized matrices
  n.size <- nrow(mat)
  n.splits <- min(n.splits,n.size)

  base.size <- floor(n.size/n.splits)
  remainder <- n.size%%n.splits

  list.of.mats <- vector("list",length=n.splits)
  for(j in 1:n.splits){
    if(j != n.splits) list.of.mats[[j]] <- mat[c((base.size*(j-1) + 1):(base.size*j)),]
    if(j == n.splits) list.of.mats[[j]] <- mat[c((base.size*(j-1) + 1):n.size),]
  }
  names(list.of.mats) <- paste0("mat",1:n.splits)
  return(list.of.mats)
}

# 2)
splitVectorIntoList <- function(vec,n.splits){
  ### Takes a matrix and splits into n.splits evenly sized matrices
  n.size <- length(vec)
  n.splits <- min(n.splits,n.size)

  base.size <- floor(n.size/n.splits)
  remainder <- n.size%%n.splits

  list.of.vecs <- vector("list",length=n.splits)
  for(j in 1:n.splits){
    if(j != n.splits) list.of.vecs[[j]] <- vec[c((base.size*(j-1) + 1):(base.size*j))]
    if(j == n.splits) list.of.vecs[[j]] <- vec[c((base.size*(j-1) + 1):n.size)]
  }
  names(list.of.vecs) <- paste0("vec",1:n.splits)
  return(list.of.vecs)
}




