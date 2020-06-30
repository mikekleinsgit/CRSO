

### Convenience library: these functions can be used in general contexts
### 1. orderMatByFreq(D)
### 2. getRulesAsStrings(rm)
### 3. myHeatMap(mat) ### heatmap with row order preserved, columns clustered
### 4. getDuos(rules)


### 4. ls.nofun()
### 5. splitMatIntoList(mat,n.splits)
### 6. splitVectorIntoList(vec,n.splits)
### 8. matchMatrixNames(mat.1,mat.2)


# 
# library(gplots)
# library(RColorBrewer)

# 1.
# orderMatByFreq <- function(D){
#   freqs <- rowSums(D)
#   D <- D[order(freqs,decreasing=TRUE),]
#   return(D)
# }

# # 2. 
# getRulesAsStrings <- function(rm){
#   if(nrow(t(rm))==1){
#     return(paste0(sort(names(which(rm==1))),collapse = "+"))
#   }
#   rm.strings <- rep(NA,nrow(rm))
#   names(rm.strings) <- rownames(rm)
#   for(j in 1:length(rm.strings)){
#     rm.strings[j] <- paste0(sort(names(which(rm[j,]==1))),collapse = "+")
#   }
#   return(rm.strings)
# }

# 3) make heatmap, preserve row order
myHeatMap <- function(mat,rowv,colv,cexrow,my_palette){
  if(missing(my_palette)) my_palette <- colorRampPalette(c("grey","blue"))(n = 1000)
  if(missing(rowv)) rowv = FALSE
  if(missing(colv)) colv = TRUE
  if(missing(cexrow)) cexrow = 1
  heatmap.2(mat,Rowv=rowv,Colv = colv,trace="none",dendrogram="none",col=my_palette,
            key = FALSE, labCol = "",lhei = c(0.1,2),lwid=c(0.1,2),margins=c(2,10),cexRow = cexrow)
}


# ### 4. getDuos(rules)
# getDuos <- function(rules){
#   duos.list <- vector("list")
#   for(rule in rules){
#     events <- strsplit(rule,"\\+")[[1]]
#     temp <- t(combn(events,2,simplify = FALSE))
#     duos.list <- c(duos.list,temp)
#   }
#   duos <- unlist(lapply(duos.list,function(x)paste0(x,collapse="+")))
#   return(unique(duos))
# }
# 
# # 4) ls without functions ...
# ################################################
# # improved list of objects
# .ls.objects <- function (pos = 1, pattern, order.by,
#                          decreasing=FALSE, head=FALSE, n=5) {
#   napply <- function(names, fn) sapply(names, function(x)
#     fn(get(x, pos = pos)))
#   names <- ls(pos = pos, pattern = pattern)
#   obj.class <- napply(names, function(x) as.character(class(x))[1])
#   obj.mode <- napply(names, mode)
#   obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
#   obj.size <- napply(names, object.size)
#   obj.dim <- t(napply(names, function(x)
#     as.numeric(dim(x))[1:2]))
#   vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
#   obj.dim[vec, 1] <- napply(names, length)[vec]
#   out <- data.frame(obj.type, obj.size, obj.dim)
#   names(out) <- c("Type", "Size", "Rows", "Columns")
#   if (!missing(order.by))
#     out <- out[order(out[[order.by]], decreasing=decreasing), ]
#   if (head)
#     out <- head(out, n)
#   out
# }
# # shorthand
# lsos <- function(..., n=10) {
#   .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
# }
# 
# ls.nofun <- function(){
#   L <- lsos(n=Inf)
#   print(L[L$Type != "function",])
# }
################################################

########################################################################
### Parallel friendly functions
########################################################################
# 
# # 5)
# splitMatIntoList <- function(mat,n.splits){
#   ### Takes a matrix and splits into n.splits evenly sized matrices
#   n.size <- nrow(mat)
#   n.splits <- min(n.splits,n.size)
#   
#   base.size <- floor(n.size/n.splits)
#   remainder <- n.size%%n.splits
#   
#   list.of.mats <- vector("list",length=n.splits)
#   for(j in 1:n.splits){
#     if(j != n.splits) list.of.mats[[j]] <- mat[c((base.size*(j-1) + 1):(base.size*j)),]
#     if(j == n.splits) list.of.mats[[j]] <- mat[c((base.size*(j-1) + 1):n.size),]
#   }
#   names(list.of.mats) <- paste0("mat",1:n.splits)
#   return(list.of.mats)
# }
# 
# # 6)
# splitVectorIntoList <- function(vec,n.splits){
#   ### Takes a matrix and splits into n.splits evenly sized matrices
#   n.size <- length(vec)
#   n.splits <- min(n.splits,n.size)
#   
#   base.size <- floor(n.size/n.splits)
#   remainder <- n.size%%n.splits
#   
#   list.of.vecs <- vector("list",length=n.splits)
#   for(j in 1:n.splits){
#     if(j != n.splits) list.of.vecs[[j]] <- vec[c((base.size*(j-1) + 1):(base.size*j))]
#     if(j == n.splits) list.of.vecs[[j]] <- vec[c((base.size*(j-1) + 1):n.size)]
#   }
#   names(list.of.vecs) <- paste0("vec",1:n.splits)
#   return(list.of.vecs)
# }



# 8)
matchMatrixNames <- function(mat.1,mat.2){
  if(!identical(colnames(mat.1),colnames(mat.2)))return(FALSE)
  if(!identical(rownames(mat.1),rownames(mat.2)))return(FALSE)
  return(TRUE)
}




