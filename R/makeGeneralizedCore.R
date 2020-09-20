
### This file will contain functions for making GCRs, GCDs and GCEs

### 1. oneSubCoreIteration(D,Q,rm,tpl,til,subset.size,num.evaluated,msa)
### 2. makeSubCoreList(D,Q,rm,tpl,til,num.subsets,num.evaluated,msa,shouldPrint) - Exported
### 3. makeConfLevel(conf)
### 4. getGCRs(list.subsets.cores) - Exported
### 5. getDuosOneRule(rule)
### 6. getTriosOneRule(rule)
### 7. getDuos(rules) # Get unique duos from multiple rules
### 8. getTrios(rules) # Get unique duos from multiple rules

### 8. getGCDs(list.subsets.cores) - Exported
### 9. getGCTs(list.subsets.cores) - Exported
### 10. getGCEs(list.subsets.cores) - Exported



### 1.

#' @importFrom stats runif
oneSubCoreIteration <- function(D,P,rm,til,subset.size,num.evaluated,msa){
  # if(missing(subset.size)) subset.size <- 0.8
  # if(missing(num.evaluated)) num.evaluated <- 100
  msa.sub <- ceiling(subset.size*msa)
  msa.sub <- max(msa.sub,3) ### never less than 3
  idc <- sample.int(ncol(D),floor(ncol(D)*subset.size),replace = FALSE)
  D.sub <- D[,idc]
  P.sub <- P[,idc]
  rm.cm.sub <- makeRSCoverageMat(D.sub,rm)

  sub.rs.idc.list <- vector("list",length=length(til))
  names(sub.rs.idc.list) <- names(til)
  K <- 1
  sub.rs.idc.list[[K]] <- which.max(getSingleRuleWs(D.sub,rm,P.sub,msa.sub))
  mini.beg <- Sys.time()
  for(K in 2:length(til)){
    im <- til[[K]]
    if(nrow(t(im))!=1) im <- im[1:min(nrow(im),num.evaluated),]

    if(nrow(t(im))==1) im <- as.matrix(t(im))
    perfs <- evaluateIM(D.sub,rm,rm.cm.sub,im,P.sub,msa.sub)
    sub.rs.idc.list[[K]] <- im[which.max(perfs),]
  }
  core.cov.thresh <- runif(1,min=85,max=99)
  core.perf.thresh <- runif(1,min=85,max=99)

  ### Get best perfs and best covs
  best.sub.perfs <- rep(NA,length=length(sub.rs.idc.list))
  best.sub.covs <- best.sub.perfs
  for(K in 1:length(sub.rs.idc.list)){
    rs.idc <- sub.rs.idc.list[[K]]
    best.sub.perfs[K] <- getWofRS(D.sub,rm,rm.cm.sub,rs.idc,P.sub,msa.sub)
    if(K==1) best.sub.covs[K] <- mean(rm.cm.sub[rs.idc,])
    if(K>1) best.sub.covs[K] <- mean(colSums(rm.cm.sub[rs.idc,])>0)
  }

  best.perf.percents <- 100*best.sub.perfs/max(best.sub.perfs)
  best.covs.percents <- 100*best.sub.covs/max(best.sub.covs)
  ### Determine core K
  core.idc <- which(best.covs.percents>=core.cov.thresh)
  core.idc <- intersect(core.idc,which(best.perf.percents>=core.perf.thresh))
  core.K <- min(core.idc)
  if(length(core.idc)==0) {
    #print("no K satisfies core")
    core.K <- length(best.sub.perfs)
  }
  sub.core.rules <- getRulesAsStrings(rm[sub.rs.idc.list[[core.K]],])
  return(sub.core.rules)
}

### 2.
#' @title Get list of core rules from random subsets of samples
#'
#' @param D input matrix D
#' @param Q input matrix Q
#' @param rm binary rule matrix
#' @param til list of top rule set index matrices
#' @param num.subsets number of subset iterations, default is 100
#' @param num.evaluated number of top rs considered per k per iteration, default is 1000
#' @param msa msa
#' @param shouldPrint Print progress updates? Default is TRUE
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.full,
#'           pool.sizes=c(60,20,20),max.stored=10,msa=15,shouldPrint = FALSE)
#' subcore.list <- makeSubCoreList(D,Q,rm.full,til.p2,num.subsets=3,num.evaluated=10,msa=15)
#' @export
makeSubCoreList <- function(D,Q,rm,til,num.subsets,num.evaluated,msa,shouldPrint){
  if(missing(num.subsets)) num.subsets <- 100
  if(missing(num.evaluated)) num.evaluated <- 1000
  if(missing(shouldPrint)) shouldPrint <- TRUE
  subset.sizes <- runif(num.subsets,0.67,0.85)
  k.max <- length(til)
  list.subsets.cores <- vector("list",length=num.subsets)

  beg <- Sys.time()
  for(j in 1:num.subsets){
    list.subsets.cores[[j]] <- oneSubCoreIteration(D,Q,rm,til,subset.sizes[j],num.evaluated,msa)
    if(shouldPrint) {
      if(j%%5==0) print(paste0("Subset Core Iteration = ",j))
    }
  }
  if(shouldPrint) print(Sys.time()-beg)
  return(list.subsets.cores)
}


### 3. Decide confidence level thresholds in here
makeConfLevel <- function(conf){
  #med.thresh <- 40
  #high.thresh <- 80
  noncon.thresh <- 1
  con.thresh <- 51
  conf.level <- rep("Low",length(conf))
  #conf.level[which(conf>=med.thresh)] <- "Medium"#paste0("Intermediate: [",med.thresh,", ",high.thresh,")")
  #conf.level[which(conf>=high.thresh)] <- "Consensus GCR"#paste0("High: >= ",high.thresh)
  
  conf.level[which(conf>=noncon.thresh)] <- "Not Consensus"
  conf.level[which(conf>=con.thresh)] <- "Consensus"
  
  return(conf.level)
}

### 4
#' @title Get Generalized Core Rules
#'
#' @param list.subset.cores list of subset cores
#' @examples
#' list.subset.cores <- list(c("A.B.C","D.E","A.D"),c("A.C","B.C.D","D.E"),
#' c("A.B.C","D.E"),c("A.B.C","D.E","B.C.D"))
#' getGCRs(list.subset.cores) # Confidence column should be 100, 75, 50, 25, 25
#' @export
getGCRs <- function(list.subset.cores){
  conf <- 100*sort(table(unlist(list.subset.cores)))/length(list.subset.cores) ### gen core confidence level
  conf.level <- makeConfLevel(conf)
  rules <- names(conf) ### gen core rules
  ### Make df.r (r for rules)
  df.r <- data.frame(GCR=rules,Confidence=as.numeric(conf),Confidence.Level=conf.level)
  df.r$GCR <- factor(df.r$GCR, levels = df.r$GCR[order(df.r$Confidence)])
  df.r <- df.r[nrow(df.r):1,]
  return(df.r)
}



### 5. Get duos
#' Title
#'
#' @param rule rule
#'
#' @return
#' @export
#'
#' @examples
getDuosOneRule <- function(rule){
  all.duos <- c()
  events <- strsplit(rule,"\\.")[[1]]
  if(length(events)==2) all.duos <- c(all.duos, paste0(events[1],".",events[2]))
  if(length(events)>2){
    duos.mat <- t(combn(events,2))
    for(j in 1:nrow(duos.mat)) all.duos <- c(all.duos,paste0(duos.mat[j,],collapse = "."))
  }
  all.duos <- sort(all.duos)
  return(all.duos)
}

### 6. Get unique duos from multiple rules
#' Title
#'
#' @param rules rules
#'
#' @return
#' @export
#'
#' @examples
getDuos <- function(rules){
  all.duos <- c()
  for(j in 1:length(rules)){
    rule <- rules[j]
    all.duos <- c(all.duos,getDuosOneRule(rule))
  }
  all.duos <- sort(unique(all.duos))
  return(all.duos)
}

#' Title
#'
#' @param rule rule
#'
#' @return
#' @export
#'
#' @examples
getTriosOneRule <- function(rule){
  all.trios <- c()
  events <- strsplit(rule,"\\.")[[1]]
  if(length(events)==2) return(NULL)
  if(length(events)==3) return(rule)
  if(length(events)>3){
    trios.mat <- t(combn(events,3))
    for(j in 1:nrow(trios.mat)) all.trios <- c(all.trios,paste0(trios.mat[j,],collapse = "."))
  }
  all.trios <- sort(all.trios)
  return(all.trios)
}

#' Title
#'
#' @param rules rules
#'
#' @return
#' @export
#'
#' @examples
getTrios <- function(rules){
  all.trios <- c()
  for(j in 1:length(rules)){
    rule <- rules[j]
    all.trios <- c(all.trios,getTriosOneRule(rule))
  }
  all.trios <- sort(unique(all.trios))
  return(all.trios)
}

### 7.
#' @title Get Generalized Core Duos
#'
#' @param list.subset.cores list of subset cores
#' @examples
#' list.subset.cores <- list(c("A.B.C","D.E","A.D"),c("A.C","B.C.D","D.E"),
#' c("A.B.C","D.E"),c("A.B.C","D.E","B.C.D"))
#' getGCDs(list.subset.cores) # Confidence column should be 100, 100, 100, 75, 50, 25, 25
#' @export
getGCDs <- function(list.subset.cores){
  list.subset.duos <- list.subset.cores
  for(j in 1:length(list.subset.duos)) list.subset.duos[[j]] <- getDuos(list.subset.cores[[j]])
  temp <- getGCRs(list.subset.duos)
  colnames(temp)[1] <- "GCD"
  return(temp)
}


#' Title
#'
#' @param list.subset.cores list of subset cores
#'
#' @return
#' @export
#'
#' @examples
getGCTs <- function(list.subset.cores){
  list.subset.trios <- list.subset.cores
  #alljs <- c()
  for(j in 1:length(list.subset.trios)){
    trios <- getTrios(list.subset.cores[[j]])
    list.subset.trios[[j]] <- getTrios(list.subset.cores[[j]])
    if(is.null(trios)) list.subset.trios[[j]] <- ""
    #if("KCNN3-A.NOTCH2-A.NRAS-M" %in% list.subset.trios[[j]]) alljs <- c(alljs,j)
  } 
  # conf <- 100*sort(table(unlist(list.subset.trios)))/length(list.subset.trios) ### gen core confidence level
  # conf.level <- makeConfLevel(conf)
  # rules <- names(conf) ### gen core rules
  # ### Make df.r (r for rules)
  # df.r <- data.frame(GCR=rules,Confidence=as.numeric(conf),Confidence.Level=conf.level)
  # df.r$GCR <- factor(df.r$GCR, levels = df.r$GCR[order(df.r$Confidence)])
  # df.r <- df.r[nrow(df.r):1,]
  # 
  temp <- getGCRs(list.subset.trios)
  remove.idc <- which(temp$GCR=="")
  temp[setdiff(c(1:nrow(temp)),remove.idc),]
  colnames(temp)[1] <- "GCT"
  return(temp)
}

### 9.
#' @title Get Generalized Core Events
#'
#' @param list.subset.cores list of subset cores
#' @examples
#' list.subset.cores <- list(c("A.B.C","D.E","A.D"),
#' c("A.C","B.C.D","D.E"),c("A.B.C","D.E"),c("A.B.C","D.E","B.C.D"))
#' getGCEs(list.subset.cores) # Confidence column should be 100, 100, 100, 100, 100
#' @export
getGCEs <- function(list.subset.cores){
  list.subset.events <- list.subset.cores
  for(j in 1:length(list.subset.cores)) list.subset.events[[j]] <- unique(as.character(unlist(sapply(list.subset.cores[[j]],getRuleEvents))))
  temp <- getGCRs(list.subset.events)
  colnames(temp)[1] <- "GCE"
  return(temp)
}

# ### 5. Decide confidence level thresholds in here
# makeConfLevel <- function(conf,med.thresh,high.thresh){
#   if(missing(med.thresh)) med.thresh <- 40
#   if(missing(high.thresh)) high.thresh <- 80
#   
#   conf.level <- rep(paste0("Low: < ",med.thresh),length(conf)) 
#   conf.level[which(conf>=med.thresh)] <- paste0("Intermediate: [",med.thresh,", ",high.thresh,")")
#   conf.level[which(conf>=high.thresh)] <- paste0("High: >= ",high.thresh)
#   return(conf.level)
# }

### 6
#Title

#@param list.subsets.cores list.subsets.cores
#@param tissue tissue
#@param max.genes.shown max.genes.shown

#@return
#@export

#@examples



# makeDf.Rules <- function(list.subsets.cores,tissue,max.genes.shown){
#   if(missing(max.genes.shown)) max.genes.shown <- 4
#   conf <- 100*sort(table(unlist(list.subsets.cores)))/length(list.subsets.cores) ### gen core confidence level
#   conf.level <- makeConfLevel(conf)
#   rules <- names(conf) ### gen core rules
#   rules <- as.character(sapply(rules,function(x)makeRuleNice(x,tissue,max.genes.shown,sep=" + ")))
#   ### Make df.r (r for rules)
#   df.r <- data.frame(conf=as.numeric(conf),rules=rules,conf.level=conf.level)
#   df.r$rules <- factor(df.r$rules, levels = df.r$rules[order(df.r$conf)])
#   #df.r <- df.r[which(conf>=25),]
#   df.r <- df.r[nrow(df.r):1,]
#   return(df.r)
# }
#' # 
# ### 7 
# getConsensusRules <- function(tissue,min.conf.thresh,max.genes.shown){
#   if(missing(min.conf.thresh)) min.conf.thresh <- 40
#   if(missing(max.genes.shown)) max.genes.shown <- 2
#   
#   load(paste0(result.dir,tissue,"_filtered_results.RData"))
#   D <- filtered.results.list$D
#   rm.imp <- filtered.results.list$rm.imp
#   original.names <- rownames(D)
#   event.names <- makeEventsNice(original.names,tissue,max.genes.shown)
#   event.names <- gsub("\\(","<",event.names)
#   event.names <- gsub("\\)",">",event.names)
#   rownames(D) <- event.names
#   
#   colnames(rm.imp) <- event.names
#   #fm <- getFamMat(rm.imp)
#   #rm.cm.imp <- makeRSCoverageMat(D,rm.imp)
#   
#   load(paste0("~/Dropbox/MAKE.CRO.RESULTS/SUBSET.CORES/",tissue,"_subset_cores.RData"))
#   df.rules <- makeDf.Rules(list.subsets.cores,tissue,max.genes.shown)
#   
#   df.rules$rules <- gsub(" \\+ ",".",df.rules$rules)
#   df.rules$rules <- gsub("\\(","<",df.rules$rules)
#   df.rules$rules <- gsub("\\)",">",df.rules$rules)
#   is.valid <- c(1)
#   for(j in 2:nrow(df.rules)){
#     rs <- df.rules$rules[1:j]
#     rs.idc <- which(getRulesAsStrings(rm.imp)%in%rs)
#     #if(isRSvalid(rm.imp,rs.idc,fm)) is.valid <- c(is.valid,j)
#     
#     if(sum(getFamMat(rm.imp[rs.idc,]))==length(rs.idc)) is.valid <- c(is.valid,j)
#   }
#   consensus.idc <- intersect(is.valid,which(df.rules$conf>=min.conf.thresh))
#   return(df.rules$rules[consensus.idc])
# }
# 
# 
