# CIBERSORT R script.
# TODO: 根据参考基因集解析基因表达矩阵中不同细胞亚型的相对比例.
# 
# Author: xiangjun.
###############################################################################
 
#' 该script主要包含三个函数,分别为CIBERSORT(),CoreAlg()和doPerm().
#' CIBERSORT(): 运行CIBERSORT的主要函数,嵌入后面两个函数,评估基因表达矩阵中的细胞比例,计算p值.
#' CoreAlg():   核心算法,用于计算样本中的细胞比例.
#' doPerm():    置换检验函数,采用蒙特卡罗采样方法计算p值.
#'
#' Input data:  mixture.file和sig.matrix,mixture.file格式如下,sig.matrix为官方提供的参考基因集文件.
#' Output data: 为包含所有结果的矩阵对象,共25列,分别为22种免疫细胞亚型在每个样本中的比例,p值,RMSE(均方根标准差)和R(相关系数).
#' 
#' 在R中的使用:
#'       source('CIBERSORT.R')
#'       TME.results <- CIBERSORT(mixture.file, sig.matrix, perm, QN)

#' 主要函数.
#' @param mixture.file 待分析的基因表达矩阵,为.txt文件,制表符分隔,第一列为基因名,第一行为样本名,名称不要有空格,会报错.
#' @param sig.matrix   参考基因矩阵,可使用的参考数据矩阵有LM22.txt,ABIS_seq和ABIS_microarray.txt.
#' @param perm         置换检验的次数,推荐perm=500.
#' @param QN           是否要进行分位数标准化(TRUE/FALSE).
CIBERSORT <- function(mixture.file, sig.matrix, perm=500, QN=TRUE){
  
  library(preprocessCore)
  #读入文件.
  X <- read.table(sig.matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- read.table(mixture.file, header=T, sep="\t", row.names=1,check.names=F)
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #排序.
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm ##置换检验次数.
  
  #如果mixture文件中最大值小于50，则去log.
  if(max(Y) < 50) {Y <- 2^Y}
  
  #mixture文件的分位数标准化.
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #两个表达矩阵的基因取交集.
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #对sig.matrix进行zscore转换.
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #相关系数的null distribution.
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
  
  print(nulldist)
  
  header <- c('Mixture',colnames(X),"P.value","Correlation","RMSE")
  print(header)
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #遍历基因表达矩阵.
  while(itor <= mixtures){
    
    y <- Y[,itor]
    
    #标准化样本数据集.
    y <- (y - mean(y)) / sd(y)
    
    #运行CIBERSORT核心算法.
    result <- CoreAlg(X, y)
    
    #获得结果.
    w <- result$w
    mix.r <- result$mix.r
    mix.rmse <- result$mix.rmse
    
    #计算p-value.
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix.r)) / length(nulldist))}
    
    #输出output
    out <- c(colnames(Y)[itor],w,pval,mix.r,mix.rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}
    
    itor <- itor + 1
    
  }
  
  #保存结果.
  #write.table(rbind(header,output), file="CIBERSORT.Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  
  #返回包含所有结果的矩阵对象.
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P.value","Correlation","RMSE")
  obj
}

#' 核心算法.
#' @param X 细胞特异性表达的基因.
#' @param y 每个样本的混合表达(mixed expression per sample).
CoreAlg <- function(X, y){
  
  library(preprocessCore)
  library(parallel)
  library(e1071)
  
  #尝试不同的nu值.
  svn.itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn.itor, res, mc.cores=1) else
    out <- mclapply(1:svn.itor, res, mc.cores=svn.itor)
  
  nusvm <- rep(0,svn.itor)
  corrv <- rep(0,svn.itor)
  
  #进行cibersort.
  t <- 1
  while(t <= svn.itor) {
    mySupportVectors <- out[[t]]$SV ##连续变量.
    myCoefficients <- out[[t]]$coefs ##系数定义.
    
    weights = t(myCoefficients) %*% mySupportVectors ##%*%为矩阵乘法.
    
    weights[which(weights<0)]<-0 ##设置权重.
    w<-weights/sum(weights) ##设置相关性.
    
    u <- sweep(X,MARGIN=2,w,'*') ##根据对应的权重与参考集相乘.
    
    k <- apply(u, 1, sum) ##统计每行总和.
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #选择最佳模型，从nus为0.25，0.5，0.75的3个模型里面挑选一个即可.
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #获得系数并将其标准化.
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q)) ##为计算后的22种免疫细胞的比例.
  
  mix.rmse <- rmses[mn]
  mix.r <- corrv[mn]
  
  #返回随机的y的免疫细胞组成情况，就是权重w.
  newList <- list("w" = w, "mix.rmse" = mix.rmse, "mix.r" = mix.r)
  
}

#' 置换检验函数.
#' @param X 细胞特异性表达的基因.
#' @param Y 每个样本的混合表达(mixed expression per sample).
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y)) ##从表达矩阵Y里面,随机挑选LM22矩阵基因数量的表达量值.
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor) ##打印进度.
    
    #random mixture.
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #标准化mixture,也就是scale函数.
    yr <- (yr - mean(yr)) / sd(yr)
    
    #运行CIBERSORT核心算法.
    result <- CoreAlg(X, yr)
    
    mix.r <- result$mix.r
    
    #store correlation.
    if(itor == 1) {dist <- mix.r}
    else {dist <- rbind(dist, mix.r)}
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist) ##计算p值.
}


########2022年9月12日15:22########
#' -----------测试案例-----------

# setwd("/pub6/Temp/xiangjun/R_exp/CIBERSORT")                                ## 设置工作路径.
# load("/pub6/Temp/xiangjun/R_exp/CIBERSORT/coadread_mRNA_tpm.RData")         ## 加载TPM格式的RNA-seq数据.
# library(tibble)
# mixture.file <- tibble::rownames_to_column(coadread_mRNA_tpm,"Gene_symbol") ## 把行名变成一列,名称不要有空格,会报错.
# write.table(mixture.file,file = "mixture.file.txt",row.names = F,quote = F,sep = "\t") ## mixture.file格式要正确,不能有NA或负值,tab分隔.
# source("/pub6/Temp/xiangjun/R_exp/CIBERSORT/CIBERSORT.R") 
# sig.matrix <- "LM22.txt"             ## 可以换为ABIS_seq和ABIS_microarray.
# TME.results <- CIBERSORT("mixture.file.txt", 
#                          sig.matrix, 
#                          perm = 500, ## perm是置换检验的次数.
#                          QN = F)     ## 芯片数据设置为T,测序数据设置为F. 
# TME.results

#' -----------结束测试-----------
#####################
