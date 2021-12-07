process_gibbs <- function(gibbsout, n.burn) {

S <- nrow(gibbsout$sb)
p <- ncol(gibbsout$sb)
mygibbsout <- gibbsout$sb[(n.burn+2):S,]
mygibbsmodel <- gibbsout$models
mymodel <- names(sort(table(mygibbsmodel),decreasing=TRUE)[1])
varlist0 <- as.vector(as.numeric(unlist(strsplit(mymodel,""))))
if (length(varlist0)==1) {
  varlist <- rep(0,p)
} else {
  varlist <- varlist0
}

hpm <- which(varlist!=0)+1

nzperiter <- c()
for (i in 1:nrow(mygibbsout)){
  nzperiter[i] <- sum(mygibbsout[i,-1]!=0)
}
avgnzv <- mean(nzperiter)
avgzv <- (p-2) - avgnzv

m1 <- aperm(mygibbsout,c(2,1))

mbeta <- apply(m1,1,mean)
msbeta <- apply(m1,1,sd)
length(mbeta)
cvbeta <- abs(mbeta)/msbeta
cvbeta[cvbeta=="NaN"]=0
kk <- 2
mmm <- kmeans(cbind(abs(mbeta[-1]),cvbeta[-1]),kk,algorithm=c("Lloyd"),iter.max=100)
l <- c()
for (m in 1:kk){l[m] <- length(which(mmm$cluster==m))}
cind <- which(l<max(l))
ss <- list()
for (cc in 1:length(cind)){ss[[cc]] <- which(mmm$cluster==cind[cc])}

nzvar <- unlist(ss)+1

post_prob <- c()

for (i in 1:p){ post_prob[i] <- length(which(m1[i,]!=0))/length(m1[1,])}

return(list(selected_model=nzvar, Coef_est=mbeta, posterior_prob= post_prob, high.pr.model=hpm, avgnzv=avgnzv, avgzv= avgzv))

}
