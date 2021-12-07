n <- 300
p <- 50
set.seed(123)
Sigma <- diag(p)
full <- matrix(c(rep(0.5, p*p)), ncol=p)
Sigma <- full + 0.5*Sigma
cholS <- chol(Sigma)
Beta <- c(-0.9,1.4,-1.8)
XX <- matrix(rnorm(n*p), ncol=p)
XX <- XX%*%cholS
beta <- numeric(p+1)
beta[c(1:length(Beta))] <- Beta
X <- cbind(rep(1,n),XX)
XB <- X%*%beta
probs <- as.vector(exp(XB)/(1+exp(XB)))
y <- rbinom(n,1,probs)


results <- pmmlogit(X, y, 5000, 10^3, 200, 0)




n <- 300
p <- 50
set.seed(123)
Sigma <- diag(p)
full <- matrix(c(rep(0.5, p*p)), ncol=p)
Sigma <- full + 0.5*Sigma
cholS <- chol(Sigma)
Beta <- c(-0.9,1.4,-1.8)
XX <- matrix(rnorm(n*p), ncol=p)
XX <- XX%*%cholS
beta <- numeric(p+1)
beta[c(1:length(Beta))] <- Beta
X <- cbind(rep(1,n),XX)
XB <- X%*%beta
probs <- as.vector(exp(XB)/(1+exp(XB)))
y <- rbinom(n,1,probs)

gibbsout <- pmmlogit(X, y, 5000, 10^3, 200, 0, 1, log(p), 0.05, 1)

finalresults <- process_gibbs(gibbsout, 100)

