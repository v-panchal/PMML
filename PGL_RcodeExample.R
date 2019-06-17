
#==============================================
# Required Packages
#==============================================

library(mvtnorm)
library(VGAM) #for rinv.gaussian
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(LaplacesDemon)

pmml <- function(x, y, S, sigma, trunc) {

cppFunction(depends=c("RcppArmadillo","RcppDist"),' List pglexactC(arma::mat x, arma::vec y, int S, double sigma, int trunc){
            // SOME USEFUL INFORMATION
            int n = x.n_rows;
            int p = x.n_cols;
            // INITIALIZE STORAGE
            arma::mat res_beta(S,p);
            arma::mat res_invtau(S,p);
            arma::mat res_invdelta(S,p);
            arma::mat res_weights(S,p);
            arma::vec res_gamma(S);
            arma::mat res_q(S,n);
            CharacterVector models(S);
            models(0)="0";
            // ASSIGN STARTING VALUES
            for(int k=0; k<p;k++){
            res_beta(0,k) = 0;
            res_invtau(0,k) = 1;
            res_weights(0,k)=0.5;
            }
            for(int k=0; k<n;k++){
            res_q(0,k) = 1;
            }
            res_gamma(0) = 1;
            double b0= 0;
            double vb= 2;
            arma::vec kvec(n);
            for (int dd=0; dd<n; dd++){
            kvec[dd] = y[dd] - 0.5;
            }
            //weights = rep(myw,p-1);
            // GIBBS SAMPLER -- BIG LOOP
            for (int s=1; s < S; s++){
            arma::vec beta(p);
            arma::vec invTau(p);
            for (int pp=0; pp<p; pp++){
            beta(pp) = res_beta(s-1,pp);
            invTau(pp) = res_invtau(s-1,pp);
            }
            double gamsq = res_gamma(s-1);
            arma::vec weights = rep(r4beta(1,(log(p)/2)+sum(res_beta.row(0)!=0),1+sum(res_beta.row(0)==0),0,1),p-1);
            
            //Full conditional for q
            arma::vec q(n);
            for (int aa=0; aa <n; aa++){
            double qsum =0;
            double xbeta = 0;
            for(int cc=0; cc<p; cc++){xbeta+=x(aa,cc)*beta(cc);}
            for (int bb=0; bb<trunc; bb++){
            qsum += (R::rgamma(1,1))/ ( (pow((bb+1)-0.5,2.0)*39.4784) + pow(xbeta,2.0));
            }
            q[aa] = 2*qsum;
            res_q(s,aa)=q(aa);
            }
            // Full conditional for beta
            arma::vec bmean(p-1);
            arma::vec bvar(p-1);
            arma::vec sumother(n);
            arma::vec wone(p-1);
            arma::vec wtwo(p-1);
            double gamsq_sigma=gamsq*sigma;
            std::string current;
            for (int gg=0; gg<p-1; gg++){
            for (int ee=0; ee<n; ee++){
            double temp=0;
            for (int ggg=0; ggg<p; ggg++){
            temp = temp+ x(ee,ggg)*beta(ggg);
            }
            sumother(ee) = temp - x(ee,gg+1)*beta(gg+1);
            }
            double sum1=sum(kvec%x.col(gg+1));
            double sum2=sum(q%sumother%x.col(gg+1));
            double sum3=sum(q%pow(x.col(gg+1),2.0));
            bmean(gg) = (sum1 - sum2 ) / (  invTau(gg+1)*(1.0/(gamsq_sigma)) 
            + sum3  );    
            bvar(gg) = 1.0/(  invTau(gg+1)*(1.0/(gamsq_sigma)) + sum3     );
            wone(gg) = sqrt((invTau(gg+1))/((gamsq_sigma)*((1.0/(gamsq_sigma))*(invTau(gg+1)) 
            + sum3)))  ;
            wtwo(gg) = exp(0.5*  pow((sum1 - sum2),2.0) / 
            (((invTau(gg+1)) * 1.0/(gamsq_sigma)) + sum3 ) );
            double sweights = 1- ((1-weights(gg))/(1-weights(gg) + weights(gg)*(wone(gg)*wtwo(gg)) )    );
            res_weights(s,gg)=sweights;
            double u = R::runif(0,1);
            if(u < sweights) { beta(gg+1) = R::rnorm(bmean(gg),sqrt(bvar(gg))); 
            current += "1";
            invTau(gg+1) = myrig( fabs(sqrt(gamsq_sigma)/(beta(gg+1))) ,1);} 
            else {beta(gg+1) = 0; current += "0";}	
            }
            models(s)=current;
            // Full conditional for intercept
            invTau(0) = myrig( fabs(sqrt(sigma*gamsq)/(beta(0))) ,1);
            arma::vec sotwo(n);
            for (int hh=0; hh<n; hh++){
            double temp2=0;
            for (int ii=0; ii<p; ii++){
            temp2 = temp2 + x(hh,ii)*beta(ii);
            }
            sotwo(hh) = temp2 - x(hh,0)*beta(0);
            }
            double bomean = (sum(kvec%x.col(0)) - sum(q%(sotwo%x.col(0)) )+ b0*(1.0/vb) )/( (1.0/vb) + sum(q%pow(x.col(0),2.0)) ) ;
            double bovar = 1.0/((1.0/vb) + sum(q%pow(x.col(0),2.0))) ;
            beta(0) = R::rnorm(bomean,bovar) ;
            // Full conditional for gamsq
            gamsq = myig( 3,accu(  pow(beta,2.0) / (2*sigma*(1.0/(invTau))) )+2 );
            res_gamma(s)=1;
            // SAVING VALUES
            for (int ll=0; ll<p; ll++){
            res_beta(s,ll)=beta(ll);
            res_invtau(s,ll)=invTau(ll);
            }
            }
            List result;
            result["sb"] = res_beta;
            result["gamma"] = res_gamma;
            result["tau"] = res_invtau;
            result["wts"] = res_weights;
            result["models"]=models;
            return result;
            }  ', includes = c(' double myrig(double mu, double lambda){
                               double nu=R::rnorm(0,1);
                               double y= pow(nu, 2.0);
                               double x= mu + ((0.5*pow(mu,2.0)*y)/lambda) - ((0.5*mu)/lambda)*sqrt( (4*mu*lambda*y) + pow(mu*y,2.0) );
                               if (R::runif(0,1) > (mu/(mu+x))) {return x = pow(mu,2.0)/x ;}
                               else {return x;}
            }
                               ', ' double myig(double shape, double scale){
                               return double (1.0/R::rgamma(shape, 1.0/scale));
                               } ','#include <4beta.h>') )
 
  mygibbsout0 <- pglexactC(X, y, S, sigma, trunc)
  mygibbsout <- mygibbsout0$sb
  mygibbsmodel <- mygibbsout0$models
  mymodel <- names(sort(table(mygibbsmodel),decreasing=TRUE)[1])
  varlist0 <- as.vector(as.numeric(unlist(strsplit(mymodel,""))))
  if (length(varlist0)==1) {
    varlist <- rep(0,p)
  } else {
    varlist <- varlist0
  }
  
  hpm <- which(varlist!=0)+1
  
  m1 <- aperm(mygibbsout,c(2,1))
  
  mbeta <- apply(m1,1,mean)
  msbeta <- apply(m1,1,sd)
  length(mbeta)
  cvbeta <- abs(mbeta)/msbeta
  cvbeta[cvbeta=="NaN"]=0
  kk <- 2
  mmm <- kmeans(cbind(abs(mbeta[-1]),cvbeta[-1]),kk,algorithm=c("Lloyd"),iter.max=1000)
  length(which(mmm$cluster==1))
  length(which(mmm$cluster==2))
  
  l <- c()
  for (m in 1:kk){l[m] <- length(which(mmm$cluster==m))}
  cind <- which(l<max(l))
  ss <- list()
  for (cc in 1:length(cind)){ss[[cc]] <- which(mmm$cluster==cind[cc])}
  sgenes <- c(unlist(ss)+1)
 
  post_prob <- c()
  
  for (i in 1:(p+1)){ post_prob[i] <- length(which(m1[i,]!=0))/S}
  
  return(list(selected_model=sgenes, Coef_est=mbeta, posterior_prob= post_prob, hpm=hpm))
   
}




