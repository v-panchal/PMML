#####################################
# Gibbs Sampler
#####################################

require("Rcpp")
  

cppFunction(depends=c("RcppArmadillo","RcppDist"),' List pmmlogit(arma::mat x, arma::vec y, int S, 
  double sigmasq, int trunc, double wp, double a, double b, double omega, double omtype){
    // SOME USEFUL INFORMATION
            
            int n = x.n_rows;
            int p = x.n_cols;
            
            // INITIALIZE STORAGE
               arma::mat all_sumother(n,p-1);
            arma::mat res_beta(S,p);
            arma::mat res_invtau(S,p);
            arma::mat res_invdelta(S,p);
            arma::mat res_weights(S,p);
            arma::mat res_q(S,n);
            CharacterVector models(S);
            models(0)="0";
            
            arma::mat res_omega(S,p-1);
            arma::vec omvec = rep(omtype, S);
            
            // ASSIGN STARTING VALUES
            
            for(int k=0; k<p;k++){
            res_beta(0,k) = 0;
            res_invtau(0,k) = 1;
            res_weights(0,k)=0.5;
            }
            
            for(int k=0; k<p-1; k++){
            res_omega(0,k)=0.01;
            }
            
            for(int k=0; k<n;k++){
            res_q(0,k) = 1;
            }
            
           
            double b0= 0;
            double vb= 2;
            arma::vec kvec(n);
         
            for (int dd=0; dd<n; dd++){
            kvec[dd] = y[dd] - 0.5;
            }
            
        
            
            /////////////////////////////////////////////
            // GIBBS SAMPLER -- BIG LOOP
            
            for (int s=1; s < S; s++){
            arma::vec beta(p);
            arma::vec invTau(p);
            
            for (int pp=0; pp<p; pp++){
            beta(pp) = res_beta(s-1,pp);
            invTau(pp) = res_invtau(s-1,pp);
            }
            
            double wt;
            if(omvec(s)==1){
                            wt = omega;
            }else{ 
                            wt = r4beta(1,a+sum(res_beta.row(s-1)!=0),b+sum(res_beta.row(s-1)==0),wp,1)(0);
            }
            arma::vec weights = rep(wt, p-1);
          
            
            
            //Full conditional for q
            
            arma::vec q(n);
            for (int aa=0; aa <n; aa++){
            double xbeta = 0;
  
            for(int cc=0; cc<p; cc++){xbeta+=x(aa,cc)*beta(cc);}
            
            q[aa] = PG_6(xbeta);
            res_q(s,aa)=q(aa);
            
            }
            
            //  Full conditional for beta
              
              arma::vec bmean(p-1);
              arma::vec bvar(p-1);
              arma::vec sumother(n);
              arma::vec wone(p-1);
              arma::vec wtwo(p-1);
             
              
              std::string current;
            sumother=x*beta;
            for (int gg=0; gg<p-1; gg++){
 
            if (gg==0){
            sumother = sumother-x.col(gg+1)*beta(gg+1);
            
            }
            else{
            sumother= sumother-x.col(gg+1)*beta(gg+1)+x.col(gg)*beta(gg)   ;    
            }
            all_sumother.col(gg)=sumother;
          
                 double sum1=sum(kvec%x.col(gg+1));
                 double sum2=sum(q%sumother%x.col(gg+1));
                 double sum3=sum(q%pow(x.col(gg+1),2.0));
                 
               
                 
                 bmean(gg) = (sum1 - sum2 ) / (  invTau(gg+1)*(1.0/(sigmasq)) 
                 + sum3  );    
                 bvar(gg) = 1.0/(  invTau(gg+1)*(1.0/(sigmasq)) + sum3     );
                 wone(gg) = sqrt((invTau(gg+1))/((sigmasq)*((1.0/(sigmasq))*(invTau(gg+1)) 
                 + sum3)))  ;
                 wtwo(gg) = exp(0.5*  pow((sum1 - sum2),2.0) / 
                 (((invTau(gg+1)) * 1.0/(sigmasq)) + sum3 ) );
                 double sweights = 1- ((1-weights(gg))/(1-weights(gg) + weights(gg)*(wone(gg)*wtwo(gg)) )    );
                 res_weights(s,gg)=sweights;
                 res_omega(s,gg) = weights(gg);
                 double u = R::runif(0,1);
                 if(u < sweights) { beta(gg+1) = R::rnorm(bmean(gg),sqrt(bvar(gg))); 
                 
              current += "1";
               invTau(gg+1) = myrig( fabs(sqrt(sigmasq)/(beta(gg+1))) ,1);} 
               else {beta(gg+1) = 0; current += "0";}
              
                }  
              
              models(s)=current;
            // Full conditional for intercept
            
            invTau(0) = myrig( fabs(sqrt(sigmasq)/(beta(0))) ,1);
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
            beta(0) = R::rnorm(bomean,sqrt(bovar)) ;
            
          
            
            // SAVING VALUES
            
            for (int ll=0; ll<p-1; ll++){ 
            res_omega(s,ll) = weights(ll);
            }
            
            for (int ll=0; ll<p; ll++){
            res_beta(s,ll)=beta(ll);
            res_invtau(s,ll)=invTau(ll);
            }
            
           }
            
            List result;
            
            result["sb"] = res_beta;
            result["models"]=models;
            result["om"] = res_omega;
            return result;
            
            }  ', includes = c(' double myrig(double mu, double lambda){
                               double nu=R::rnorm(0,1);
                               double y= pow(nu, 2.0);
                               double x= mu + ((0.5*pow(mu,2.0)*y)/lambda) - ((0.5*mu)/lambda)*sqrt( (4*mu*lambda*y) + pow(mu*y,2.0) );
                               if (R::runif(0,1) > (mu/(mu+x))) {return x = pow(mu,2.0)/x ;}
                               else {return x;}
            }
                               ','
                                  double pinvg(double x, double mu, double lambda){
                                  return R::pnorm(sqrt(lambda/x)*(x/mu-1),0,1,1,0)+exp(2*lambda/mu)*R::pnorm(-sqrt(lambda/x)*(x/mu+1),0,1,1,0);
                               
                               } ','double seq_an(double x,int n){
                                 double Pi = 3.14159265358979323846;
                               double  t =0.64;
                               if ((0<x) && (x<=t)){
                               x = Pi*(n+.5)*pow(2/(Pi*x),1.5)*exp(-(2*pow(n+.5,2))/x);
                               }else{
                               x = Pi*(n+.5)*exp(-(Pi*Pi*pow(n+.5,2)*x)/2);
                               }
                               return(x);
                               
                               }','double test_un(double x){
                               double G= 0;
                               G=R::runif(0,x);
                               return(G);
                               }','double trun_igauss(double mu, double t){
                               double z=1/mu;  
                               double E=0;
                               double Ep=0;
                               double X= 1+t;
                               if (t>z){  
                               double alpha = 0;
                               
                               while (R::runif(0,1) > alpha){
                               E=R::rexp(1);
                               Ep=R::rexp(1);
                               while (E*E>2*Ep/t){
                               E=R::rexp(1);
                               Ep=R::rexp(1);
                               X=(1+t*E);
                               X= t/(X*X);
                               alpha=exp(-0.5*z*z*X);
                               }
                               }
                               }else {
                               double Y=0;
                               while (X>t){
                               Y=R::rnorm(0,1);
                               Y *= Y;
                               X = mu + .5*mu*mu*Y-.5*mu*sqrt(4*mu*Y+mu*Y*mu*Y);
                               
                               if (R::runif(0,1) > mu/(mu+X)){
                               X=(mu*mu)/X;
                               
                               }  
                               
                               } 
                               
                               
                               }
                               return(X);
                               }','double PG_6( double Z2){
  
                               double Z=fabs(Z2)/2; 
                               double t=0.64; 
                               double Pi= 3.14159265358979323846;
                               double K = (Pi*Pi)/8+(Z*Z)/2;
                               double p = Pi/(2*K)*exp((-K*t));
                               double q = 2*exp(-fabs(Z))*pinvg(t,1/Z,1);
                               double U = R::runif(0,1);
                               //double V = R::runif(0,1);
                               double X=0;
                               double mu= 1/Z;
                               if (U<p/(p+q)){
                               double E = R::rexp(1);
                               X = t + E/K; 
                               }else{
                               
                               X = trun_igauss(mu,t);
                               
                               }
                               
                               
                               double S = seq_an(X,0);  
                               double Y = R::runif(0,1)*S;
                               //int n = 0;
                               for (int n = 1; n++;) {
                               if (n % 2==1){
                               S = S-seq_an(X,n);
                               if (Y<=S){
                               X=X/4;
                               break;
                               }
                               }else{
                               S = S+seq_an(X,n);
                               if (Y>S){
                               break;
                               }
                               }
                               }
                               return X;
                               }
                               
                               ','#include <4beta.h>') ) 
  
