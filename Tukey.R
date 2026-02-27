 library(MASS)
 library(tmvtnorm)
 library(geoR)
 library(GIGrvg)
 library(truncnorm)
 library(FNN)
 library(LearnBayes)

#**************************************************************************
# Tukey transformation function: h parameter controls tail heaviness
#**************************************************************************

Tu<-function(w,h) w*exp(h*w^2/2)

#**************************************************************************
#                         **************************
#                         *  P R E D I C T I O N   *
#                         **************************
#**************************************************************************

# Function for spatial prediction at new locations
 predic<-function(n,n.p,DIST.total,ICORR,phi,alpha,h,sig,V,U){

     CORR.p<-exp(-DIST.total/phi)

     C.oo<-CORR.p[seq(1,n,1),seq(1,n,1)]
     C.pp<-CORR.p[seq(n+1,n+n.p,1),seq(n+1,n+n.p,1)]
     C.po<-CORR.p[seq(n+1,n+n.p,1),seq(1,n,1)]
	
     CIC<-C.po%*%ICORR
     Spp<-C.pp-CIC%*%t(C.po)
     U0<-as.vector(rmvnorm(1,rep(0,n.p)+as.vector(CIC%*%U),Spp))
     V0<-as.vector(rmvnorm(1,rep(0,n.p)+as.vector(CIC%*%V),Spp))

      return(alpha*abs(U0)+sig*Tu(V0,h))

 }

#**************************************************************************
#                   Generate from truncated inverse Gamma
#**************************************************************************

rtigamma<-function(s,r,L,U){
  g<-0
  while(g<L | g>U)   g<-rigamma(1,s,r)
  return(g)
}

#**************************************************************************
# Likelihood functions
#**************************************************************************

# Likelihood for a single observation
L<-function(ep,tau2,sig,alpha,h,U,V) exp(-.5/tau2*(ep-alpha*abs(U)-sig*Tu(V,h))^2)

# Log-likelihood (sum over all observations)
Ls<-function(ep,tau2,sig,alpha,h,U,V) -.5/tau2*sum( (ep-alpha*abs(U)-sig*Tu(V,h))^2 )


#**************************************************************************
#                 Compute conditional mean (br) and sd (dr) 
#                  for nearest neighbor Gaussian process
#**************************************************************************

b.d.i<-function(ii,phi,Neighbor.i,Dist){
     nu<-exp(-1/phi)
     Rsr.i<-nu^Dist[ii,Neighbor.i]
     Rrr.i<-nu^Dist[Neighbor.i,Neighbor.i]
     br<-solve(Rrr.i,Rsr.i)
     dr<-sqrt(1-sum(br*Rsr.i))
     return(list(br=br,dr=dr))
}

#==========================================================================

#              *******************************************            
#=======                ADDITIVE TUKEY MODEL                      =============
#              *******************************************            

#==========================================================================

Tukey<-function(simulate,burnin,Break,tun){

# Load data - replace with actual data path
Data<-read.csv("your_data.csv")  # Example: columns: x_coord, y_coord, response

coord<-Data[,1:2]
Y<-Data[,3]; n<-nrow(coord)

#**************************************************************************
#                         Set nearest neighbors
#**************************************************************************

 n.n<-  # Number of neighbors
 Neighbor.dis<-get.knnx(data=coord,query=coord,k=n.n+1)
 Neighbor<-Neighbor.dis$nn.index[, -1]  # Remove self-neighbor
 
#**************************************************************************
#                      Distance matrix 
#**************************************************************************

 DIST<-as.matrix(dist(coord,method="euclidean"))

#**************************************************************************
#      Initial values for parameters and latent variables; Hyperparameters
#**************************************************************************

# Initial phi (spatial decay parameter)
phi<-
m.p<-phi   # Prior mean for phi
c.p<-      # Prior scale for phi
Lp<-       # Lower bound for phi
Up<-       # Upper bound for phi

# Initial sigma (scale parameter)
sig<-
m.s<-sig   # Prior mean for sig
c.s<-      # Prior scale for sig
Lsi<-      # Lower bound for sig
Usi<-      # Upper bound for sig

# Initial nugget effect (error variance)
tau2<-
Lt<-       # Lower bound for tau2
Ut<-       # Upper bound for tau2

# Initial h (tail parameter)
h<-
Eh<-       # Prior mean for h
Vh<-       # Prior variance for h
Sh<-       # Grid for sampling h

# Initial alpha (skewness parameter)
alpha<-
m.a<-alpha # Prior mean for alpha
c.a<-      # Prior scale for alpha
La<-       # Lower bound for alpha
Ua<-       # Upper bound for alpha
Sa<-       # Grid for sampling alpha

# Initial beta (regression coefficient)
beta<-0
m.b<-beta  # Prior mean for beta
c.b<-      # Prior scale for beta


# Design matrix and linear predictor
X<-cbind(rep(1,n))
XB<-as.vector(X%*%beta)
ep<-Y-XB  # Residuals

# Initial latent variables
U<-rep(sqrt(2/pi),n)  # Skewness component (positive)
V<-rep(0,n)           # Heavy-tail component

#**************************************************************************
#                  Correlation matrix for latent variables
#**************************************************************************

  CORR<-exp(-DIST/phi)
  CHOL.C<-chol(CORR)
  ICORR<-chol2inv(CHOL.C)
  deter<-prod(diag(CHOL.C)^2)

  # Compute conditional moments for each location given neighbors
  br<-list(); dr<-c();  
  for(i in 1:n){
      br.dr<-b.d.i(i,phi,Neighbor[i,],DIST)
      br[[i]]<-br.dr$br;    dr[i]<-br.dr$dr
  }

#===========================================================================================

#==============    MCMC SIMULATION: UPDATE LATENT VARIABLES AND PARAMETERS    =============

#===========================================================================================

  # Storage for posterior samples
  D.tau<-c();  D.sig<-c(); D.beta<-c(); D.phi<-c(); D.h<-c(); D.alpha<-c();
  E.Y0<-c()

  R1<-0; sim<-0; ac<-0;  # Counters
  start_time <- Sys.time() 

  for(R in 1:simulate){

#**************************************************************************
#                      FULL CONDITIONAL FOR V AND U
#                   (Sampling using grid approximation)
#**************************************************************************

  # Update V for each location
  for(i in 1:n){
    m.v<-sum(br[[i]]*V[Neighbor[i,]])
    s.v<-dr[i]
    SV<-seq(m.v-20*s.v,m.v+20*s.v,.01)
    p.y<-L(ep[i],tau2,sig,alpha,h,U[i],SV)
    p.v<-dnorm(SV,m.v,s.v)
    pp0<-p.y*p.v
   if(all(pp0==0)){
    print(i); print("===========V=============");
   }else{
    V[i]<-sample(SV,1,prob=pp0/sum(pp0))
   }
  }
  
  # Update U for each location
  for(i in 1:n){
    m.u<-sum(br[[i]]*U[Neighbor[i,]])
    s.u<-dr[i]
    SU<-seq(m.u-20*s.u,m.u+20*s.u,.01)
    p.y<-L(ep[i],tau2,sig,alpha,h,SU,V[i])
    p.u<-dnorm(SU,m.u,s.u)
    pp0<-p.y*p.u
   if(all(pp0==0)){
    print(i); print("===========U=============");
   }else{
    U[i]<-sample(SU,1,prob=pp0/sum(pp0))
   }
  }

#**************************************************************************
#                      FULL CONDITIONAL FOR h
#                       (Grid-based sampling)
#**************************************************************************

   l.p.y<-c()
   for(j in 1:length(Sh)) l.p.y[j]<-Ls(ep,tau2,sig,alpha,Sh[j],U,V)
   p.h<-dtruncnorm(Sh,0,2,Eh,Vh)
   pp0<-exp(l.p.y-max(l.p.y))*p.h
   if(all(pp0==0)){
     print("===========h=============");    print(l.p.y)
   }else{
     h<-sample(Sh,1,prob=pp0/sum(pp0))
   }
   
#**************************************************************************
#                       FULL CONDITIONAL FOR alpha
#                         (Grid-based sampling)
#**************************************************************************
  
    l.p.y<-c()
    for(j in 1:length(Sa)) l.p.y[j]<-Ls(ep,tau2,sig,Sa[j],h,U,V)
    p.a<-dtruncnorm(Sa,a=La,b=Ua,m.a,c.a)

   pp0<-exp(l.p.y-max(l.p.y))*p.a
   if(all(pp0==0)){
     print("=========alpha===========");   print(l.p.y)
   }else{
     alpha<-sample(Sa,1,prob=pp0/sum(pp0))
   }

#**************************************************************************
#                       FULL CONDITIONAL FOR tau2
#                       (Inverse Gamma)
#**************************************************************************

   eta<-alpha*abs(U)+sig*Tu(V,h)

   yy<-.5*sum((ep-eta)^2)
   tau2<-rtigamma((n-1)/2,yy,Lt,Ut); 

#**************************************************************************
#                       FULL CONDITIONAL FOR beta
#                       (Normal)
#**************************************************************************
 
   v.b<-(n/tau2+1/c.b^2)^(-1)
   beta<-rnorm(1,v.b*(sum(Y-eta)/tau2+m.b/c.b^2),sqrt(v.b))

   XB<-as.vector(X%*%beta)
   ep<-Y-XB

#**************************************************************************
#                       FULL CONDITIONAL FOR sig
#                       (Truncated Normal)
#**************************************************************************

    TV<-Tu(V,h)
    v.s<-(sum(TV^2)/tau2+1/c.s^2)^(-1)
    sig<-rtruncnorm(1,Lsi,Usi,v.s*(sum(TV*( ep-alpha*abs(U) ))/tau2+m.s/c.s^2),sqrt(v.s))

#**************************************************************************
#                       METROPOLIS UPDATE FOR phi
#                       (Random walk proposal)
#**************************************************************************

   c.phi<-rnorm(1,phi,tun)   # tun=the tunning paramter
 if(c.phi>Lp & c.phi<Up){

  CORR.c<-exp(-DIST/c.phi)
  CHOL.C<-chol(CORR.c)
  ICORR.c<-chol2inv(CHOL.C)
  deter.c<-prod(diag(CHOL.C)^2) 

    # Likelihood ratio for both latent processes
    ratio.lik1<-((deter/deter.c)^.5)*exp(.5*t(V)%*%(ICORR-ICORR.c)%*%V)
    ratio.lik2<-((deter/deter.c)^.5)*exp(.5*t(U)%*%(ICORR-ICORR.c)%*%U)
    r.prior<-dtruncnorm(c.phi,a=Lp,b=Up,m.p,c.p)/dtruncnorm(phi,a=Lp,b=Up,m.p,c.p)

    r<-ratio.lik1*ratio.lik2*r.prior

      u<-runif(1,0,1)
      if(u<r){
	  phi<-c.phi
	  deter<-deter.c
	  CORR<-CORR.c
	  ICORR<-ICORR.c
        # Update conditional moments for new phi
        for(i in 1:n){
          br.dr<-b.d.i(i,phi,Neighbor[i,],DIST)
          br[[i]]<-br.dr$br;    dr[i]<-br.dr$dr
        }
        ac<-ac+1
      }
    }

#**************************************************************************
#                             SAVE POSTERIOR SAMPLES
#**************************************************************************

   if(R>burnin){
	  R1<-R1+1
      if(R1==Break){
	  sim<-sim+1		

#**************************************************************************
#                       Store simulated values
#**************************************************************************

  D.beta<-c(D.beta,beta)
  D.sig<-c(D.sig,sig)
  D.tau<-c(D.tau,tau2)
  D.phi<-c(D.phi,phi)
  D.h<-c(D.h,h)
  D.alpha<-c(D.alpha,alpha)

#****************************************************************************
  R1<-0
            }  #End if R1==break  
        }      #End if R>burnin

  }            #End for R.simulate

#****************************************************************************
#                         **************************
#                         *     POSTERIOR SUMMARY  *
#                         **************************
#****************************************************************************

print("E.beta - Posterior summary for beta")                                    
E.beta<-summary(D.beta);     print(E.beta);    
V.beta<-sd(D.beta);          print(V.beta); 

print("E.sig - Posterior summary for sig")                                    
E.sig<-summary(D.sig);  print(E.sig);    
V.sig<-sd(D.sig);       print(V.sig);

print("E.tau2 - Posterior summary for tau2")                                    
E.tau2<-summary(D.tau);   print(E.tau2);    
V.tau2<-sd(D.tau);        print(V.tau2);

print("E.phi - Posterior summary for phi")                                    
E.phi<-summary(D.phi);     print(E.phi);    
V.phi<-sd(D.phi);          print(V.phi); 

print("E.h - Posterior summary for h (tail parameter)")                                    
E.h<-summary(D.h);     print(E.h);    
V.h<-sd(D.h);          print(V.h); 

print("E.alpha - Posterior summary for alpha (skewness parameter)")                                    
E.alpha<-summary(D.alpha);     print(E.alpha);    
V.alpha<-sd(D.alpha);          print(V.alpha);

end_time <- Sys.time()
time<-end_time - start_time
cat("time====>",time,"\n")

print("accept - Acceptance rate for phi"); print(ac/simulate)                                    
  
}     # END FUNCTION