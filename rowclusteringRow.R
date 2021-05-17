########################################################################
#### Row Clustering with COLUMN THE SAME ,and COVARIATE term. #########
### set theta to have an array (i,r,k) ###

#load("data_cov.Rda")

###data frame to matrix form y.mat###
df2mat <- function(data,y,subject,question){
      row <- length(levels(subject))
      col<- length(levels(question))
      my.mat <- matrix(NA,row,col,byrow=T)
      for (i in 1:row) for (j in 1:col){
        leveli <- levels(subject)[i]
        levelj <- levels(question)[j]
        temp.df <- data[(subject==leveli)&(question==levelj),]
        if (length(temp.df$y)>0) my.mat[i,j] <- temp.df$y
      }
      return(my.mat)
    }

####Read y.mat###
    y.mat<-df2mat(data,data$y,as.factor(data$subject),as.factor(data$question))
    n<-nrow(y.mat)
    p<-ncol(y.mat)
    numberofcategory <- length(unique(c(data)))
    q <- numberofcategory

    numberofcluster <- 3 ## R
	 numpar <-  q + 3*numberofcluster - 3
	 set.seed(12345)
   	param <- runif(numpar,-1,1)
	 proportion <-  rep(1/numberofcluster,numberofcluster)
	# pomformula <- "Y~row + covariate"
    invect <- param
    RG <- numberofcluster
    pi.v <- proportion[1:RG]



####function 1: row clustering with column the same ###
    X1<-function(invect,pi.v,y.mat,RG){

      ##The log-likelihood in R, not used in optimisation##
      Rcluster.ll=function(theta,ppr.m,pi.v){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=numberofcategory
        theta[theta==0]=0.00001
        theta[theta<0]=0.00001
        pi.v[pi.v==0]=0.00001
        llc=0
        for(i in 1:n){
          for(j in 1:p){
            if(!is.na(y.mat[i,j])){
              for(r in 1:RG){
                llc=llc+ppr.m[i,r]*log(theta[i,r,y.mat[i,j]])
              }
            }
          }
        }
        for(i in 1:n){
          for(r in 1:RG){
            llc=llc+ppr.m[i,r]*log(pi.v[r])
          }
        }
        -llc
      }

      ###The imcomplete log-likelihhod, used in model selection###
      Rcluster.Incll=function(theta,pi.v)
      {
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=numberofcategory
        theta[theta==0]=0.00001
        theta[theta<0]=0.00001
        pi.v[pi.v==0]=0.00001
        logl = 0
        for(i in 1:n){
          sumoverR=0
          for(r in 1:RG){
            prodoverp=1
            for(j in 1:p)
              if(!is.na(y.mat[i,j])){prodoverp=prodoverp*theta[i,r,y.mat[i,j]]}
            sumoverR=sumoverR+pi.v[r]*prodoverp
          }
          logl=logl+log(sumoverR)
        }
        logl
      }


      ###rs Function ###
      POFM.rs=function(invect,ppr.m,pi.v){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=numberofcategory
        mu.in=(invect[1:(q-1)])
        #alpha.in=c(0,invect[(q):(q+RG-2)]) # set to zero constraints #
        #delta.in=c(invect[(q+RG-1):(q+2*RG-2)])
        alpha.in=c(invect[(q):(q+RG-1)]) # sum to zero constraints #
        alpha.in=c(alpha.in[1:(RG-1)],-sum(alpha.in[1:(RG-1)]))
	  	delta.in=c(invect[(q+RG):(q+2*RG-1)])
        this.theta=array(NA,c(n,RG,q))
	  	cov <- gender

	for(i in 1:n){
        for(r in 1:RG){
          this.theta[i,r,1]=exp(mu.in[1]-alpha.in[r]-delta.in[r]*cov[i])/(1+exp(mu.in[1]-alpha.in[r]-delta.in[r]*cov[i]))
        }}

	for (i in 1:n){
        for(r in 1:RG){
          for(k in 2:(q-1)){
            this.theta[i,r,k]=exp(mu.in[k]-alpha.in[r]-delta.in[r]*cov[i])/(1+exp(mu.in[k]-alpha.in[r]-delta.in[r]*cov[i]))-exp(mu.in[k-1]-alpha.in[r]-delta.in[r]*cov[i])/(1+exp(mu.in[k-1]-alpha.in[r]-delta.in[r]*cov[i]))
          }
        }}

	for(i in 1:n){
        for(r in 1:RG){
          this.theta[i,r,q]=1-sum(this.theta[i,r,1:(q-1)])
        }
        }


        this.theta[this.theta==0]=0.00001
        this.theta[this.theta<0]=0.00001
        pi.v[pi.v==0]=0.00001
        Rcluster.ll(this.theta,ppr.m,pi.v)
      }

      # Mian function used ##
      n=nrow(y.mat)
      p=ncol(y.mat)
      q=numberofcategory
      #plot(rep(0,RG),pi.v,xlim=c(0,500),ylim=c(0,1))
      ppr.m=matrix(NA,n,RG)
      theta.arr=array(1, c(n,RG,q))
      mu.in=(invect[1:(q-1)])
      #alpha.in=c(0,invect[(q):(q+RG-2)]) # set to zero constraints #
	  #delta.in=c(invect[(q+RG-1):(q+2*RG-2)])
	  alpha.in=c(invect[(q):(q+RG-1)]) # sum to zero constraints #
      alpha.in=c(alpha.in[1:(RG-1)],-sum(alpha.in[1:(RG-1)]))
	  delta.in=c(invect[(q+RG):(q+2*RG-1)])

	  cov <- gender

	for(i in 1:n){
      	for(r in 1:RG){
       	 theta.arr[i,r,1]=exp(mu.in[1]-alpha.in[r]-delta.in[r]*cov[i])/(1+exp(mu.in[1]-alpha.in[r]-delta.in[r]*cov[i]))
      }}

	for(i in 1:n){
      	for(r in 1:RG){
       	 for(k in 2:(q-1)){
          theta.arr[i,r,k]=exp(mu.in[k]-alpha.in[r]-delta.in[r]*cov[i])/(1+exp(mu.in[k]-alpha.in[r]-delta.in[r]*cov[i]))-exp(mu.in[k-1]-alpha.in[r]-delta.in[r]*cov[i])/(1+exp(mu.in[k-1]-alpha.in[r]-delta.in[r]*cov[i]))
       	 }
     		}
	  }


	for(i in 1:n){
      for(r in 1:RG){
        theta.arr[i,r,q]=1-sum(theta.arr[i,r,1:(q-1)])
      }
}
      outvect=invect
      # Run the EM cycle:
      iter=1

      while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>1e-04)))&(iter<500))
      {

        # E-step - Update posterior probabilities
        num.r=matrix(log(pi.v),n,RG,byrow=T)

        for(i in 1:n){
          for(r in 1:RG){
            for(j in 1:p){
              if(!is.na(y.mat[i,j])){
                num.r[i,r]=num.r[i,r]+log(theta.arr[i,r,y.mat[i,j]])
              }
            }
          }
        }
        for(i in 1:n) ppr.m[i,]=num.r[i,]-log(sum(exp(num.r[i,] + min(abs(num.r[i,]))))) + min(abs(num.r[i,]))

        ppr.m <- exp(ppr.m)

        pi.v=apply(ppr.m,2,mean)

        #point(rep(iter,RG),pi.v,pch=1,col="black")
        invect=outvect
        # M-step:
        #use numerical maximisation
        temp=optim(par=invect,fn=POFM.rs,
                   ppr.m=ppr.m,pi.v=pi.v,
                   method="L-BFGS-B",
                   hessian=F,control=list(maxit=10000))
        #print(temp)
        outvect=temp$par
        #print(abs(invect-outvect))
        #print(outvect)

        mu.out=(outvect[1:(q-1)])
       # alpha.out=c(0,outvect[(q):(q+RG-2)])
	   # delta.out=c(outvect[(q+RG-1):(q+2*RG-2)])
	   alpha.out=c(outvect[(q):(q+RG-1)])
	   alpha.out=c(alpha.out[1:(RG-1)],-sum(alpha.out[1:(RG-1)]))
	   delta.out=c(outvect[(q+RG):(q+2*RG-1)])

        theta.arr=array(1, c(n,RG,q))
	    cov = gender

	for (i in 1:n) {
        for(r in 1:RG){

          theta.arr[i,r,1]=exp(mu.out[1]-alpha.out[r]-delta.out[r]*cov[i])/(1+exp(mu.out[1]-alpha.out[r]-delta.out[r]*cov[i]))

        }}

	for (i in 1:n){
        for(r in 1:RG){
          for(k in 2:(q-1)){
            theta.arr[i,r,k]=exp(mu.out[k]-alpha.out[r]-delta.out[r]*cov[i])/(1+exp(mu.out[k]-alpha.out[r]-delta.out[r]*cov[i]))-exp(mu.out[k-1]-alpha.out[r]-delta.out[r]*cov[i])/(1+exp(mu.out[k-1]-alpha.out[r]-delta.out[r]*cov[i]))
          }}}

	for (i in 1:n){
        for(r in 1:RG){
         for(j in 1:p){
            theta.arr[i,r,q]=1-sum(theta.arr[i,r,1:(q-1)])
          }
        }
	}
        iter=iter+1
        if (floor(iter/5)==(iter/5))cat('iter=',iter, ' log.like=', temp$value ,'\n')
        #print(iter)
      }
      # Find cluster groupings:
      Rclus = vector("list",RG)
      for (rr in 1:RG) Rclus[[rr]] =
        (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
      # Save results:
      logl=Rcluster.Incll(theta.arr,pi.v)
      res.dev = -2*logl
      npar = numpar
      aic  = -2*logl + 2*npar
      aicc = aic + (2*npar*(npar+1)/(n*p - npar - 1))
      aic3 = -2*logl + 3*npar
      bic = -2*logl + npar*log(n*p)
     # icl = 2*temp$value + npar*log(n*p)
      out1 = round(c(n,p,q,logl,npar,aic,aicc,aic3,bic,RG),3)
      names(out1) = c("n","p","q","LogL","npar","AIC","AICc","AIC3","BIC","R")
      list("info"=out1,
           "pi"=round(pi.v,3),
           "theta"=round(theta.arr,3),
           "mu"= c(0,mu.out),
           "alpha"=alpha.out,
		   "delta"=delta.out,
           "ppr"=round(ppr.m,3),
           "Row Clusters"=Rclus)
    }

#sink("output_Row_R10.txt")
#X1(invect,pi.v,y.mat,RG)

