### ker.eq.R                   
### Implements the Kernel Method of Test Equating
###
### Copyright: Jorge Gonzalez, 2012.
### Last modification: 25-05-2012.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### Author contact information:
###
###      Jorge Gonzalez B.
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22 
###      Santiago
###      Chile
###      Voice: +56-2-3545467  URL  : http://www.mat.puc.cl/~jgonzale
###      Fax  : +56-2-3547729  Email: jgonzale@mat.puc.cl
###

ker.eq<-function(scores,kert,hx=NULL,hy=NULL,degree,design,Kp=1,scores2,degreeXA,degreeYA,J,K,L,wx,wy,w) 
UseMethod("ker.eq")

ker.eq.default<-function(scores,kert,hx=NULL,hy=NULL,degree,design,Kp=1,scores2,degreeXA,degreeYA,J,K,L,wx,wy,w)
{
	###########################
	#Call parameters
	###########################
	cl<-match.call()
	ni<-100
	
	#################################
	#Design specific data structure
	#################################
	if(design=="EG"){
	if(!is.matrix(scores)) stop("'scores' must be a matrix")
	if(dim(scores)[2]!=2) stop("'scores' must be a two column matrix")
		J<-K<-dim(scores)[1]
		xj=0:(J-1)
		yk=0:(K-1)	
		x=rep(xj,scores[,1])
		y=rep(yk,scores[,2])
				}
	else if(design=="SG"){
	if(!is.matrix(scores)) stop("'scores' must be a matrix")
	if(dim(scores)[2]<4) stop("'scores' must be a joint frequency matrix")
		J<-dim(scores)[1]
		K<-dim(scores)[2]
		xj=0:(J-1)
		yk=0:(K-1)
		x=rep(xj,apply(scores,1,sum))
		y=rep(yk,apply(scores,2,sum))
				}
	else if(design=="CB"){
		dat12<-scores
		dat21<-scores2
		J<-J
		K<-K
		N12<-dim(dat12)[1]
		N21<-dim(dat21)[1]
		xj<-0:(J-1)
		yk<-0:(K-1)
		
		x1=dat12[,1]
		x2=dat21[,1]
		y1=dat12[,2]
		y2=dat21[,2]
		y<-c(y1,y2)
		x<-c(x1,x2)
				}
	else if(design=="NEAT_CE" | design=="NEAT_PSE"){
		J<-J
		K<-K
		L<-L
		Np<-dim(scores)[1]
		Nq<-dim(scores2)[1]
		xj<-0:(J-1)
		yk<-0:(K-1)
		al<-0:(L-1)
		x<-scores[,1]
		ax<-scores[,2]
		y<-scores2[,1]
		ay<-scores2[,2]
				}
	Kp<-Kp

	##################################
	#loglinear presmoothing
	##################################
	if(design=="EG"){
		modelj<-loglin.smooth(scores[,1],degree[1],design) 
		modelk<-loglin.smooth(scores[,2],degree[2],design) 

	 	rj<-modelj$sp.est
 		sk<-modelk$sp.est
 			
		###########################
		#Forming C matrices
		###########################

		C_r=modelj$C
		C_s=modelk$C
				}
	else if(design=="SG"){
		model<-loglin.smooth(scores=scores,degree=degree,design=design)

		rj<-model$sp.est[,1]
		sk<-model$sp.est[,2]

		Cp<-model$C
				}
	else if(design=="CB"){
		model<-loglin.smooth(scores=scores,degree=degree,design=design,
		scores2=scores2,J=J,K=K,wx=wx,wy=wy)
		
		rj<-model$sp.est$rj
		sk<-model$sp.est$sk

		C<-model$C
				}
	else if(design=="NEAT_CE"){
		model<-loglin.smooth(scores=scores,degreeXA=degreeXA,degreeYA=degreeYA,
		design=design,scores2=scores2,K=K,J=J,L=L)

		rp<-rj<-model$sp.est$rp
		tp<-model$sp.est$tp
		tq<-model$sp.est$tq
		sq<-sk<-model$sp.est$sq

		C<-model$C
		D<-model$D
				}
	else if(design=="NEAT_PSE"){
		model<-loglin.smooth(scores=scores,degreeXA=degreeXA,degreeYA=degreeYA,
		design=design,scores2=scores2,K=K,J=J,L=L,w=w)

		rw<-rj<-model$sp.est$rw
		sw<-sk<-model$sp.est$sw

		C<-model$C
		D<-model$D
				}

	##################################
	#Summary statistics
	##################################

	### "EG" and "SG" designs ###

	if(design=="EG" | design=="SG"){
	mux.hat<-mean(x)
	muy.hat<-mean(y)
	
	sdx<-sd(x)
	sdy<-sd(y)
	
	nx<-length(x)
	ny<-length(y)

	skewx<-(sum((x-mux.hat)^3)/nx)/(sum((x-mux.hat)^2)/nx)^(3/2)
	skewy<-(sum((y-muy.hat)^3)/ny)/(sum((y-muy.hat)^2)/ny)^(3/2)

	kurtx<-nx*sum( (x-mux.hat)^4 )/(sum( (x-mux.hat)^2 )^2)
	kurty<-ny*sum( (y-muy.hat)^4 )/(sum( (y-muy.hat)^2 )^2)
	}

	### "CB" design ###
	else if(design=="CB"){
	mux.hat<-mean(x)
	muy.hat<-mean(y)
	
	sdx<-sd(x)
	sdy<-sd(y)

	mux1.hat<-mean(x1)
	muy1.hat<-mean(y1)
	mux2.hat<-mean(x2)
	muy2.hat<-mean(y2)
	
	sdx1<-sd(x1)
	sdy1<-sd(y1)
	sdx2<-sd(x2)
	sdy2<-sd(y2)
	
	nx1<-length(x1)
	ny1<-length(y1)
	nx2<-length(x2)
	ny2<-length(y2)

	skewx1<-(sum((x1-mux1.hat)^3)/nx1)/(sum((x1-mux1.hat)^2)/nx1)^(3/2)
	skewy1<-(sum((y1-muy1.hat)^3)/ny1)/(sum((y1-muy1.hat)^2)/ny1)^(3/2)
	skewx2<-(sum((x2-mux2.hat)^3)/nx2)/(sum((x2-mux2.hat)^2)/nx2)^(3/2)
	skewy2<-(sum((y2-muy2.hat)^3)/ny2)/(sum((y2-muy2.hat)^2)/ny2)^(3/2)

	kurtx1<-nx1*sum( (x1-mux1.hat)^4 )/(sum( (x1-mux1.hat)^2 )^2)
	kurty1<-ny1*sum( (y1-muy1.hat)^4 )/(sum( (y1-muy1.hat)^2 )^2)
	kurtx2<-nx2*sum( (x2-mux2.hat)^4 )/(sum( (x2-mux2.hat)^2 )^2)
	kurty2<-ny2*sum( (y2-muy2.hat)^4 )/(sum( (y2-muy2.hat)^2 )^2)
	}

	### "NEAT_CE" and "NEAT_PSE" designs ###
	else if(design=="NEAT_CE" | design=="NEAT_PSE"){
	mux.hat<-mean(x)
	muy.hat<-mean(y)
	
	sdx<-sd(x)
	sdy<-sd(y)
	
	nx<-length(x)
	ny<-length(y)

	skewx<-(sum((x-mux.hat)^3)/nx)/(sum((x-mux.hat)^2)/nx)^(3/2)
	skewy<-(sum((y-muy.hat)^3)/ny)/(sum((y-muy.hat)^2)/ny)^(3/2)

	kurtx<-nx*sum( (x-mux.hat)^4 )/(sum( (x-mux.hat)^2 )^2)
	kurty<-ny*sum( (y-muy.hat)^4 )/(sum( (y-muy.hat)^2 )^2)
	
	muax.hat<-mean(ax)
	muay.hat<-mean(ay)   

	sdax<-sd(ax)
	sday<-sd(ay)
		
	np<-length(x)
	nq<-length(y) 

	skewax<-(sum((ax-muax.hat)^3)/np)/(sum((ax-muax.hat)^2)/np)^(3/2)
	skeway<-(sum((y-muy.hat)^3)/nq)/(sum((ay-muay.hat)^2)/nq)^(3/2)

	kurtax<-np*sum( (ax-muax.hat)^4 )/(sum( (ax-muax.hat)^2 )^2)
	kurtay<-nq*sum( (ay-muay.hat)^4 )/(sum( (ay-muay.hat)^2 )^2)

	minx<-min(x)
	miny<-min(y)
	minax<-min(ax)
	minay<-min(ay)

	maxx<-max(x)
	maxy<-max(y)
	maxax<-max(ax)
	maxay<-max(ay)
	}

	#############################################################
	#Automatic bandwidth selection according to specific equating
        #design and kernel type
	###########################################################
	if(design=="EG"){
		if(is.null(hx) & is.null(hy)){
	  	h.x<-bandwidth(scores[,1],kert,degree[1],design)$h
	  	h.y<-bandwidth(scores[,2],kert,degree[2],design)$h
			                 }
	      else{ 
	  	h.x<-hx
	  	h.y<-hy
           		}
				}
	else if(design=="SG"){
		if(is.null(hx) & is.null(hy)){
		ban<-bandwidth(scores,kert,degree,design)
	  	h.x<-ban$hx
	  	h.y<-ban$hy
			                 }
	      else{ 
	  	h.x<-hx
	  	h.y<-hy
           		}
				}
	else if(design=="CB"){
		if(is.null(hx) & is.null(hy)){
		ban<-bandwidth(scores=scores,kert=kert,degree=degree,design=design,
			Kp=Kp,scores2=scores2,J=J,K=K,wx=wx,wy=wy) 

	  	h.x<-ban$hx
	  	h.y<-ban$hy
			                 }
	      else{ 
	  	h.x<-hx
	  	h.y<-hy
           		}
				}

	else if(design=="NEAT_CE"){
		if(is.null(hx) & is.null(hy)){
		ban<-bandwidth(scores=scores,kert=kert,degreeXA=degreeXA,
		degreeYA=degreeYA,design=design,Kp=Kp,scores2=scores2,J=J,K=K,
		L=L) 

	  	h.x<-ban$hx
	  	h.y<-ban$hy
		h.ap<-ban$hap
		h.aq<-ban$haq
			                 }
	      else{ 
	  	h.x<-hx
	  	h.y<-hy
           		}
				}

	else if(design=="NEAT_PSE"){
		if(is.null(hx) & is.null(hy)){
		ban<-bandwidth(scores=scores,kert=kert,degreeXA=degreeXA,
		degreeYA=degreeYA,design=design,Kp=Kp,scores2=scores2,J=J,K=K,
		L=L,w=w) 

	  	h.x<-ban$hx
	  	h.y<-ban$hy
							}
		else{ 
	  	h.x<-hx
	  	h.y<-hy
           		}
				}


		#############################
		#kernel specific a parameters
		#############################

		if(kert=="gauss"){
			a.x<-sqrt(var(x)/(var(x)+h.x^2))
      		a.y<-sqrt(var(y)/(var(y)+h.y^2))
      		if(design=="NEAT_CE"){
				a.p<-sqrt(var(ax)/(var(ax)+h.ap^2))
      			a.q<-sqrt(var(ay)/(var(ay)+h.aq^2))
				  	         }
					}
      		else if(kert=="logis"){     	
		      a.x<-sqrt(var(x)/(var(x)+(pi^2/3)*h.x^2))
     			a.y<-sqrt(var(y)/(var(y)+(pi^2/3)*h.y^2))
			if(design=="NEAT_CE"){
			      a.p<-sqrt(var(ax)/(var(ax)+(pi^2/3)*h.ap^2))
     				a.q<-sqrt(var(ay)/(var(ay)+(pi^2/3)*h.aq^2))
	     		                     }
					}
		else if(kert=="unif"){
			a.x<-sqrt(var(x)/(var(x)+(1/12)*h.x^2))
     			a.y<-sqrt(var(y)/(var(y)+(1/12)*h.y^2))
			if(design=="NEAT_CE"){
				a.p<-sqrt(var(ax)/(var(ax)+(1/12)*h.ap^2))
				a.q<-sqrt(var(ay)/(var(ay)+(1/12)*h.aq^2))
  	   				         }
					}

	################################################
	#Cumulative distribution functions (kernel type)
	################################################

	F.x.1=function(x,kert)
	{
		if(kert=="gauss"){
		aux=pnorm((x-a.x*xj-(1-a.x)*mux.hat)/(a.x*h.x))
            }
            else if(kert=="logis"){
		aux=plogis((x-a.x*xj-(1-a.x)*mux.hat)/(a.x*h.x))
            }
            else if(kert=="unif"){
            aux=punif((x-a.x*xj-(1-a.x)*mux.hat)/(a.x*h.x),-1/2,1/2)
            }
        aux2=rj*aux
	  return(aux2)
       }

	F.x=function(z,kert){sum(F.x.1(z,kert))}
	
	F.inv=function(u,kert)
	{
		a=function(z)F.x(z,kert)-u
		uniroot(a,c(-1,ni+1))$root  
	}


	G.y.1=function(y,kert)
	{
		if(kert=="gauss"){
		aux=pnorm((y-a.y*yk-(1-a.y)*muy.hat)/(a.y*h.y))
            }
            else if(kert=="logis"){
            aux=plogis((y-a.y*yk-(1-a.y)*muy.hat)/(a.y*h.y))
            }
            else if(kert=="unif"){
            aux=punif((y-a.y*yk-(1-a.y)*muy.hat)/(a.y*h.y),-1/2,1/2)
            }
            
	  aux2=sk*aux
	  return(aux2)
	}

	G.y=function(z,kert){sum(G.y.1(z,kert))}

	G.inv=function(u,kert)
	{
		a=function(z)G.y(z,kert)-u
		uniroot(a,c(-1,ni+1))$root  ## consequences of the interval -1:ni+1?
	}

	### More CDF for the NEAT_CE design ###
	if(design=="NEAT_CE"){
		H.p.1=function(ax,kert)
	{
		if(kert=="gauss"){
		aux=pnorm((ax-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap))
            }
            else if(kert=="logis"){
		aux=plogis((ax-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap))
            }
            else if(kert=="unif"){
            aux=punif((ax-a.p*al-(1-a.p)*muax.hat)/(a.p*h.ap),-1/2,1/2)
            }
        aux2=tp*aux
	  return(aux2)
       }

	H.p=function(z,kert){sum(H.p.1(z,kert))}
	
	Hp.inv=function(u,kert)
	{
		a=function(z)H.p(z,kert)-u
		uniroot(a,c(-1,ni+1))$root  	
	}



	H.q.1=function(ay,kert)
	{
		if(kert=="gauss"){
		aux=pnorm((ay-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq))
            }
            else if(kert=="logis"){
		aux=plogis((ay-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq))
            }
            else if(kert=="unif"){
            aux=punif((ay-a.q*al-(1-a.q)*muay.hat)/(a.q*h.aq),-1/2,1/2)
            }
        aux2=tq*aux
	  return(aux2)
       }

	H.q=function(z,kert){sum(H.q.1(z,kert))}
	
	Hq.inv=function(u,kert)
	{
		a=function(z)H.q(z,kert)-u
		uniroot(a,c(-1,ni+1))$root  	
	}

				}

	Fx=sapply(xj,F.x,kert=kert)
	G.i=sapply(Fx,G.inv,kert=kert)
	Gy=sapply(yk,G.y,kert=kert)
	F.i=sapply(Gy,F.inv,kert=kert)

	if(design=="NEAT_CE"){
	eqAx<-sapply(Fx,Hp.inv,kert=kert)
	alfa.p<-sapply(eqAx,H.q,kert=kert)
	G.iCE<-sapply(alfa.p,G.inv,kert=kert)

	eqAy<-sapply(Gy,Hq.inv,kert=kert)
	alfa.q<-sapply(eqAy,H.p,kert=kert)
	F.iCE<-sapply(alfa.q,F.inv,kert=kert)

	Hq=sapply(al,H.q,kert=kert)
	eqYa=sapply(Hq,G.inv,kert=kert)
				}

	#################################
	#Results (equated values)
	#################################

        eqYx=G.i
        eqXy=F.i

	if(design=="NEAT_CE"){
	eqCEYx<-G.iCE
	eqCEXy<-F.iCE				
				}
	
	##################################
	#Calculating SEE
	##################################

		############################
		#Functions derivatives
		############################

			################
			#Gaussian Kernel
			################

			if(kert=="gauss"){
			g=function(t)
			{
			 den=c()
   			for(i in 1:length(t)){
		         cc=sum((sk*dnorm((t[i]-a.y*yk-(1-a.y)*mean(y))/(a.y*h.y)))/(a.y*h.y))
	   		den[i]=cc
   	 		}
			return(den)
			}

			dGds=function(t)
			{
		     	d=matrix(0,nrow=length(t),ncol=length(yk))
			 for(i in 1:length(t)){
		   	d[i,]=pnorm((t[i]-a.y*yk-(1-a.y)*mean(y))/(a.y*h.y))-(((t[i]-mean(y))*(1-a.y^2)*((yk-mean(y))/sd(y))^2)/2+(1-a.y)*yk)*(sum(sk*dnorm((t[i]-a.y*yk-(1-a.y)*mean(y))/(a.y*h.y)))/(a.y*h.y)) 
	     		}
		   	return(d)
			}

			f=function(t)
			{
		   	den=c()
		   	for(i in 1:length(t)){
		     	cc=sum((rj*dnorm((t[i]-a.x*xj-(1-a.x)*mean(x))/(a.x*h.x)))/(a.x*h.x))
			   den[i]=cc
		    	}
			return(den)
			}

			dFdr=function(t){
		     	d=matrix(0,nrow=length(t),ncol=length(xj))
		     	for(i in 1:length(t)){
		   	d[i,]=pnorm((t[i]-a.x*xj-(1-a.x)*mean(x))/(a.x*h.x))-(((t[i]-mean(x))*(1-a.x^2)*((xj-mean(x))/sd(x))^2)/2+(1-a.x)*xj)*(sum(rj*dnorm((t[i]-a.x*xj-(1-a.x)*mean(x))/(a.x*h.x)))/(a.x*h.x))
			     }
		  	 return(d)
			}
				if(design=="NEAT_CE"){
				
				 hp=function(t){
				 den=c()
		   	       for(i in 1:length(t)){
				 cc=sum((tp*dnorm((t[i]-a.p*al-(1-a.p)*mean(ax))/(a.p*h.ap)))/(a.p*h.ap))
	                   den[i]=cc
   	 						}
					return(den)
						}
				
				 hq=function(t){
				 den=c()
		   	       for(i in 1:length(t)){
				 cc=sum((tq*dnorm((t[i]-a.q*al-(1-a.q)*mean(ay))/(a.q*h.aq)))/(a.q*h.aq))
	                   den[i]=cc
   	 						}
					return(den)
						}

				 dHpdtp=function(t){
			     	 d=matrix(0,nrow=length(t),ncol=length(al))
				 for(i in 1:length(t)){
			   	 d[i,]=pnorm((t[i]-a.p*al-(1-a.p)*mean(ax))/(a.p*h.ap))-(((t[i]-mean(ax))*(1-a.p^2)*((al-mean(ax))/sd(ax))^2)/2+(1-a.p)*al)*(sum(tp*dnorm((t[i]-a.p*al-(1-a.p)*mean(ax))/(a.p*h.ap)))/(a.p*h.ap)) 
		     		 }
			   	 return(d)
							}

				 dHqdtq=function(t){
			     	 d=matrix(0,nrow=length(t),ncol=length(al))
				 for(i in 1:length(t)){
			   	 d[i,]=pnorm((t[i]-a.q*al-(1-a.q)*mean(ay))/(a.q*h.aq))-(((t[i]-mean(ay))*(1-a.q^2)*((al-mean(ay))/sd(ay))^2)/2+(1-a.q)*al)*(sum(tq*dnorm((t[i]-a.q*al-(1-a.q)*mean(ay))/(a.q*h.aq)))/(a.q*h.aq)) 
		     		 }
			   	 return(d)
							}

							   }
}
	else if(kert=="logis"){

			################
			#Logistic Kernel
			################

	g=function(t)
	{
	 den=c()
   		for(i in 1:length(t)){
         cc=sum((sk*dlogis((t[i]-a.y*yk-(1-a.y)*mean(y))/(a.y*h.y)))/(a.y*h.y))
	   den[i]=cc
   	 }
	return(den)
	}

	dGds=function(t)
	{
     	d=matrix(0,nrow=length(t),ncol=length(yk))
	 for(i in 1:length(t)){
   	d[i,]=plogis((t[i]-a.y*yk-(1-a.y)*mean(y))/(a.y*h.y))-(((t[i]-mean(y))*(1-a.y^2)*((yk-mean(y))/sd(y))^2)/2+(1-a.y)*yk)*(sum(sk*dlogis((t[i]-a.y*yk-(1-a.y)*mean(y))/(a.y*h.y)))/(a.y*h.y)) 
	     }
   	return(d)
	}

	f=function(t)
	{
   	den=c()
   	for(i in 1:length(t)){
     	cc=sum((rj*dlogis((t[i]-a.x*xj-(1-a.x)*mean(x))/(a.x*h.x)))/(a.x*h.x))
	   den[i]=cc
    	}
	return(den)
	}

	dFdr=function(t){
     	d=matrix(0,nrow=length(t),ncol=length(xj))
     	for(i in 1:length(t)){
   	d[i,]=plogis((t[i]-a.x*xj-(1-a.x)*mean(x))/(a.x*h.x))-(((t[i]-mean(x))*(1-a.x^2)*((xj-mean(x))/sd(x))^2)/2+(1-a.x)*xj)*(sum(rj*dlogis((t[i]-a.x*xj-(1-a.x)*mean(x))/(a.x*h.x)))/(a.x*h.x))
	     }
  	 return(d)
	}
	}

			################
			#Uniform Kernel
			################

	else if(kert=="unif"){
	g=function(t)
	{
	 den=c()
   		for(i in 1:length(t)){
         cc=sum((sk*dunif((t[i]-a.y*yk-(1-a.y)*mean(y))/(a.y*h.y),-1/2,1/2))/(a.y*h.y))
	   den[i]=cc
   	 }
	return(den)
	}

	dGds=function(t)
	{
     	d=matrix(0,nrow=length(t),ncol=length(yk))
	 for(i in 1:length(t)){
      d[i,]=punif((t[i]-a.y*yk-(1-a.y)*mean(y))/(a.y*h.y),-1/2,1/2)-(((t[i]-mean(y))*(1-a.y^2)*((yk-mean(y))/sd(y))^2)/2+(1-a.y)*yk)*(sum(sk*dunif((t[i]-a.y*yk-(1-a.y)*mean(y))/(a.y*h.y),-1/2,1/2))/(a.y*h.y))
	     }
   	return(d)
	}

	f=function(t)
	{
   	den=c()
   	for(i in 1:length(t)){
     	cc=sum((rj*dunif((t[i]-a.x*xj-(1-a.x)*mean(x))/(a.x*h.x),-1/2,1/2))/(a.x*h.x))
	   den[i]=cc
    	}
	return(den)
	}

	dFdr=function(t){
     	d=matrix(0,nrow=length(t),ncol=length(xj))
     	for(i in 1:length(t)){
   	d[i,]=punif((t[i]-a.x*xj-(1-a.x)*mean(x))/(a.x*h.x),-1/2,1/2)-(((t[i]-mean(x))*(1-a.x^2)*((xj-mean(x))/sd(x))^2)/2+(1-a.x)*xj)*(sum(rj*dunif((t[i]-a.x*xj-(1-a.x)*mean(x))/(a.x*h.x)),-1/2,1/2)/(a.x*h.x))
	     }
  	 return(d)
	}
	}


	###################
	#Forming SEE vector

	Jey<-function(val) (1/g(eqYx)[val+1])*c(   dFdr(val)         ,-1*dGds(eqYx)[val+1,])
        Jex<-function(val) (1/f(eqXy)[val+1])*c(-1*dFdr(eqXy)[val+1,],   dGds(val)         )

	if(design=="NEAT_CE"){
	JeyCE<-function(val) c(hq(eqAx)[val+1]/g(eqCEYx)[val+1]*(1/hp(eqAx)[val+1])*c(dFdr(val),-1*dHpdtp(eqAx)[val+1,]),(1/g(eqAx)[val+1])*c(dHqdtq(val),-1*dGds(eqAx)[val+1,]))
				}

	####################################
	#Setting design-specific parameters
	####################################

	if(design=="EG"){
	Tr=degree[1]
	Ts=degree[2]
	T=Tr+Ts
	drdR=diag(J)
	dsdS=diag(K)
	CS=C_s
	CR=C_r
	JDF<-adiag(drdR,dsdS)
        C<-adiag(CR,CS)
	}

	else if(design=="SG"){
	Tp=dim(Cp)[2]
	T=Tp	
	M<-do.call(cbind, rep(list(diag(J)), K))
	N<-do.call(adiag,rep(list(t(rep(1,J))),K))
	JDF<-rbind(M,N)
        C<-Cp
	}

	else if(design=="CB"){
	T<-dim(C)[2]
	M<-do.call(cbind, rep(list(diag(J)), K))
	N<-do.call(adiag,rep(list(t(rep(1,J))),K))
	JDF<-rbind(cbind(wx*M,(1-wx)*M),cbind((1-wy)*N,wy*N))
	C<-C
	}	
	else if(design=="NEAT_CE"){
	T<-dim(C)[2]
	MP<-do.call(cbind,rep(list(diag(J)), L))
	NP<-do.call(adiag,rep(list(t(rep(1,J))),L))
	MQ<-do.call(cbind,rep(list(diag(K)), L))
	NQ<-do.call(adiag,rep(list(t(rep(1,K))),L))
	JDF<-adiag(rbind(MP,NP),rbind(NQ,MQ))
        C<-C
	}
#	else if(design=="NEAT_PSE"){
	T<-dim(C)[2]
#	jdf11<-
#	jdf12<-
#	jdf21<-
#	jdf22<-
#	JDF<-rbind(cbind(jfd11,jdf12),cbind(jdf21.jdf22))
#	JDF<-matrix(runif((J*L+K*L)*(J*L+K*L)),ncol=(J*L+K*L))
#	C<-C
#	}

	###########################################
	#SEE-vector and calculating SEEYx and SEEXy
	###########################################

	if(design=="NEAT_CE"){
	sevecYx<-matrix(0,nrow=J,ncol=T)
	SEEYx<-c()
	for(i in 0:(J-1)){
	sevecYx[i+1,]<-(JeyCE(i)%*%JDF%*%C)
	SEEYx[i+1]<-sqrt(sum((JeyCE(i)%*%JDF%*%C)^2))
				}
	}
	else if(design=="CB"){
	sevecYx<-matrix(0,nrow=J,ncol=T)
        sevecXy<-matrix(0,nrow=K,ncol=T)
	SEEYx<-c()
	SEEXy<-c()
	for(i in 0:(J-1)){
	sevecYx[i+1,]<-(Jey(i)%*%JDF%*%C)
	SEEYx[i+1]<-sqrt(sum((Jey(i)%*%JDF%*%C)^2))
				}
	for(i in 0:(K-1)){
	sevecXy[i+1,]<-(Jex(i)%*%JDF%*%C) 
	SEEXy[i+1]<-sqrt(sum((Jex(i)%*%JDF%*%C)^2))
				   }
	}
	else if(design=="NEAT_PSE"){
	SEEYx<-runif(J)
	SEEXy<-runif(K)
	sevecYx<-matrix(runif(J*T),ncol=T)
	sevecXy<-matrix(runif(K*T),ncol=T)
	}
	else{
	sevecYx<-matrix(0,nrow=J,ncol=T)
        sevecXy<-matrix(0,nrow=K,ncol=T)
	SEEYx<-c()
	SEEXy<-c()
	for(i in 0:(J-1)){
	sevecYx[i+1,]<-(Jey(i)%*%JDF%*%C)
	sevecXy[i+1,]<-(Jex(i)%*%JDF%*%C) 
	SEEXy[i+1]<-sqrt(sum((Jex(i)%*%JDF%*%C)^2))
	SEEYx[i+1]<-sqrt(sum((Jey(i)%*%JDF%*%C)^2))
				}
	}

	#############
	#Output
	#############
	if(design=="EG" | design=="SG"){
	res<-list(call=cl,kert=kert,design=design,eqYx=eqYx,eqXy=eqXy,h.x=h.x,
	h.y=h.y,SEEYx=SEEYx,SEEXy=SEEXy,sevecYx=sevecYx,sevecXy=sevecXy,
	score=0:(J-1),rj=rj,sk=sk,nx=nx,ny=ny,meanx=mux.hat,meany=muy.hat,
	sdx=sdx,sdy=sdy,kurtx=kurtx,kurty=kurty,skewx=skewx,skewy=skewy)
	}
	else if(design=="CB"){
	res<-list(call=cl,kert=kert,design=design,eqYx=eqYx,eqXy=eqXy,h.x=h.x,
	h.y=h.y,SEEYx=SEEYx,SEEXy=SEEXy,sevecYx=sevecYx,sevecXy=sevecXy,
	score.x=0:(J-1),score.y=0:(K-1),rj=rj,sk=sk,meanx1=mux1.hat,
	meany1=muy1.hat,meanx2=mux2.hat,meany2=muy2.hat,sdx1=sdx1,sdy1=sdy1,
	sdx2=sdx2,sdy2=sdy2,N12=N12,N21=N21,skewx1=skewx1,skewy1=skewy1,
	skewx2=skewx2,skewy2=skewy2,kurtx1=kurtx1,kurty1=kurty1,kurtx2=kurtx2,
	kurty2=kurty2)
	}
	else if(design=="NEAT_CE"){
	res<-list(call=cl,kert=kert,design=design,eqCEYx=eqCEYx,h.x=h.x,h.y=h.y,
	h.ap=h.ap,h.aq=h.aq,SEEYx=SEEYx,sevecYx=sevecYx,score.x=0:(J-1),
	score.y=0:(K-1),score.a=0:(L-1),rp=rp,sq=sq,tp=tp,tq=tq,np=np,nq=nq,
	meanx=mux.hat,meany=muy.hat,meanax=muax.hat,meanay=muay.hat,sdx=sdx,
	sdy=sdy,sdax=sdax,sday=sday,kurtx=kurtx,kurty=kurty,kurtax=kurtax,
	kurtay=kurtay,skewx=skewx,skewy=skewy,skewax=skewax,skeway=skeway,
	minx=minx,miny=miny,minax=minax,minay=minay,maxx=maxx,maxy=maxy,
	maxax=maxax,maxay=maxay)
					}
	else if(design=="NEAT_PSE"){
	res<-list(call=cl,kert=kert,design=design,eqYx=eqYx,eqXy=eqXy,h.x=h.x,
	h.y=h.y,SEEYx=SEEYx,SEEXy=SEEXy,sevecYx=sevecYx,sevecXy=sevecXy,
	score.x=0:(J-1),score.y=0:(K-1),score.a=0:(L-1),rj=rj,sk=sk,np=np,nq=nq,
	meanx=mux.hat,meany=muy.hat,meanax=muax.hat,meanay=muay.hat,sdx=sdx,
	sdy=sdy,sdax=sdax,sday=sday,kurtx=kurtx,kurty=kurty,kurtax=kurtax,
	kurtay=kurtay,skewx=skewx,skewy=skewy,skewax=skewax,skeway=skeway,
	minx=minx,miny=miny,minax=minax,minay=minay,maxx=maxx,maxy=maxy,
	maxax=maxax,maxay=maxay)
					}
       class(res)<-"ker.eq"
		   return(res)
}


print.ker.eq<-function(x,...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\nEquated values under the",x$design,"design:\n")	
	cat("\n")
	if(x$design=="EG" | x$design=="SG"){	
	print(data.frame(Score=x$score,eqYx=x$eqYx,eqXy=x$eqXy))
	}
	else if(x$design=="CB"){
	print(data.frame(Score_X=x$score.x,eqYx=x$eqYx))
	cat("\n")
	print(data.frame(Score_Y=x$score.y,eqXy=x$eqXy))
	cat("\n")
	}
	else if(x$design=="NEAT_CE"){
	print(data.frame(Score=x$score.x,eqCEYx=x$eqCEYx))
	cat("\n")
	}
	else if(x$design=="NEAT_PSE"){
	print(data.frame(Score=x$score.x,eqYx=x$eqYx,eqXy=x$eqXy))
	cat("\n")
	}
}	


summary.ker.eq<-function(object,...) 
{
	if(object$design=="EG" | object$design=="SG"){	

	descriptives<-rbind(c(object$nx,object$ny),
                          c(object$meanx,object$meany),
				  c(object$sdx,object$sdy),
				  c(object$skewx,object$skewy),
				  c(object$kurtx,object$kurty)) 
	dimnames(descriptives)<-list(c("Total","Mean","SD","Skewness","Kurtosis"),
					     c("X","Y"))
	kert<-object$kert
	design<-object$design
	SEEYx<-object$SEEYx
	SEEXy<-object$SEEXy
	equatedVal<-data.frame(Score=object$score,eqYx=object$eqYx,eqXy=
				     object$eqXy,SEEYx=SEEYx,SEEXy=SEEXy)
	bandwidthVal<-data.frame(hx=object$h.x,hy=object$h.y)
	res<-list(call=object$call,equatedVal=equatedVal,
		    bandwidthVal=bandwidthVal,descriptives=round(descriptives,4),
		    design=design,kert=kert)
	}
	else if(object$design=="CB"){
	descriptives12<-rbind(c(object$N12,object$N21),
                          c(object$meanx1,object$meany2),
				  c(object$sdx1,object$sdy2),
				  c(object$skewx1,object$skewy2),
				  c(object$kurtx1,object$kurty2))
	descriptives21<-rbind(c(object$N21,object$N12),
                          c(object$meanx2,object$meany1),
				  c(object$sdx2,object$sdy1),
				  c(object$skewx2,object$skewy1),
				  c(object$kurtx2,object$kurty1))
 
	dimnames(descriptives12)<-list(c("Total","Mean","SD","Skewness","Kurtosis"),
					     c("X1","Y2"))
	dimnames(descriptives21)<-list(c("Total","Mean","SD","Skewness","Kurtosis"),
					     c("X2","Y1"))

	kert<-object$kert
	design<-object$design
	SEEYx<-object$SEEYx
	SEEXy<-object$SEEXy
	equatedValYx<-data.frame(Score=object$score.x,eqYx=object$eqYx,SEEYx=SEEYx)
	equatedValXy<-data.frame(Score=object$score.y,eqXy=object$eqXy,SEEXy=SEEXy)

	bandwidthVal<-data.frame(hx=object$h.x,hy=object$h.y)
	res<-list(call=object$call,equatedValYx=equatedValYx,equatedValXy=equatedValXy,
		    bandwidthVal=bandwidthVal,descriptives12=round(descriptives12,4),
		    descriptives21=round(descriptives21,4),design=design,kert=kert)
	}
	else if(object$design=="NEAT_CE"){
	descriptivesP<-rbind(c(object$np,object$np),
                          c(object$meanx,object$meanax),
				  c(object$sdx,object$sdax),
				  c(object$skewx,object$skewax),
				  c(object$kurtx,object$kurtax),
				  c(object$minx,object$minax),
				  c(object$maxx,object$maxax))
	descriptivesQ<-rbind(c(object$nq,object$nq),
                          c(object$meany,object$meanay),
				  c(object$sdy,object$sday),
				  c(object$skewy,object$skeway),
				  c(object$kurty,object$kurtay),
				  c(object$miny,object$minay),
				  c(object$maxy,object$maxay))

	dimnames(descriptivesP)<-list(c("Total","Mean","SD","Skewness","Kurtosis",
						  "Min","Max"), c("X","A"))
	dimnames(descriptivesQ)<-list(c("Total","Mean","SD","Skewness","Kurtosis",
					        "Min","Max"), c("Y","A"))

	kert<-object$kert
	design<-object$design
	SEEYx<-object$SEEYx
	equatedVal<-data.frame(Score=object$score.x,eqCEYx=object$eqCEYx,SEEYx=SEEYx)

	bandwidthVal<-data.frame(hx=object$h.x,hy=object$h.y,h_AX=object$h.ap,
					 h_AY=object$h.aq)
	res<-list(call=object$call,equatedVal=equatedVal,bandwidthVal=bandwidthVal,
		    descriptivesP=round(descriptivesP,4),
		    descriptivesQ=round(descriptivesQ,4),design=design,kert=kert)
	}
	else if(object$design=="NEAT_PSE"){
	descriptivesP<-rbind(c(object$np,object$np),
                          c(object$meanx,object$meanax),
				  c(object$sdx,object$sdax),
				  c(object$skewx,object$skewax),
				  c(object$kurtx,object$kurtax),
				  c(object$minx,object$minax),
				  c(object$maxx,object$maxax))
	descriptivesQ<-rbind(c(object$nq,object$nq),
                          c(object$meany,object$meanay),
				  c(object$sdy,object$sday),
				  c(object$skewy,object$skeway),
				  c(object$kurty,object$kurtay),
				  c(object$miny,object$minay),
				  c(object$maxy,object$maxay))

	dimnames(descriptivesP)<-list(c("Total","Mean","SD","Skewness","Kurtosis",
						  "Min","Max"), c("X","A"))
	dimnames(descriptivesQ)<-list(c("Total","Mean","SD","Skewness","Kurtosis",
					        "Min","Max"), c("Y","A"))

	kert<-object$kert
	design<-object$design
	SEEYx<-object$SEEYx
	SEEXy<-object$SEEXy

	equatedVal<-data.frame(Score=object$score.x,eqYx=object$eqYx,
	eqXy=object$eqXy,SEEYx=SEEYx,SEEXy=SEEXy)

	bandwidthVal<-data.frame(hx=object$h.x,hy=object$h.y)
	res<-list(call=object$call,equatedVal=equatedVal,bandwidthVal=bandwidthVal,
		    descriptivesP=round(descriptivesP,4),
		    descriptivesQ=round(descriptivesQ,4),design=design,kert=kert)
	}

	class(res)<-"summary.ker.eq"
res
}


print.summary.ker.eq<-function(x,...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\nSummary statistics\n")
		if(x$design=="EG" | x$design=="SG"){
		print(x$descriptives)
		}
		else if(x$design=="CB"){
		print(x$descriptives12)
		cat("\n")
		print(x$descriptives21)
		}
		else if(x$design=="NEAT_CE" | x$design=="NEAT_PSE"){
		print(x$descriptivesP)
		cat("\n")
		print(x$descriptivesQ)
		}
	cat("\nBandwidth parameters used\n")
	print(x$bandwidthVal)
	cat("\nKernel type used\n")
	if(x$kert=="gauss"){
	print("Gaussian")}
	else if(x$kert=="logis"){
	print("Logistic")}
	else if (x$kert=="unif"){
	print("Uniform")}
	cat("\nEquated values and SEE under the",x$design,"design\n")
		if(x$design=="EG" | x$design=="SG"){
		print(x$equatedVal)
		}
		else if(x$design=="CB"){
		print(x$equatedValYx)
		cat("\n")
		print(x$equatedValXy)
		}
		else if(x$design=="NEAT_CE"){
		print(x$equatedVal)
		}
		else if(x$design=="NEAT_PSE"){
		print(x$equatedVal)
		}
	
}


