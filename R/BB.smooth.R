###
### This is a modified version of the functions provided in ER, by IOWA group
### We have translated the C code into R. The original ER softwre is distributed 
### Under the XXX License. 

### The C function Smooth_BB() originally written in C has changed the name to 
### BB.smooth() 

### Future work include to make BB.smooth independent of equate package, specia
### lly when using the function scales(). So far, we just ask the user to load
### the equate library before using BB.smooth().

### Hanson, B. A. (1991) Method of moments estimates for the four
### parameter beta compund binomial model and the calculation of
### classification consistency estimates. (ACT Research Report 91-5). 
### Iowa City, IA: American College Testing

BB.smooth<-function(x,nparm=4,rel) 
  UseMethod("BB.smooth")

#library(moments)
BetaMoments = function(n, k, rmoment){
  #Given raw score mean, s.d., skewness and kurtosis,
  #compute true score mean, s.d., skewness and kurtosis
  
  #Input
  #n = number of items
  #k = Lord k
  #rmoment = raw score mean, s.d., skewness and kurtosis
  
  #Output
  #tmoment = True score mean, s.d., skewness and kurtosis. 
  #If true score variance is negative,
  #then s.d., skew and kurt are set to zero.
  #nctmoment = Non-central true score moments
  cmoment=numeric(4)		                        # central moments 
  ncmoment=numeric(4)	                          # non-central moments 
  fmoment=numeric(4)      		                  # factoral moments
  nctmoment=numeric(4)
  tmoment=numeric(4)
  prod2=prod3=n2=dn=dnm2=kr=0
  
  #compute raw score central moments
  
  cmoment[3] = rmoment[3] * rmoment[2] * rmoment[2] * rmoment[2] #tercer momento centrado
  cmoment[2] = rmoment[2] * rmoment[2]          # segundo momento centrado o varianza
  cmoment[4] = rmoment[4] * cmoment[2] * cmoment[2] #cuarto momento centrado
  
  #compute raw score non-central moments
  
  prod2 = rmoment[1] * rmoment[1];
  ncmoment[2] = cmoment[2] + prod2;
  prod3 = prod2 * rmoment[1];
  ncmoment[3] = cmoment[3] + 3.0 * cmoment[2] * rmoment[1] + prod3;
  ncmoment[4] = cmoment[4] + 4.0 * rmoment[1] * cmoment[3] +
    6.0 * prod2 * cmoment[2] + prod3 * rmoment[1];
  
  # compute raw factoral moments 
  
  fmoment[2] = ncmoment[2] - rmoment[1];
  fmoment[3] = ncmoment[3] - 3.0 * ncmoment[2] + 2.0 * rmoment[1];
  fmoment[4] = ncmoment[4] - 6.0*ncmoment[3] + 11.0*ncmoment[2] - 6.0*rmoment[1];
  
  # compute true score non-central moments
  
  # first moment 
  dn = n;
  n2 = n*(n-1);
  # second moment 
  dnm2 = 1.0;
  kr = k * 2.0;
  ncmoment[1] = rmoment[1] / dn;
  nctmoment[1] = ncmoment[1];
  ncmoment[2] = (fmoment[2]/dnm2 + kr*ncmoment[1]) / (n2+kr);
  nctmoment[2] = ncmoment[2];
  # third moment
  dnm2 = dn-2.0;
  kr = k * 6.0;
  ncmoment[3] = (fmoment[3]/dnm2 + kr*ncmoment[2]) / (n2+kr);
  nctmoment[3] = ncmoment[3];
  # fourth moment 
  dnm2 = dnm2*(dn-3.0)
  kr = k * 12.0;
  ncmoment[4] = (fmoment[4]/dnm2 + kr*ncmoment[3]) / (n2+kr);
  nctmoment[4] = ncmoment[4];
  
  # put true score moments on observed score scale */
  
  dn = dn*dn;
  ncmoment[2] = ncmoment[2]*dn;
  dn = dn*n;
  ncmoment[3] = ncmoment[3]*dn;
  dn = dn*n;
  ncmoment[4] = ncmoment[4]*dn;
  
  # compute true score central moments
  
  prod2 = rmoment[1] * rmoment[1];
  cmoment[2] = ncmoment[2] - prod2;
  prod3 = prod2 * rmoment[1];
  cmoment[3] = ncmoment[3] - 3.0 * ncmoment[2] * rmoment[1] + 2.0*prod3;
  cmoment[4] = ncmoment[4] - 4.0 * rmoment[1] * ncmoment[3] +
    6.0 * prod2 * ncmoment[2] - 3.0 * prod3 * rmoment[1];
  
  # compute true score mean, s.d., skewness and kurtosis 
  
  tmoment[1] = rmoment[1];	
  #The central true score moment may be negative, if so assign
  #the true score s.d., skewness and kurtosis to zero. */
  if (cmoment[2] > 0.0){
    tmoment[2] = sqrt(cmoment[2]);
    tmoment[3] = cmoment[3] / (tmoment[2] * tmoment[2] * tmoment[2]);
    tmoment[4] = cmoment[4] / (cmoment[2] * cmoment[2]);
  }else{
    tmoment[2] = 0.0;
    tmoment[3] = 0.0;
    tmoment[4] = 0.0;
  }
  return(list(t=tmoment,nct=nctmoment))
} 

#calculamos los momentos para las ecuaciones
CalcLordk=function(kr20,nitems,rmoment){
  
  # Calculate value of Lord's k given KR20, number of items,
  #and first two moments of observed raw score distribution
  
  #Input
  #kr20 = KR20 reliability
  #nitems = number of items 
  #rmoment = raw score mean and standard deviation 
  
  #Return k = Lord's k
  n = nitems;
  
  mnm = rmoment[1] * (n - rmoment[1]);
  varr = rmoment[2]*rmoment[2];
  
  #calculate variance of item difficulties
  
  varp = mnm/n^2 - (varr/n) * (1.0 - ((n - 1.0)/n) * kr20)
  # calculate k
  
  k =(n^2*(n-1.0)*varp)/(2*(mnm - varr - n*varp))
  if(k <0){k=0}
  return(k)
} 

#determinamos el valor de k
CalcBetaParaMM=function(nitems,moment){
  
  
  # calculate alpha and beta
  y4 = moment[4];
  y3.2 = moment[3] * moment[3];
  
  r = (6.0 * (y4-y3.2-1.0))/(6.0 + 3.0*y3.2 - 2.0*y4)
  
  rad=(24.0 * (r+1.0))/((r+2.0)*(r+3.0)*y4-3.0*(r + 1.0)*(r - 6.0))
  
  
  if(rad > 1.0){return("cannot compute alpha and beta")}else{rad = sqrt(1.0-rad)}
  
  if(moment[2] < 0.0)  	             #compute values of alpha and beta
  {
    a = (r/2)*(1.0 + rad);
    b = (r/2)*(1.0 - rad);
  }else
  {
    a = (r/2)*(1.0 - rad);
    b = (r/2)*(1.0 + rad);
  }
  
  if (a <= 0.0 || b <= 0.0){return("cannot compute alpha and beta")}
  
  #compute values of lower and upper limit parameters
  
  bma = (a+b) * sqrt(a+b+1)
  bma = bma*moment[2]
  bma = bma/sqrt(a*b)
  
  l = -bma * (a/(a+b))
  l = l+moment[1]		#	                            /* compute lower limit */
  
  u = bma + l
  
  # assign values for output
  para=numeric(4)      
  para[3] = l/nitems;  
  para[4]= u/nitems; 
  para[1]=a; 
  para[2]=b;
  
  if (para[3] < 0.0 || para[4] > 1.0) return("cannot compute u and l ");
  #if (para[3] < 0.0){para[3]=0}
  #if (para[4] > 1.0){para[4]=1}
  return(para)
} 

# Calculamos el valor de los par?metros a,b,l y u
Beta4Smooth=function(x,kr20=0,klord){ 
  xy=c() 
  for(i in 1:length(scales(x))){xy=c(xy,rep(scales(x)[i],x[[i]]))} 
  #Par?metros 
  #sd1=sqrt((1/length(xy))*sum((xy-mean(xy))^2))
  momentos=c(mean(xy),sd(xy),skewness(xy),kurtosis(xy))
  #1.Primero determinamos el k de lord
  #klord=lordk
  #klord=CalcLordk(kr20=kr20,nitems=length(scales(x))-1,rmoment = momentos)
  
  #2.Obtenemos los nuevos momentos
  #tmoment=BetaMoments(n=length(scales(x))-1,k=klord,rmoment = momentos) 
  tmoment=BetaMoments(n=length(scales(x))-1,k=klord,rmoment = momentos) 
  #3. Calculamos los par?metros
  theta=CalcBetaPara(nitems =length(scales(x))-1,moment = tmoment$t,nctmoment = tmoment$nct)  #que pasa si consideramos los mismo parame que antes?
  #4. 
  return(theta)
} 

# Sirve para obtener los 4 parametros de diversas maneras
BB.smooth.default = function(x,nparm=4,rel){
  
  cl <- match.call()
  
  xy=c() 
  for(i in 1:length(scales(x))){xy=c(xy,rep(scales(x)[i],x[[i]]))} 
  #Par?metros 
  #sd1=sqrt((1/length(xy))*sum((xy-mean(xy))^2))
  rmoments=c(mean(xy),sd(xy),skewness(xy),kurtosis(xy))
  
  #  Performs beta binomial smoothing.  There are three types:
  #  2 parameter beta, binomial errors
  #  4 parameter beta, binomial errors
  #  4 parameter bete, compound binomial errors
  
  #Input
  n = sum(x)#number of persons
  nitems = length(scales(x))-1#number of items
  #fd[] = frequency distribution
  #mts[] = raw score moments
  #nparm = number of parameters for beta (2 or 4)
  #rel = reliability (usually KR20); 
  #if rel == 0, binomial errors;
  #otherwise compund binomial errors
  
  #Output
  if(nparm==2){
    lordk=0
    theta=Beta2Smooth(x)
  }else{
    lordk = CalcLordk(kr20=rel,nitems,rmoment=rmoments)
    theta=Beta4Smooth(x,rel,lordk)
  }
  K =length(scales(x))-1;prop.fit=numeric(K+1) 
  for(i in 0:K){
    density_i2=function(eta){
      g=((-theta[3]+eta)^(theta[1]-1)*(theta[4]-eta)^(theta[2]-1))/((theta[4]-theta[3])^(theta[1]+theta[2]-1)*(beta(theta[1],theta[2])))
      Pxi=dbinom(i,K,eta)-lordk*eta*(dbinom(i,K-2,eta)-2*dbinom(i-1,K-2,eta)+dbinom(i-2,K-2,eta))
      return(Pxi*g)
    }
    prop.fit[i+1]=integrate(density_i2,lower=theta[3],upper=theta[4],rel.tol=0.0000001)$value
  }
  prop.fit=prop.fit/sum(prop.fit)
  res<-list(call=cl,freq.est=(prop.fit)*length(xy),prob.est=prop.fit,parameters=theta)
  class(res)<-"BB.smooth"
  return(res)
#  return(list(freq.est=(prop.fit)*length(xy),prob.est=prop.fit,parameters=theta))
} 

Beta2Smooth=function(x){
  
  #Calculate 2 parameter beta binomial distribution.
  
  #Input
  
  #rfreq			         = Raw (# correct) score frequencies
  #  betafit->num_items = Number of items
  
  #  Output
  #  betafit->density	 = fitted raw score distribution (density)
  #  (space allocated before function called)
  
  #  Returns zero if no error, nonzero if error.
  
  #  Function calls other than C or NR utilities:
  #    BetaMoments() 
  #  EstNegHypGeo()
  
  xy=c() 
  for(i in 1:length(scales(x))){xy=c(xy,rep(scales(x)[i],x[[i]]))} 
  #Par?metros 
  #sd1=sqrt((1/length(xy))*sum((xy-mean(xy))^2))
  momentos=c(mean(xy),sd(xy),skewness(xy),kurtosis(xy))      
  
  nsamp = sum(x)
  nitems = length(scales(x))-1
  lordk = 0.0 #;	                 /* binomial error distribution */		
  
  # calculate true score moments */
  
  tmoment=BetaMoments(n=length(scales(x))-1,k=lordk,rmoment = momentos)
  
  # calculate parameters of beta true score distribution  */
  
  para = numeric(4)
  para[3] = 0.0;  para[4] = 1.0
  if(is.numeric(EstNegHypGeo(nitems,tmoment$t,para))==FALSE){
    # if two moments not fit, set alpha=1.0 and fit mean */
    para[3] = 0.0;  para[4] = 1.0;
    para[1] = 1.0;
    para[2] = (nitems/momentos[1]) - para[1];
    if (para[2] < 0.0){para[2]=1}
    param=para
    
  }else{param=EstNegHypGeo(nitems,tmoment$t,para)}
  
  return(param)
} 

#retorna par?metros a y b--> l = 0 y u = 1
CalcBetaPara= function(nitems,moment,nctmoment){
  
  #Estimate parameters of 4-parameter beta binomial model
  
  #Input
  #nitems    = number of items on test
  #moment    = true score moments (mean, s.d., skewness, kurtosis)
  #nctmoment = non-central true score moments
  
  #Output
  #para = method of moments estimates of four parameters of beta
  #distribution (alpha, beta, lower limit, upper limit)
  
  #Returns number of moments fit to obtain estimates
  
  #Function calls other than C or NR utilities:
  # CalcBetaParaMM() 
  #CalcBetaParaLS()
  #EstNegHypGeo()
  
  para=numeric(4)
  #Only try to fit more than one moment if the true score
  #s.d. is greater than zero */
  if (moment[2] > 0.0){
    # try to fit all four moments */
    
    if (is.numeric(CalcBetaParaMM(nitems,moment))==TRUE){return(CalcBetaParaMM(nitems,moment))}
    
    #/* if four moments not fit try to fit three moments 
    #such that the squared difference in the estimated 
    #and observed kurtosis is as small as possible */
    
    if(is.numeric(CalcBetaParaLS(nitems,moment,nctmoment))==TRUE){return(CalcBetaParaLS(nitems,moment,nctmoment))}
    
    #/* if three moments not fit set lower = 0.0 and upper = 1.0 
    #and fit first two moments */
    para[3] = 0.0;  para[4] = 1.0
    if (is.numeric(EstNegHypGeo(nitems,moment,para))==TRUE){return(EstNegHypGeo(nitems,moment,para))}
  }
  
  #if two moments not fit, set alpha=1.0 and fit mean */
  
  para[3] = 0.0;  para[4] = 1.0;
  para[1] = 1.0;
  para[2] = (nitems / moment[1]) - para[1];
  if (para[2] > 0.0){
    return(para)
  }else{
    # if mean not fit return uniform distribution on 0 to 1.0 */
    para[2]=1 
    return(para)
  }
}

EstNegHypGeo=function(nitems,moment,para){
  #Estimate parameters alpha and beta of negative 
  #hypergeometric distribution, given mean, sd, min, and max.
  
  #Input
  #nitems = number of items 
  #moment = true score moments (only first and second used)
  #para = last two elements contain minimum and maximum
  
  #Output
  #para = first two elements contain alpha and beta
  
  
  smean=svar=0	          
  
  dn = nitems;
  maxmin = dn * (para[4] - para[3])  #difference between minimum and maximum
  smean = (moment[1] - dn*para[3])/maxmin # standardized mean and variance
  svar = (moment[2]*moment[2])/(maxmin*maxmin)
  
  #ppara = para
  ppara = smean*smean*(1.0 - smean);
  ppara = ppara/svar;
  ppara= ppara-smean; #alpha?
  para[1] = ppara
  
  #ppara = para+1;
  ppara = smean * (1.0-smean);
  ppara = ppara/svar;
  ppara = ppara-1.0;
  ppara = ppara-para[1];
  para[2]=ppara
  if (para[1] > 0.0 && para[2] > 0.0){return(para)}
}

CalcKurt=function(para){
  
  a = para[1];
  b = para[2];
  k1 = 3.0 * (a+b+1.0);
  k2 = 2.0 * (a+b) * (a+b);
  k3 = a*b;
  k4 = a+b-6.0;
  k5 = a + b + 2.0;
  k6 = a+b+3.0;
  
  kurt = (k1*(k2 + k3*k4))/(k3*k5*k6)
  return(kurt)	
}

FindLower=function(para,nctmoment){
  #Find lower limit of 4-parameter beta-binomial model
  #that for input value of upper limit produces moments
  #that match up to third moment. 
  
  #Input
  #para = parameters of beta binomial dist (para[4] is input)
  #tmoment = non-central true score moments
  
  #Output
  #para = lower limit  (para[3] is output)
  
  #Returns 1 if successful; otherwise returns 0.
  
  
  #register	double	fr1,	                             /* temporary */
  #m1,m2,m3,  /* first through third central moments */
  #upper;	/* upper limit of true score distribution */
  #double	numerator;	       /* numerator of expression for upper limit */
  
  upper = para[4];
  m1 = nctmoment[1];
  m2 = nctmoment[2];
  m3 = nctmoment[3];
  
  # calculate lower limit
  
  #numerator
  fr1 = m1 * m1 * (m2*upper - 2.0 * m3);
  fr1 = fr1+m2*m2 * (m1 - 2.0 * upper);
  fr1 = fr1+m3 * (m2 + m1 * upper);
  numerator = fr1;
  
  #denominator
  fr1 = m1*m1 * (2.0 * m1 * upper - m2);
  fr1 = fr1+m2 * (2.0*m2 - 3.0*m1*upper);
  fr1 = fr1+m3 * (upper - m1);
  
  para[3] = numerator / fr1;
  
  if (para[4] <= para[3] || para[3] < 0.0){ return("Invalid Values")}else{para[3]} # invalid */
} #retorna solo el valor de "l" dado los momentos y dado "u" (duda si es t o nct el imput)

FindUpper=function(para,nctmoment){
  
  # Find upper limit of 4-parameter beta-binomial model
  #that for input value of lower limit produces moments
  #that match up to third moment.  
  
  #Input
  #para = parameters of beta binomial dist. (para[2] is input)
  #tmoment = non-central true score moments
  
  #Ooutput
  #para = upper limit of beta binomial dist. (para[3] is output)  
  
  lower = para[3];
  m1 = nctmoment[1];
  m2 = nctmoment[2];
  m3 = nctmoment[3];
  
  #calculate upper limit */
  
  # numerator */
  fr1 = m1 * m1 * (m2*lower - 2.0 * m3);
  fr1 = fr1+m2*m2 * (m1 - 2.0 * lower);
  fr1 = fr1+m3 * (m2 + m1 * lower);
  numerator = fr1;
  
  # denominator */
  fr1 = m1*m1 * (2.0 * m1 * lower - m2);
  fr1 = fr1+m2 * (2.0*m2 - 3.0*m1*lower);
  fr1 = fr1+m3 * (lower - m1);
  
  para[4] = numerator / fr1;
  
  if (para[4] <= para[3] || para[4] > 1.0){return("Invalid Values")}else{para[4]} # invalid */
} 

#retorna solo el valor de "u" dado los momentos y dado "l" (duda si es t o nct el imput)
Kurtfuncd=function(l,nitems,tmoment,nctmoment){
  # Using input lower limit find upper limit that matches skewness,
  #and using these parameters compute squared difference of
  #predicted and observed kurtosis.
  
  #Input
  #x = lower limit (0 - 1 scale)
  #nitems = number of items 
  #tmoment = true score moments (mean, s.d., skewness, kurtosis)
  #nctmoment = non-central true score moments
  
  #Output
  #kurt = squared difference in observed and predicted kurtosis
  #para = parameters of four-parameter beta binomial
  
  
  #register double tkurt;
  #
  para=numeric(4)
  kurt = -1.0;#	              /* initialize in case of return of 0 */
  para[3] = l;
  if (is.numeric(FindUpper(para,nctmoment))==FALSE){return("Inv?lido FindUpper")}else{para[4]=FindUpper(para,nctmoment)};
  if (is.numeric(EstNegHypGeo(nitems,tmoment,para))==FALSE) return("Ivalido NegHypGeo");
  para=EstNegHypGeo(nitems,tmoment,para)
  tkurt = CalcKurt(para);
  tkurt = tmoment[4] - tkurt;
  tkurt = tkurt^2;
  kurt = tkurt;
  return(list(kurtdif=kurt,parametros=para));	
} 

#retorna diferencia en kurtosis con u ajustado
KurtfuncUpper=function(u, nitems, tmoment,nctmoment){
  #Using input upper limit find lower limit that matches skewness,
  #and using these parameters compute squared difference of
  #predicted and observed kurtosis.
  
  #Input
  #u = upper limit (0 - 1 scale)
  #nitems = number of items 
  #para = beta-binomial parameters (lower limit used for input)
  #tmoment = true score moments (mean, s.d., skewness, kurtosis)
  #nctmoment = non-central true score moments
  
  #Output
  #kurt = squared difference in observed and predicted kurtosis
  
  
  #register double tkurt;
  para=numeric(4)
  kurt = -1.0;	      #        /* initialize in case of return of 0 */
  para[4] = u;
  if (is.numeric(FindLower(para,nctmoment))==FALSE){return("Inv?lido FindLower")}else{para[3]=FindLower(para,nctmoment)};
  if (is.numeric(EstNegHypGeo(nitems,tmoment,para))==FALSE) return("Ivalido NegHypGeo");
  para=EstNegHypGeo(nitems,tmoment,para)
  tkurt = CalcKurt(para);
  tkurt = tmoment[4] - tkurt;
  tkurt = tkurt^2;
  kurt = tkurt;
  return(list(kurtdif=kurt,parametros=para));
} 

#retorna diferencia cuadr?tica entre kurtosis obs y fit
CalcBetaParaLS = function(nitems,tmoment,nctmoment){
  
  #Find value of lower limit of 4-parameter beta compound binomial model
  #that minimizes squared difference between predicted and observed kurtosis, 
  #such that observed and predicted mean, variance and skewness are equal.
  #Used when method of moments fails to produce acceptable values
  #of all four parameters.
  
  #When the solution that minimizes the squared difference in kurtosis
  #is not the solution with the lower limit = 0 nor the solution with
  #the upper limit = 1 it may be missed since an initial solution in this
  #case is found by a rough grid search.
  
  #Input
  #nitems    = number of items 
  #tmoment   = true score moments (mean, s.d., skewness, kurtosis)
  #nctmoment = non-central true score moments
  
  #Output
  #para = method of moments estimates of parameters
  
  #Returns 1 if solution found, otherwise returns 0.
  #Accuracy of solution is +-KACC.
  
  
  #  kurt=0 #,	/* squared difference of fitted and observed kurtosis */
  #  ll,tll;		                            /* holds lower limits */
  #  double tpara[4];	           /* temporary space for beta parameters */
  #  double delta,rkurt,lkurt;
  #  int existr, existl, i;
  
  
  #/* Initialize squared kurtosis difference to be infinity.
  #If no valid solution is found kurt will never be less than
  #this. HUGE_VAL is a macro that should be declared in <math.h>
  # for any compiler following the ANSI C standard. It represents
  #the largest representable floating point number or the floating
  #point representation of infinity. */
  
  kurt = Inf;
  
  # find solution for lower limit of 0 */
  
  existl = Kurtfuncd(0.0,nitems,tmoment,nctmoment);
  if (is.list(existl)==TRUE){kurt = existl$kurtdif;para=existl$parametros}
  
  #/* find solution for upper limit of 1 */
  existr = KurtfuncUpper(1.0,nitems,tmoment,nctmoment);
  if (is.list(existr)==TRUE)
  {
    if(existr$kurtdif < kurt){
      kurt = existr$kurtdif
      para = existr$parametros
    }
  }
  
  #/* attempt to find a solution with a lower squared difference
  #in kurtosis than two solutions found above */
  tll=0  
  delta = .05;
  ll = delta;
  existl = 0;	#	 /* existl will flag whether a better solution is found */
  while (ll < .55){
    exist.p=Kurtfuncd(ll,nitems,tmoment,nctmoment)
  #  if (class(exist.p)=="list") ppp
    if (is(exist.p,"list"))
    {
      if (exist.p$kurtdif < kurt) #/* a better solution is found if kurtl < kurt */
      {
        kurt = exist.p$kurtdif
        para = exist.p$parametros
        tll = ll;
      }
    }				
    ll = ll+delta;
  }
  
  if (kurt == Inf){return("No valid solution has been found")}
  if(tll==0){
    return(para)
  }
  #  /* if no better solution than lower limit = 0 or
  #upper limit = 1 has been found then return */
  
  #if (!existl) return(1);
  
  #/* loop to find solution to accuracy of KACC */
  delta = .01;
  ll = tll;
  while (delta >= 1e-07)
  {
    #/* evaluate function at points to left and right 
    #of current solution "ll" */
    
    existr = Kurtfuncd(ll+delta,nitems,tmoment,nctmoment);
    existl = Kurtfuncd(ll-delta,nitems,tmoment,nctmoment);
    if (is.list(existr)==TRUE) #	               /* step to the right */
    {
      if(existr$kurtdif < kurt){
        ll = ll-delta;
        exist.aux=Kurtfuncd(ll+delta,nitems,tmoment,nctmoment)
        while(((ll+delta) < 1.0) & exist.aux$kurtdif < kurt & is.numeric(exist.aux$kurtdif)==TRUE)
        {
          
          ll = ll+delta;
          kurt = exist.aux$kurtdif
          exist.aux=Kurtfuncd(ll+delta,nitems,tmoment,nctmoment)
        }
      }
    }else{
      
      if (is.list(existl)==TRUE){#	         /* step to the left */
        if (existl$kurtdif < kurt){
          ll = ll-delta
          exist.aux2=Kurtfuncd(ll-delta,nitems,tmoment,nctmoment)
          while(((ll-delta) > 0.0) & exist.aux2$kurtdif < kurt & is.numeric(exist.aux2$kurtdif)==TRUE)
          {
            ll = ll-delta;
            kurt = exist.aux2$kurtdif
            exist.aux2=Kurtfuncd(ll-delta,nitems,tmoment,nctmoment)
          }
        }
      }
    }
    delta = delta/10.0;
  }   #                                                /* end while loop */
  
  # /* calculate final values to return */
  Fit=Kurtfuncd(ll,nitems,tmoment,nctmoment)
  if (is.list(Fit)==TRUE){return(Fit$parametros)}else{return("Invalid parameters")}
}

#### Print function

print.BB.smooth<-function(x,...)
{
  cat("\nCall:\n")
  print(x$call)
    cat("\nEstimated frequencies and score probabilities:\n")
    cat("\n")
    print(data.frame(Est.Freq.=x$freq.est,Est.Score.Prob.=x$prob.est))
    cat("\nParameters:\n")
    cat("\n")
    print(x$parameters)
 
}

