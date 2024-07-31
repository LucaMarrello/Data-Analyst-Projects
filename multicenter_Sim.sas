options nocenter nodate nonotes ls=150 ps=1000 formdlim=' ' 
     formchar= '|----|+|---+=|-/\<>*'; *  mprint ;
%let alpha=0.05;
%let crit = probit(1-&alpha/2);
%LET NULL= .5;
*%let runs= 1000;

%macro ODSOff(); /* Call prior to BY-group processing */
ods graphics off;
ods exclude all;
ods noresults;
%mend;

%macro ODSOn(); /* Call after BY-group processing */
ods graphics on;
ods exclude none;
ods results;
%mend;



%macro bell;
*plays the trumpet call, useful to put at end of batch program to
   know when the batch file has ended;
data _null_;
call sound(659.25,100); 
call sound(659.25,100);
call sound(192.00,70); *first argument is frequency, second is duration;
call sound(523.25,200);
call sound(392.00,70); 
call sound(523.25,100); 
call sound(523.25,100); 
call sound(523.25,100); 
call sound(659.25,200);
call sound(1059.25,200);
run;
%mend;

/***********************************************************************
 Programs for 
 Wicklin, Rick, 2013, Simulating Data with SAS, SAS Institute Inc., Cary NC.

 Appendix B: Generating Multivariate Ordinal Variables

 This program saves functions to the SAS/IML storage library.
 To use these functions in a SAS/IML program, run this program
 and then load the functions, as follows:

 %include "C:\<path>\RandMVOrd.sas";
 proc iml;
 load module=_all_;
 <...call the functions...>

 ***********************************************************************/

proc iml;

    /* P1   P2    P3  */
/* example matrix of PMFs */
/*
P = {0.25  0.50  0.20 ,
     0.75  0.20  0.15 ,
      .    0.30  0.25 ,
      .     .    0.40 };
*/

/* OrdN: number of values for each variable */
start OrdN(P);
   return( countn(P, "col") );
finish;

/* OrdMean: Expected value for each variable is Sigma_i (i*p[i])    */
start OrdMean(P);
   x = T(1:nrow(P));                 /* values of ordinal vars      */
   return( (x#P)[+,] );              /* expected values E(X)        */
finish;

/* OrdVar: variance for each variable */
start OrdVar(P);
   d = ncol(P);   m = OrdMean(P);
   x = T(1:nrow(P));                 /* values of ordinal vars      */
   var = j(1, d, 0);
   do i = 1 to d;
      var[i] = sum( (x - m[i])##2 # P[,i] );    /* defn of variance */
   end;
   return( var );
finish;

/* OrdCDF: Given PMF, compute CDF = cusum(PDF) */
start OrdCDF(P);
   cdf = j(nrow(P), ncol(P));        /* cumulative probabilities    */
   do i = 1 to ncol(P);
      cdf[,i] = cusum(P[,i]);
   end;
   return( choose(P=., ., cdf) );    /* missing vals for short cols */
finish;


/* Function that returns ordered pairs on a uniform grid of points.
   Return value is an (Nx*Ny x 2) matrix */
start Expand2DGrid( _x, _y );
   x  = colvec(_x); y  = colvec(_y);
   Nx = nrow(x);    Ny = nrow(y);
   x = repeat(x, Ny);
   y = shape( repeat(y, 1, Nx), 0, 1 );
   return ( x || y );
finish;


/* OrdQuant: Compute normal quantiles for CDF(P) */
start OrdQuant(P);
   N = OrdN(P);
   CDF = OrdCDF(P);
   /* QUANTILE function in SAS/IML 12.1 does not accept 1 as parameter */
   /* Replace 1 with missing value to prevent error */
   CDF = choose(CDF > 1 - 2e-6, ., CDF);
   quant = quantile( "Normal", cdf );
   do j = 1 to ncol(P);      /* set upper quantile to .I = infinity */
      quant[N[j],j] = .I;    /* .I has special meaning to BIN func  */
   end;                      
   return( quant );
finish;

/* OrdFindRoot: Use bisection to find the MV normal correlation that 
   produces a specified MV ordinal correlation. */
start OrdFindRoot(P1, P2,  target);
   N1 = countn(P1);   N2 = countn(P2);
   q1 = OrdQuant(P1); q2 = OrdQuant(P2);
   v1 = q1[1:N1-1];   v2 = q2[1:N2-1];
   g = Expand2DGrid(v1, v2);
   /* find value of rho so that sum(probbnrm(g[,1], g[,2], rho))=target */
   /* Bisection: find root on bracketing interval [a,b] */
   a = -1; b = 1;                 /* look for correlation in [-1,1] */
   dx = 1e-8; dy = 1e-5;
   do i = 1 to 100;               /* iterate until convergence      */
      c = (a+b)/2;
      Fc = sum( probbnrm(g[,1], g[,2], c) ) - target;
      if (abs(Fc) < dy) | (b-a)/2 < dx then 
         return(c);
      Fa = sum( probbnrm(g[,1], g[,2], a) ) - target;
      if Fa#Fc > 0 then a = c;
      else b = c;
   end;
   return (.);                    /* no convergence                 */
finish;

/* alternative root-finding algorithm that uses FROOT (SAS 12.1) */
/*
start MMRoot(x) global(_grid, _target);
   return( sum( probbnrm(_grid[,1], _grid[,2], x) ) - _target );
finish;

start AltOrdFindRoot(P1, P2,  target) global(_grid, _target);
   N1 = countn(P1);   N2 = countn(P2);
   q1 = OrdQuant(P1); q2 = OrdQuant(P2);
   v1 = q1[1:N1-1];   v2 = q2[1:N2-1];
   _grid = Expand2DGrid(v1, v2);
   _target = target;
   return( froot("MMRoot", {-1 1}) );
finish;
*/

/* OrdMVCorr: Compute a MVN correlation matrix from the PMF and 
   the target correlation matrix for the ordinal variables. */
start OrdMVCorr(P, Corr);
   d = ncol(P);
   N = OrdN(P);
   mean = OrdMean(P);
   var  = OrdVar(P);
   cdf  = OrdCDF(P);
   R = I(d);
   do i = 1 to d-1;
      sumCDFi = sum(cdf[1:N[i]-1, i]); 
      do j = i+1 to d;
         sumCDFj = sum(cdf[1:N[j]-1, j]); 
         hStar = Corr[i,j] * sqrt(var[i]*var[j]) + mean[i]*mean[j] 
                 - N[i]*N[j] + N[i]*sumCDFj + N[j]*sumCDFi;
         R[i,j] = OrdFindRoot(P[,i], P[,j], hStar);
         R[j,i] = R[i,j];
      end;
   end;
   return(R);
finish;

/* RandMVOrdinal: 
   N     Number of desired observations from MV ordinal distribution, 
   P     Matrix of PMF for ordinal vars. The j_th col is the j_th PMF.
         Use missing vals if some vars have fewer values than others.
   Corr  Desired correlation matrix for ordinal variables. Not every
         matrix is a valid as the correlation of ordinal variables. */
start RandMVOrdinal(N, P, Corr);
   d = ncol(P);
   C = OrdMVCorr(P, Corr);     /* 1. compute correlation matrix, C  */
   mu = j(1, d, 0);
   X = RandNormal(N, mu, C);   /* 2. simulate X ~ MVN(0,C)          */
   Ni = OrdN(P);
   quant = OrdQuant(P);        /* compute normal quantiles for PMFs */
   do j = 1 to d;              /* 3. convert to ordinal             */
      X[,j] = bin(X[, j], .M // quant[1:Ni[j],j]);
   end;
   return(X);  * other prog used X` ;
finish;
 


store module=_all_;


quit;






ods graphics off;

*%inc "C:\Users\gy.zou\Box Sync\Wilkin\RandMVord.sas";


%macro oneRun(run=, effect=, N=, cmin = 4, cmax  = 10,
        miss =, rho1=, rho0=, r12 = );
*cmin, min for center size, cmax, max for center size;

ods listing;
 

proc iml;
 
*d0 = date();
*t0 = time();

load module=_all_;

call randseed(&run*123987); ** don not forget seeds;


rho0= &rho0; *0.1;
rho1= &rho1; *0.05;

     if &r12=1 then rho12=rho1;
else if &r12=2 then rho12=rho1/2;
else if &r12=3 then rho12=0;

/*
upper5 = {.5 1 1 1 1,
          0 .5 1 1 1,
          0 0 .5 1 1, 
          0 0 0 .5 1,
          0 0 0 0 .5};


p0 = pdf('Binom', 0:4, .1, 4);  * binomial(n, p); 

     if &effect=1 then p1=p0;
else if &effect=2 then p1 = pdf('Binom', 0:4, .14, 4); *for param=0.56;
else if &effect=3 then p1 = pdf('Binom', 0:4, .2, 4); * for param=0.639;
else if &effect=4 then p1 = pdf('Binom', 0:4, .26, 4); * param=0.71;
 


 *munzel pharmaceutical Stat;
 p0= {.098 .455 .196 .188 .063};

     if &effect =1 then p1 = p0;
else if &effect =2 then p1 = {.2 .2 .2 .2 .2};
else if &effect =3 then p1 = {.063 .188 .196 .455 .098};


* Boos;

p0 ={.01 .30 .48 .20 .01};
  
     if &effect =1 then p1 = p0; *{.01 .30 .48 .20 .01 } ;
else if &effect =2 then p1 =      {.02 .15 .34 .40 .09}; 
else if &effect =3 then p1 =      {.02 .10 .25 .50 .13};
 */






p0  =  pdf('Binom', 0:6, .5, 6);  * binomial(n, p); 

      if &effect=1 then p1=p0;
      if &effect=2 then p1 =  pdf('Binom', 0:6, .55, 6) ;   
      if &effect=3 then p1 =  pdf('Binom', 0:6, .61, 6) ;  
      if &effect=4 then p1 =  pdf('Binom', 0:6, .66, 6) ;  

upper  =
    {.5 1 1 1 1 1 1,
     0 .5 1 1 1 1 1,
     0 0 .5 1 1 1 1, 
     0 0 0 .5 1 1 1,
     0 0 0 0 .5 1 1,
     0 0 0 0 0 .5 1,
     0 0 0 0 0 0 .5}; 


param =  p0 * upper * p1`; 

*categories = ncol(p1);
*print  p1, p0, param 'cat=' categories;
 

 
call symputx("para", param);  



totCenter = &N;

*free simOne;

free dataI;

do center =1 to totCenter;

 halfsize = ( &cmin + floor((1+ &cmax-&cmin)*rand("uniform")));
* size for each trt within center generate from U(4,10);
 if halfsize = &cmax then ind1=J(1, 2*halfsize,1);
 else if halfsize < &cmax then ind1 = J(1, halfsize, 1)||J(1, &cmax-halfsize, 0) ||
                                      J(1, halfsize, 1)||J(1, &cmax-halfsize, 0) ;
   
 ttt=ind1;

 if &miss=1 then do;
    tt = randfun(2*&cmax, "Bernoulli", 0.8) ;
     ttt = ind1 # t(tt) ;
 end;


 dataI = dataI//ttt;


end;

* print 'data idx=' dataI;

Idx = colvec(dataI);


*print Idx;


free center;
do i =1 to totCenter;
   center = center// (repeat(i, 2*&cmax, 1)|| (repeat(1, &cmax, 1) //repeat(0, &cmax, 1))) ;
end;
 

csize1 = &cmax; 
csize0= csize1;

csize=csize0+csize1; 

corr1 =(1-rho1)* I(csize1)+ J(csize1, csize1, rho1);
corr0 =(1-rho0)* I(csize0)+ J(csize0, csize0, rho0);

matR = (corr1 || (rho12*J(csize1, csize0, 1)))//
       ( (rho12*J(csize0, csize1, 1))|| corr0 ); 
 
p11 = repeat(p1`, 1, csize1);
p00 = repeat(p0`, 1, csize0);

PP = p11 || p00;


X = RandMVOrdinal(totCenter, pp, matR); *each row is a center;

*print X;

XX = center || colvec(X)|| Idx; 

*print XX;

XXX = XX[,1:3][ loc(XX[,4]=1), ];
 
 
*simOne = simOne//YY;

*end;  *tot num cluser;
 
*print SimOne; *Y tt YY; 

create simData from XXX [colname= { "center" "trt" "outcome"}];    
     append from XXX;
close simData;


 

free XXX;
  
quit;
 


*proc print data=simData;
run;



*** Win Fractions are created within cneters,
  Boos & Brownie (1992);


title '==== Conditional approach ===';
ods listing close;
proc sort data= simData; 
   by center descending trt; run;

proc freq data= simData; by center;
   table trt/out=freqcnt (drop=PERCENT);
run; 
   
data newn; 
  set freqcnt;  by center;
   trt = 1 - trt; 
run;  


proc rank data = simData out=allrank; 
     by center;
      var  outcome;
	  ranks arank;
run;
 
proc rank data= SimData out=grank;  
    by center descending trt;
     var  outcome;
	ranks grank;
run;


data alln; 
  merge allrank grank newn; by center descending trt;
       WinF =   (arank-grank)/count ;
  run;
 



ods listing close;

title2 Stratum-specific estimates;
 
ods output ttests = t (keep=center Method df where =(Method='Pooled')) 
           ttests = t2 (keep=center Method df where =(Method='Satterthwaite'))
       Statistics = test (keep=center Method Mean StdErr   where=(Method='Satterthwaite'));

proc ttest data=alln order=data;
    class  trt;
    by center;
    var WinF;
run;
 

data estt;
   merge t2 test t(rename=(df=pooledDF));  by center;
   pooledDF=pooledDF+2;
   point = mean/2 + .5;
run;


title2 Weighted by center size;

ods output Statistics=pairedW (keep=N mean Stderr);
proc ttest data=estt h0=.5; 
    var point;
	weight pooleddf;
run;
 

data pairedW ;
   set pairedW ;

   point=mean;
   df=N-1;
  lgtPoint = log(point/(1-point));
	
	crit = tinv(1- &alpha/2, Df);
  * crit=1.96;
   lower = point - &crit*StdErr;
   upper = point + &crit*StdErr;
 
   lft0 = 100*(upper < &para);
   rgt0 = 100*(lower > &para);
   w0  =100* (upper - lower);
	cov0 = 100 - lft0 - rgt0;

	rej0 = 100*((lower > .5) + ( upper < .5)); 
 
    se1 = StdErr/(point*(1-point));
   low = lgtPoint - crit * se1;
    upp = lgtPoint + crit * se1;
    lower1 = logistic(low);
	 upper1 = logistic(upp); 

	lft1  = 100*(upper1 < &para);
   rgt1  = 100*(lower1 > &para);
   w1  = 100*(upper1 - lower1);
	cov1  = 100 - lft1  - rgt1 ;
 
	rej1  = 100*( (lower1 > .5) + (upper1 < .5)); 
 
	bias = 100*(&para-point)/&para; 

    center =&N;

    mis=&miss;

   
   rho12=&r12;
   rho1=&rho1; 
   rho0=&rho0;

param=&para;

keep rho1 rho0 rho12 mis center param point bias 
               cov0  w0 cov1 lft1 rgt1 w1   rej0 rej1 ;
run;
  

 
proc print data= pairedW;
 title PairedW;
run;



proc append base = pairedWT data = pairedW;

run;
 


ods listing close; 

title2 No ineraction;

proc mixed data=alln; 
   class  center trt/ref=first;
   model WinF = trt center;
   lsmeans trt/diff cl;
  ods output Diffs=  est_noint (keep =  Estimate StdErr DF);

run; 

data est_noint;
   set est_noint ;
   point = (Estimate+1)/2;
   lgtPoint = log(point/(1-point));
	
	crit = tinv(1- &alpha/2, Df);
  * crit=1.96;
   lower = point -  crit*StdErr;
   upper = point +  crit*StdErr;

    param = &para;

   lft0 = 100*(upper < param);
   rgt0 = 100*(lower > param);
   w0  = 100*(upper - lower);
	cov0 = 100 - lft0 - rgt0;

	rej0 = 100*((lower > .5) + ( upper < .5));


    se1 = StdErr/(point*(1-point));
   low = lgtPoint - &crit * se1;
   upp = lgtPoint + &crit * se1;
   lower1 = logistic(low);
	upper1 = logistic(upp);
  
	lft1  = 100*(upper1 < param);
   rgt1  = 100*(lower1 > param);
   w1  = 100*(upper1 - lower1);
	cov1  = 100 - lft1  - rgt1 ;
 
	rej1  = 100*( (lower1 > .5) + (upper1 < .5));
 
 
	bias = 100*( param-point)/param; 

    center=&N;
    mis=&miss;

   
   rho12=&r12;
   rho1=&rho1; 
   rho0=&rho0;

keep rho1 rho0 rho12 mis center param point bias 
               cov0  w0 cov1 lft1 rgt1 w1   rej0 rej1 ;
run;
  


proc print data= est_noint;
 title Fixed No Interaction;
run;


proc append base = FixNOINT data = est_noint;

run;

   
*ods listing close;
 * See Feaster et al 2011, p.387, Fixed effects for site 
and site-by-trt interaction; 
title2 Interaction;


proc mixed data=alln; 
   class  center trt/ref=first;
   model WinF = trt center  trt*center;
   lsmeans trt/diff cl;
  ods output Diffs=  est_int (keep =  Estimate StdErr DF);

run;
 
data est_int;
   set est_int ;

   point = (Estimate+1)/2;
   lgtPoint = log(point/(1-point));
	
	crit =  tinv(1- &alpha/2, Df);
  * crit=1.96;
   lower = point - crit*StdErr;
   upper = point + crit*StdErr;

    se1 = StdErr/(point*(1-point));
   low = lgtPoint - &crit * se1;
   upp = lgtPoint + &crit * se1;
   lower1 = logistic(low);
	upper1 = logistic(upp);
   

   param = &para;

   lft0 = 100*(upper < param);
   rgt0 = 100*(lower > param);
   w0  = 100*(upper - lower);
	cov0 = 100 - lft0 - rgt0;

	rej0 = 100*((lower > .5) + ( upper < .5));

	lft1  = 100*(upper1 < param);
   rgt1  = 100*(lower1 > param);
   w1  = 100*(upper1 - lower1);
	cov1  = 100 - lft1  - rgt1 ;
 
	rej1  = 100*( (lower1 > .5) + (upper1 < .5));
 
 
	bias = 100*( param-point)/param; 

    center=&N;
    mis=&miss;

   
   rho12=&r12;
   rho1=&rho1; 
   rho0=&rho0; 

keep   rho1 rho0 rho12 mis center param point bias 
               cov0  w0 cov1 lft1 rgt1 w1   rej0 rej1  ;
run;
  

*ods listing;
proc print data= est_int;
 title Fixed with Interaction;
run;
 

proc append base = FixINT data = est_int;

run;


proc delete data=pairedW est_noint est_int;
run;



***Marginal approach;
*** Win Fractions are created over all subjects in the study
(Brunner et al (1995);

ods listing close;

proc sort data = simData; 
    by descending trt; run;

proc freq data= SimData;  
   table trt/out=freqcnt (drop=PERCENT);
run; 
   

data new; 
  set freqcnt;   
   trt = 1 - trt; 
run;  

*proc print data=new;
run;


proc rank data = SimData out=allrank; 
      var  outcome;
	  ranks arank;
run;

proc rank data= SimData out=grank;  
     by  descending trt;
     var  outcome;
	ranks grank;
run;


data WinF; 
    merge allrank grank new; by  descending trt;
       winF =   (arank-grank)/count ;
  run;
 
*proc delete data=simData; run;

  

***======= ignore clustering =; 
proc mixed data= WinF noitprint noclprint; 
    class  trt/ref=first;
	 model WinF = trt/ddfm= sat; 
	* repeated /group=trt;
    lsmeans trt/diff cl; 
    ods output Diffs=      est0 (keep = Estimate StdErr DF); 
run;
  
data est0;
 
  set est0;  

 	point = (Estimate+1)/2;
	 

   test = (point - .5)/StdErr;
   pvalue = 2*(1-(probt(abs(test), DF)));

  

	crit =  tinv(1- &alpha/2, Df);

   lower = point - &crit*StdErr;
   upper = point + &crit*StdErr;



   param = &para;

   lft0 = 100*(upper < param);
   rgt0 = 100*(lower > param);
   w0  = 100*(upper - lower);
	cov0 = 100 - lft0 - rgt0;

	rej0 = 100*((lower > .5) + ( upper < .5));



 lgtPoint = log(point/(1-point));
    se1 = StdErr/(point*(1-point));
   low = lgtPoint - &crit * se1;
   upp = lgtPoint + &crit * se1;
   lower1 = logistic(low);
	upper1 = logistic(upp);
  
 
	lft1  = 100*(upper1 < param);
   rgt1  = 100*(lower1 > param);
   w1  = 100*(upper1 - lower1);
	cov1  = 100 - lft1  - rgt1 ;
 
	rej1  = 100*( (lower1 > .5) + (upper1 < .5));

   
 
	bias = 100*( param-point)/param; 

    center=&N;
    mis=&miss;

   rho12=&r12;
   rho1=&rho1; 
   rho0=&rho0;

keep rho1 rho0 rho12 mis center param point bias  
           cov0  w0 cov1 lft1 rgt1 w1   rej0 rej1  ;

 
run;

 
proc append base = Mixed0 data=est0;
run;




*ods listing;

proc mixed data= WinF noitprint noclprint; 
    class center trt/ref=first;
	 model WinF = trt/ddfm=betwithin;
     random intercept/subject=center;
     lsmeans trt/diff cl;
    * parms /nobound; 
    ods output Diffs=      est1 (keep = Estimate StdErr DF)
               CovParms=varcomp(keep = CovParm Estimate)  ; 
run;
 




data est1;
  set est1;

 	point = (Estimate+1)/2;
	totalV = (center + trt + resi);
   icc1 = (center + trt)/totalV;
	icc12 = (center)/totalV;


   test = (point - .5)/StdErr;
   pvalue = 2*(1-(probt(abs(test), DF)));

 

   param = &para; 
  crit =  tinv(1- &alpha/2, Df);
   lower = point -  crit*StdErr;
   upper = point +  crit*StdErr;

  lft0 = 100*(upper < param);
   rgt0 = 100*(lower > param);
    w0  = 100*(upper - lower);
	cov0 = 100 - lft0 - rgt0;
	rej0 = 100*((lower > .5) + ( upper < .5));


   lgtPoint = log(point/(1-point));
    se1 = StdErr/(point*(1-point));
   low = lgtPoint - &crit * se1;
   upp = lgtPoint + &crit * se1;
   lower1 = logistic(low);
	upper1 = logistic(upp);
   


	lft1  = 100*(upper1 < param);
   rgt1  = 100*(lower1 > param);
   w1  =100* (upper1 - lower1);
	cov1  = 100 - lft1  - rgt1 ;
 
	rej1  = 100*( (lower1 > .5) + (upper1 < .5));

    
	bias = 100*( param-point)/param; 

    center=&N;

    mis=&miss;

   
   rho12=&r12;
   rho1=&rho1; 
   rho0=&rho0;

keep rho1 rho0 rho12 mis center param point bias   cov0  w0 cov1 
      lft1 rgt1 w1   rej0 rej1  ;

run;

 
proc append base = Mixed1 data=est1;
run;

 

proc mixed data=WinF noitprint noclprint; 
    class center trt/ref=first;
	 model WinF = trt/ddfm=betwithin;

	random intercept trt/subject = center;
          **do not need nest bc center goes from 1 to totcenter;
   
    lsmeans trt/diff cl;
    *parms /nobound; 
    ods output Diffs =      est2 (keep = Estimate StdErr DF)
               CovParms=varcomp (keep = CovParm Estimate)  ; 
run;
 



data center trt resi;
    set varcomp;
	 if covParm = 'Intercept' then output center;
	 if covParm = 'trt' then output trt;
     if covParm = 'Residual' then output resi;
run;
  

data est2;
  merge est2 
            center(rename=(Estimate=center))
            trt(rename=(Estimate=trt))
            resi(rename=(Estimate=resi));
 
 	point = (Estimate+1)/2;
	totalV = (center + trt + resi);

   icc1 = (center + trt)/totalV;
	icc12 = (center)/totalV;

   test = (point - .5)/StdErr;
   pvalue = 2*(1-(probt(abs(test), DF)));

   param=&para;

    crit =  tinv(1- &alpha/2, Df);

   lower = point -   crit*StdErr;
   upper = point +   crit*StdErr;


   lft0 = 100*(upper < param);
   rgt0 = 100*(lower > param);
   w0  = (upper - lower);
	cov0 = 100 - lft0 - rgt0;

	rej0 = 100*((lower > .5) + ( upper < .5));

 lgtPoint = log(point/(1-point));
   se1 = StdErr/(point*(1-point));
   low = lgtPoint - &crit * se1;
   upp = lgtPoint + &crit * se1;
   lower1 = logistic(low);
   upper1 = logistic(upp); 

	lft1  = 100*(upper1 < param);
   rgt1  = 100*(lower1 > param);
    w1  = 100*(upper1 - lower1);
	cov1  = 100 - lft1  - rgt1 ;
 
	rej1  = 100*( (lower1 > .5) + (upper1 < .5)); 

	bias = 100*( param-point)/param; 

    center=&N;
    mis=&miss;

   
   rho12=&r12;
   rho1=&rho1; 
   rho0=&rho0;

keep rho1 rho0 rho12 mis center param point bias 
              icc1 icc12 cov0  w0 cov1 lft1 rgt1 w1  rej0 rej1  ;
run;
   

proc append base = Mixed2 data = est2;
run;


proc delete data= winF allrank grank; run;

 
*ods listing;
 

%mend oneRun;

  




%let runs=1000;

%macro ALLSenarios (out=);

  %do missI = 0 %to 1;   ***0 =1:1 randomization, 1= unequal trt assignment;

  %do r12= 1 %to 2;  * 1: rho12=rho, 2: rho12=.5 rho;

  %do NC = 1 %to 3;  

   %if &NC=1 %then %do; %let NCenter=10; %end;
   %if &NC=2 %then %do; %let NCenter=20; %end;
   %if &NC=3 %then %do; %let NCenter=50; %end; 


 %do r = 1 %to 2;

     %if &r =1 %then %do; %let rho1=0.05; %end;
     %if &r =2 %then %do; %let rho1=0.10; %end; 

 
  %do es = 1 %to 4; *1=0.5, 2=.56, 3=.64, 4=.71;

    %do run =1 %to &runs;

     %oneRun (run = &run, effect=&es, N = &NCenter, miss = &missI,
              rho1=&rho1, rho0=&rho1, r12 =&r12);
    %end; 


    %end;
   %end;
   %end;
  %end;
 %end;

ods listing close;

title Conditional;

proc sort data= pairedWT; by rho1 rho12 mis center param; run;

proc means n mean data= pairedWT   noprint;
   by rho1 rho12 mis center param;
   var point bias  cov0  w0 cov1 lft1 rgt1  w1  rej0 rej1  ;
	output out =  resultPairT (drop=_freq_ _type_) mean=/autoname; 
run;
 data resultPairT;
  mthd=1;
  set resultPairT;


proc sort data= FixNOINT; by rho1 rho12 mis center param; run;
proc means n mean data= FiXNOINT   noprint;
   by rho1 rho12 mis center param;
   var point bias  cov0  w0 cov1 lft1 rgt1 w1   rej0 rej1 ;
	output out =  resultFixNT (drop=_freq_ _type_) mean=/autoname; 
run;
 

data resultFixNT;
     mthd = 2;
  set resultFixNT;


proc sort data= FixINT; by rho1 rho12 mis center param; run;

proc means n mean data= FiXINT   noprint;
   by rho1 rho12 mis center param;
   var point bias  cov0  w0 cov1 lft1 rgt1 w1   rej0 rej1 ;
	output out =  resultFixINT (drop=_freq_ _type_) mean=/autoname; 
run;

data resultFixINT;
   mthd = 3;
  set resultFixINT;
 run;

proc print data=resultFIXINT noobs;
  title2 'Fixed  interaction';
  format _ALL_ 6.4;
run;




title Marginal;

proc sort data=Mixed0; by rho1 rho12 mis center param; run;

proc means n mean data=Mixed0  noprint;
   by rho1 rho12 mis center param;
   var point bias   cov0  w0 cov1 lft1 rgt1 w1  rej0 rej1 ;
	output out =  resultMix0(drop=_freq_ _type_) mean=/autoname; 
run;

data  resultMix0;
   mthd=4;
  set resultMix0;  
  

proc print data=resultMix0 noobs;
  title2 'ignore clustering';
  format _ALL_ 6.4;
run;
 
proc sort data=Mixed1; by rho1 rho12 mis center param; run;

proc means n mean data=Mixed1   noprint;
   by rho1 rho12 mis center param;
   var point bias cov0  w0 cov1 lft1 rgt1  w1   rej0 rej1 ;
	output out =  resultMix1(drop=_freq_ _type_) mean=/autoname; 
run;


data  resultMix1; 
   mthd=5; 
   set resultMix1;
  
proc print data=resultMix1 noobs;
  title2 'no interaction';
  format _ALL_ 6.4;
run;
 

proc sort data=Mixed2; by rho1 rho12 mis center param; run;

proc means n mean data=Mixed2   noprint;
   by rho1 rho12 mis center param;
   var point bias  cov0  w0 cov1 lft1 rgt1 w1   rej0 rej1  ;
	output out =  resultMix2(drop=_freq_ _type_) mean=/autoname; 
run;

data  resultMix2; 
   mthd=6;
   set resultMix2;  


proc print data=resultMix2 noobs;
  title2 'interaction';
  format _ALL_ 6.4;
run;

data &out;
     set resultPairT  resultFixNT  resultFixINT resultMix0 resultMix1 resultMix2;
run;

proc delete data = pairedWT FixNOINT FixINT 
                    Mixed0 Mixed1 Mixed2;
run;



%mend ALLSenarios;




 libname ll 'C:\Users\mrmar\Documents\Thesis';


%ALLSenarios(out=ll.allresults);

    
 
ods listing;
title all results; 
proc sort data=ll.allresults; by mis rho12 center  mthd  param; run;


proc print data=ll.allresults noobs round;
 *format  _ALL_ 5.3;
run; 

data table;
  set ll.allresults;
 
paramm=  put(mthd, 3.0)||' & '|| put(param, 6.3)|| ' & '; 
    * || ' & ('||put(rho1, 4.2)||','||put(rho122, 4.2)||' )';
point = ' & ' ||put(bias_mean, 5.2);
cover = ' & ' ||put(cov1_mean, 5.1) ||' & (' ||put(lft1_mean, 4.2)||', ' ||put(rgt1_mean, 4.2)||
        ') & '|| put(w1_mean, 5.2)|| ' & ' || put(rej1_mean, 5.2)||' & ';
keep rho1 rho12 paramm point cover;
run;
;

proc print data=table;
run;


data rho1;
  set table;
  if rho1=0.05;
  keep   paramm point cover;
run;


data rho2;
  set table;
  if rho1=0.10;
  point1 = point;
  cover1 = cover||' \\';
  keep   paramm point1 cover1;
run;
data all;
  merge rho1 rho2;  
run;

proc print data=all noobs;
  *var paramm point cover;
run;


proc delete data=ll.allresults table;




%bell;
