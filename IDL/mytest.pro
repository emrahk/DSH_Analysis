FUNCTION MYFUNCT, xx, inpar
    f = gaussian(xx, inpar)
RETURN, f
END


;FUNCTION SBP_MODEL, angle, D_src, 
;    
;RETURN, outpar
;END


;GOTO, skip

restore,'DATA/resparex_reid_use1.sav' ; cloud parameters
;restore,'DATA/parcloud_reidfit.sav'
restore,'DATA/parcloud_fit.sav'

wrap_simprof_cont_MC80_mde2reid,resparex_reid,  outsbp=sbpout_fit, syserr=0.04;,/ps,fname='sbpfit_fit.eps'
;wrap_simprof_cont_MC80_mde2reid,resparex_reid, cloudpar=parcloud_fit, outsbp=sbpout_fit, syserr=0.04;,/ps,fname='sbpfit_fit.eps'

GOTO, skip

;;;;; to find redchi with a different data range at the end of the wrapper:
 yyy = where(rad_im GT 20 AND rad_im LE 600)     
 mymatch = interpol(sbptot, findgen(10000L),rad_im(yyy))  
 myerr=NPROF_IM2C67FLATB[1,yyy]+NPROF_TOTC67FLATB[0,yyy]*errsys
 mychi = (NPROF_IM2C67FLATB[0,yyy]-mymatch)/myerr
 myredchi = TOTAL(mychi^2.)/(N_ELEMENTS(yyy)-3.)
 PRINT, myredchi

;;;;;;;;;

;skip:
xx = FINDGEN(100)
f = gaussian(xx, [100,50,20])
PLOT, xx, f

yy = f + RANDOMN(seed, 200) * SQRT(f)
yy_err = SQRT(yy)

OPLOT, xx, yy, psym = 3
OPLOTERR, yy, yy_err

expr = 'GAUSSIAN(X, P(0:2))'		   ; fitting function
   p0 = [80.D, 40.D, 25.]  		   ; Initial guess
   nparam = 3
   myparinfo = REPLICATE({fixed:0, limited:[1,0], limits:[0.0,0.0]}, nparam)
   
; Gaussian parameters:
;   parms[0] = maximum value (factor) of Gaussian,
;   parms[1] = mean value (center) of Gaussian,
;   parms[2] = standard deviation (sigma) of Gaussian.

;res = DIALOG_MENU(['Max = '+STRTRIM(p0[0],2), 'Mean', 'Sigma'], /index, title = 'Which parameter to fix?')
  
   
  ; p = mpfitexpr(expr, xx, yy, yy_err, p0, dof=mydof, bestnorm=mychi, parinfo=myparinfo); Fit the expression
  ; PRINT, p



p = MPFITFUN('MYFUNCT', xx, yy, yy_err, p0)
PRINT, p
oplot, xx, MYFUNCT(xx, p), color = 200	   ; Plot model

GOTO, out

;expr = 'X*(1+X)+P(0)'
;xx = FINDGEN(100)
;f = xx*(1+xx)+3
;yy = f + RANDOMN(seed, 100) * SQRT(f)
;yy_err = SQRT(yy)
;p = MPFITEXPR(expr, xx, f, y_err)



   plot, xx, yy, psym = 3				   ; Plot data
   OPLOTERR, yy, yy_err
   oplot, xx, mpevalexpr(expr, xx, p)		   ; Plot model


;expr = 'f_lor(x,P[0:2])'

;r = MPFITEXPR(expr,ff,pf,pfe,s,dof=dof,perror=perror,/quiet,bestnorm=chi,parinfo=parinfo)


;********* trial fit **********

cloud_dist = parcloud.x * parcloud.distance  ;saving actual estimated distance to MCs

     ;xn=x-j*0.01/D
     ;FOR i=0, n_elements(Foft)-1 DO BEGIN
     ;   dt=tobs-i
 ;   ;    theta=100d*sqrt(dt*(1-xn)*11./(1.4*(xn*D)))
     ;   theta=sqrt(2*c*dt*day*(1-xn)/(xn*D*kpc))*180D*3600D/!PI ;asec 
     ;   theta_sc=theta/(1-xn)
     ;   Omega=((theta_sc/1000.)^(-alpha))*meanE^(-beta)
     ;   intensity=Foft[i]*Omega/(1-xn)^2.
     ;   prof[j-1,0,i]=theta
     ;   prof[j-1,1,i]=intensity	


ang = FINDGEN(600)

c=2.998D10 ;cm/s
kpc=3.086D21 ;cm
day=86400D ;s
   alpha=parcloud.alpha[0]
   beta=parcloud.beta[0]
   meanE=parcloud.meanE[0]

D = 11.5     	    ;to be evaluated
xn = parcloud.x[0]  ;to be calculated based on D
tobs = chand-t0
dt = tobs - FINDGEN(N_ELEMENTS(Foft))
theta=sqrt(2*c*dt*day*(1-xn)/(xn*D*kpc))*180D*3600D/!PI ;asec

intensity = Foft*$
    ((sqrt(2*c*dt*day*(1-xn)/(xn*D*kpc))*180D*3600D/!PI/(1-xn)/1000.)^(-alpha))*$
    meanE^(-beta)/(1-xn)^2.

xint = INTERPOL(intensity, theta, ang)
PLOT, theta, intensity, xr = [0,600.], /NODATA
OPLOT, ang, xint, color = 250

   p0 = [80.D, 40.D, 25.]  		   ; Initial guess
   nparam = 2
   myparinfo = REPLICATE({fixed:0, limited:[1,0], limits:[0.0,0.0]}, nparam)

data_x = rad_im
data_y = TRANSPOSE(NPROF_IM2C67FLATB[0,*])
data_err = NPROF_IM2C67FLATB[1,*]+NPROF_IM2C67FLATB[0,*]*errsys

out:

skip:
END






