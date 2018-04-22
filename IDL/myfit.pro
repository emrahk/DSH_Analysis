;FUNCTION SBP_MODEL, xx, inpar
;    f = gaussian(xx, inpar)
;RETURN, f
;END

;FUNCTION SBP_MODEL, theta_arr, xn_arr, weight, time_hires, Foft_hires, tobs, p ;p0 = [D0, norm0]

FUNCTION SBP_MODEL, theta_arr, p ;p0 = [D0, norm0]

COMMON myVars, xn_arr, weight, time_hr, Foft_hr, t_obs

c = 2.998D10 ;cm/s
kpc = 3.086D21 ;cm
day = 86400D ;s
alpha=3.    ;parcloud.alpha[0]
beta=3.     ;parcloud.beta[0]
meanE=3.    ;parcloud.meanE[0]

intensity = FLTARR(N_ELEMENTS(theta_arr))
int_n = FLTARR(N_ELEMENTS(theta_arr))
intensity_abs = FLTARR(N_ELEMENTS(theta_arr))

FOR n = 0, N_ELEMENTS(xn_arr)-1 DO BEGIN

  xn = xn_arr[n]    ;of thickness 10pc
  
  ;for 1 cloud at xn
  FOR th = 0, N_ELEMENTS(theta_arr)-1 DO BEGIN
    theta = theta_arr[th]
    dt = (theta*!PI/(180D*3600D))^2 * xn*P[0]*kpc/(2*c*day*(1-xn))
    flux_t = t_obs - dt

    t_ind = MIN(WHERE(time_hr GE flux_t))

    int_n = Foft_hr[t_ind]*$
       	    	    ((theta/(1-xn)/1000.)^(-alpha))*meanE^(-beta)/(1-xn)^2.
		    
    intensity[th] += int_n
    
    intensity_abs[th] += P[1]*weight[n]*int_n

;STOP
  ENDFOR


ENDFOR	;for all clouds



RETURN, intensity_abs
END



restore,'DATA/resparex_reid_use1.sav' ; cloud parameters
restore,'DATA/parcloud_fit.sav'
restore,'DATA/prof_rgbc67mrad5_deflare_FLAT.sav'

;wrap_simprof_cont_MC80_mde2reid,resparex_reid, cloudpar=parcloud_fit, ;outsbp=sbpout_fit, syserr=0.04;,/ps,fname='sbpfit_fit.eps'

errsys = 0.04

data_x = rad_im
data_y = NPROF_IM2C67FLATB[0,*]
data_e = NPROF_IM2C67FLATB[1,*]+NPROF_IM2C67FLATB[0,*]*errsys

ploterror,data_x, data_y, data_e,$
          yr=[1e-10,0.7e-8], xr=[0,600],/nohat, xtitle='Radius (arcsec)', $
          ytitle='ph cm!E-2!N s!E-1!N asec!E-2!N',charsize=1.3,ylog=logy

;p0 = [data_y[35], 40.D, 25.]  		   ; Initial guess
p0 = 10.
nparam = 2
myparinfo = REPLICATE({fixed:0, limited:[1,0], limits:[0.0,0.0]}, nparam)

t0=57625.  ;start of outburst
tl=57761.  ;end of dates that we use maxi data
chand=57789.4 ;date of chandra (or earlier observations)
src_time = t0+indgen(floor(tl-t0)+28)-50000.
tobs = chand-t0

hires = 5000.
time_hires = FINDGEN(hires)*(MAX(src_time)-MIN(src_time))/ hires + MIN(src_time)
Foft_hires = INTERPOL(foft, src_time, time_hires)

;PLOT,src_time,Foft, psym = 1, color = 100	    
;OPLOT, time_hires, Foft_hires, psym = 4, color = 250


p = MPFITFUN('SBP_MODEL', data_x, data_y, data_e, p0)
PRINT, p
oplot, data_x, SBP_MODEL(data_x, p), color = 200	   ; Plot model

GOTO, exit




;-----




ang = FINDGEN(600)


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



exit:

END
