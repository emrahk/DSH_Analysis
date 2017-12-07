pro sim_prof, x, D, alpha, beta, meanE, Foft, tobs, totprof, clth=thcl, noplot=noplot

;This program is written to simulate SBP for given dust parameters
;
;
;INPUTS
;
; x: ratio of distance to dust cloud/distance to source (kpc)
; D: distance to the source (kpc)
; alpha: scattering cross section scattering angle dependence
; beta: scattering cross section energy dependence
; Foft: flux as a function of time, starting from t=0 (days)
; tobs: time of observation with respect to given t=0 (days)
; 
; OUTPUTS
;
; prof: output profile
;
; OPTIONAL INPUTS
;
; clth: thickness of the cloud in kpc
;
; USES
;
; NONE
;
; USED BY
;
; NONE
;
;LOGS
;
; Created by EK, November 2017
;

IF NOT keyword_set(thcl) THEN thcl_ind=1 ELSE thcl_ind=floor(thcl/0.01)
IF NOT keyword_set(noplot) THEN noplot=0

;  theta=findgen(600) ;in asec
;  theta_sc=theta/(1-x)
;  Omega=(theta_sc)^(-alpha)*meanE^(-beta)
  prof=dblarr(thcl_ind,2,n_elements(Foft))
  print, 'Cloud at ',x*D
  FOR j=1,thcl_ind DO BEGIN
     xn=x-j*0.01/D
     FOR i=0, n_elements(Foft)-1 DO BEGIN
        dt=tobs-i
        theta=100d*sqrt(dt*(1-xn)*11./(1.4*(xn*D)))
        theta_sc=theta/(1-xn)
        Omega=((theta_sc/1000.)^(-alpha))*meanE^(-beta)
        intensity=Foft[i]*Omega/(1-xn)^2.
        prof[j-1,0,i]=theta
        prof[j-1,1,i]=intensity
     ENDFOR
  ENDFOR

  print,xn,D,xn*D
  
 IF NOT noplot THEN BEGIN 
  !p.multi=[0,1,2]
  plot, prof[0,0,*],prof[0,1,*], xr=[0,600]
ENDIF
 
  IF thcl_ind GT 1 THEN BEGIN
     totprof=fltarr(10000L)
     xint=interpol(prof[0,1,*],prof[0,0,*],findgen(10000L))
     minmax=[floor(min(prof[0,0,*])),floor(max(prof[0,0,*]))]
     totprof[minmax[0]:minmax[1]]=xint[minmax[0]:minmax[1]]
     FOR k=1,thcl_ind-1 DO BEGIN
        oplot, prof[k,0,*], prof[k,1,*]
        xint=interpol(prof[k,1,*],prof[k,0,*],findgen(10000L))
        minmax=[floor(min(prof[k,0,*])),floor(max(prof[k,0,*]))]
        totprof[minmax[0]:minmax[1]]=totprof[minmax[0]:minmax[1]]+xint[minmax[0]:minmax[1]]
     ENDFOR
IF NOT noplot THEN plot,findgen(10000L),totprof/k,line=2
  ENDIF
  
IF NOT noplot THEN !p.multi=0
END

     
