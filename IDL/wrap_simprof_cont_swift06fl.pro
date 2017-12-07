pro wrap_simprof_cont_swift06fl, ps=ps, fname=namef, ylog=logy, cloudpar=parcloud, decayt=edecay, outsbp=sbpout

;This program is a wrapper that runs the single scattering flux
;accumulation for a given number of cloud and parameters for swift xrt profiles
;
; INPUTS
;
; NONE
;
; OUTPUTS
;
; A figure showing the simulated SBP over on actual data
; outsbp (optional) an array of output profiles
;
; OPTIONAL INPUTS
;
; ps: IF set plot postscript
; fname: IF set name of the postscript file
; tdif: Time of observation from the CHANDRA observation (allowing
;checking swift cases)
; ylog: IF set, plot y axis logarithmic
; cloudpar: parameter structure of clouds, also an output
; decayt: exponential decay time representing swift data
;
; USES
;
; sim_prof
;
; USED BY
;
; NONE
;
; LOGS
;
; Created by EK, Dec 2017
;

  IF NOT keyword_set(ps) THEN ps=0
  IF NOT keyword_set(namef) THEN namef='simprof_cont_swift06nfl.eps'
  IF NOT keyword_set(dift) THEN dift=0
  IF NOT keyword_set(logy) THEN logy=0
  IF NOT keyword_set(edecay) THEN edecay=0.36 ;too high
  
IF ps THEN BEGIN
   set_plot,'ps'
   device,/color
   loadct,5
   device,/encapsulated
   device, filename = namef
   device, yoffset = 2
   device, ysize = 18.
   ;device, xsize = 12.0
   !p.font=0
   device,/times
ENDIF  

!p.multi=[0,1,2]

;PLOT MAXI LIGHT CURVE 2-4 keV

sed_data=read_csv('glcbin24.0h_regbg_hv0.csv')
mjdm=sed_data.field1
m24=sed_data.field4
e24=sed_data.field5
;m410=sed_data.field6
;e410=sed_data.field7

t14d=55392.5 ;date of chandra, dift wil determine swift day

ploterror,mjdm-50000.,m24,e24,psym=4,/nohat,xr=[5180,5400],$
          xtitle='Time (MJD-50000 days)',$
          ytitle='ph cm!E-2!N s!E-1!N',charsize=1.3,/ylog,yr=[2e-4,1.]


;swift data

sdates=t14d
sflux24=[0.00033]
;sflux15=[0.00165,0.000772,0.000494,0.000124]

plotsym,0,/fill
oplot,[sdates,sdates]-50000.,[1,1]*sflux24,psym=8,color=0

;oploterror,mjdm,m410,e410,psym=5,/nohat
oplot,[t14d,t14d]-50000.,10^(!y.crange),color=0,line=2,thick=2

;arrow, chand-50000., 0.3, chand-50000., 1E-5, /data, color=0

t0=55193.  ;start of outburst
tl=55376.  ;end of dates that we use maxi data

delt=tl-t0
;oplot,[t0,t0]-50000.,!y.crange,color=0    ;show if necessary
;oplot,[tl,tl]+9-50000.,!y.crange,color=0  ;show if necessary

xx=where((mjdm GT t0) AND (mjdm LT tl))
yy=where(m24[xx]/e24[xx] GE 3.)


time=indgen(floor(tl-t0)) ;get maxi part

;m24[276]=0.02  ; average out some stupid data
;m24[279:280]=0.02


F=interpol(m24[xx[yy]],mjdm[xx[yy]],time+t0)    ; interpolate maxi data per day
;oplot,time+t0,F,color=0
;from MJD 57755 to 57780 add exponential decay

;adjust this part for swift, go 2 days before swift

edecay=0.28
exdays=floor(t14d-tl)
dd=findgen(exdays) ; decay almost all the way to Chandra data (2 days before)
Foft_add=F[n_elements(F)-1L]*exp(-dd*edecay) ; 
Foft=[F,Foft_add]
;Foft=F

;Foft=[F,0.05,0.045,0.04,0.035,0.03,0.025,0.02,0.015,0.01,0.005]
oplot,t0+indgen(floor(tl-t0)+exdays)-50000.,Foft,color=0

; parcloud is a structure that holds information up to 20 clouds
;


IF NOT keyword_set(parcloud) THEN BEGIN
   parcloud=create_struct('distance',11.58D,'x',fltarr(20),$
                          'alpha',replicate(3.,20),'beta',replicate(3.,20), $
                          'meanE',replicate(3.,20), 'thcl', fltarr(20), $
                          'weight',fltarr(20),'norm',1D-7,'ncloud',4)
ENDIF

;calculate the main cloud parameter first

;xD=11. kpc
parcloud.x[0]=0.95    ; distance to cloud/D, default 0.95
parcloud.thcl[0]=0.07 ;in kpc! default 0.07
x=parcloud.x[0]
;D=parcloud.distance
D=11D/x
print,D
alpha=parcloud.alpha[0]
beta=parcloud.beta[0]
meanE=parcloud.meanE[0]
thcl=parcloud.thcl[0]

sim_prof, x, D, alpha, beta, meanE, Foft, t14d-t0, totprof1, clth=thcl,/noplot

sbpout=dblarr(20,n_elements(totprof1))
sbpout[0,*]=totprof1
;parcloud.norm=1D-8/max(totprof1)
parcloud.norm=7.5D-9
parcloud.weight[0]=0.8 ;default 0.8          ;
sbptot=totprof1*parcloud.weight[0]*parcloud.norm

;OTHER clouds

parcloud.x[1:2]=[0.9999,0.82];,0.93,0.85,0.8,0.75,0.7,0.65]
parcloud.thcl[1:2]=[0.5,0.06];,0.05,0.05,0.05,0.05,0.05,0.3]   ;9 too thick change later
parcloud.weight[1:2]=[0.025,0.45];,0.015,0.016,0.017,0.018,0.019,0.02]

;simulate continuum
contstx=0.94
contendx=0.1
cloudn=16
;conthick=0.08
totweight=0.2
;kpc=3D21
;NH=1D22
FOR k=3,3+cloudn-1 DO BEGIN
   parcloud.x[k]=contendx+((contstx-contendx)/cloudn)*(k-3)
   parcloud.thcl[k]=((contstx-contendx)*D/cloudn)
   parcloud.weight[k]=(totweight*(1.+k*0.02))/cloudn
ENDFOR

parcloud.ncloud=cloudn+3
;x=0.99
;thcl=9.

FOR i=1, parcloud.ncloud-1 DO BEGIN

   x=parcloud.x[i]
   alpha=parcloud.alpha[i]
   beta=parcloud.beta[i]
   meanE=parcloud.meanE[i]
   thcl=parcloud.thcl[i]
   
   sim_prof, x, D, alpha, beta, meanE, Foft, t14d-t0, totprofi, clth=thcl,/noplot
   sbpout[i,*]=totprofi
   sbptot=sbptot+parcloud.norm*totprofi*parcloud.weight[i]
ENDFOR

;   norm2=.015
;   norm3=0.7
;x=0.8
;thcl=.05
;norm4=0.3

restore,'/Users/ekalemci/SWIFT/4u1630/IDL/swift06prof.sav'

rad_im=radiusverb
nprof=profverb
nprof[0,*]=profverb[0,*]/areasverb
nprof[1,*]=profverb[1,*]/areasverb


ploterror,rad_im, nflux[0,*],nflux[1,*],$
          yr=[3e-10,1e-6], xr=[0,600],/nohat, xtitle='Radius (arcsec)', $
          ytitle='ph cm!E-2!N s!E-1!N asec!E-2!N',charsize=1.3,ylog=logy

;norm=10D-9/max(totprof1)


;oplot,findgen(600),(totprof1*norm3+totprof2*norm2+totprof3*norm4)*norm,color=100,line=2,thick=4

oplot, findgen(600),sbptot, color=50, line=2, thick=4

FOR i=0, parcloud.ncloud-1 DO oplot,findgen(600),sbpout[i,*]*parcloud.norm*parcloud.weight[i],color=100,line=1,thick=2


;plotpsf

psf_norm=mintpsf*nflux[0,0]*85./(mintpsf[0]*areasverb)
oplot, rad_im, psf_norm,psym=10, line=2, color=0

psf_norm_int=interpol(psf_norm,rad_im,findgen(600))

oplot, findgen(600),sbptot+psf_norm_int, color=100, line=2, thick=4


IF ps THEN BEGIN
   device,/close
   set_plot,'x'
ENDIF


END
