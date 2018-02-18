pro wrap_simprof_cont_MC80nd_mde2reid, inpar, ps=ps, fname=namef, tdif=dift, $
                                 ylog=logy, cloudpar=parcloud, decayt=edecay, $
                                 outsbp=sbpout, syserr=errsys, totsbp=sbptot, $
                                 totchi2=chi2tot, silent=silent

;This program is a wrapper that runs the single scattering flux
;accumulation for a given number of cloud and parameters 
;
; INPUTS
;
; inpar: structure that includes cloud paremeters obtained from co_gaussfits
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
; syserr: If set, add errsys*flux to the error
; totsbp: total brightness profile model
; silent: if set do not print info on the screen
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
  IF NOT keyword_set(namef) THEN namef='simprof_cont_MC80_mde2reid.eps'
  IF NOT keyword_set(dift) THEN dift=0
  IF NOT keyword_set(logy) THEN logy=0
  IF NOT keyword_set(edecay) THEN edecay=0.34 ;too high
  IF NOT keyword_set(errsys) THEN errsys=0.
  IF NOT keyword_set(silent) THEN silent=0
  
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

chand=57789.4-dift ;date of chandra (or earlier observations)

ploterror,mjdm-50000.,m24,e24,psym=4,/nohat,xr=[7600,7800],$
          xtitle='Time (MJD-50000 days)',$
          ytitle='ph cm!E-2!N s!E-1!N',charsize=1.3,/ylog,yr=[1e-5,1.]


;swift data

sdates=[57769.845,57771.775,57773.237,57777.35]
sflux24=[0.00110, 0.000471,0.000305,0.000076]
sflux15=[0.00165,0.000772,0.000494,0.000124]

plotsym,0,/fill
oplot,sdates-50000.,sflux24,psym=8,color=0

;oploterror,mjdm,m410,e410,psym=5,/nohat
oplot,[chand,chand]-50000.,10^(!y.crange),color=0,line=2,thick=2

;arrow, chand-50000., 0.3, chand-50000., 1E-5, /data, color=0

t0=57625.  ;start of outburst
tl=57761.  ;end of dates that we use maxi data

delt=tl-t0
;oplot,[t0,t0]-50000.,!y.crange,color=0    ;show if necessary
;oplot,[tl,tl]+9-50000.,!y.crange,color=0  ;show if necessary

xx=where((mjdm GT t0) AND (mjdm LT tl))

time=indgen(floor(tl-t0)) ;get maxi part

m24[2157:2159]=0.285  ; average out some stupid data
m24[2182:2183]=0.18
F=interpol(m24[xx],mjdm[xx],time+t0) ; interpolate maxi data per day
;oplot,time+t0,F,color=0
;from MJD 57755 to 57780 add exponential decay

dd=findgen(28)+1 ; decay almost all the way to Chandra data (2 days before)
Foft_add=F[n_elements(F)-1L]*exp(-dd*edecay) ; 
Foft=[F,Foft_add]

;Foft=[F,0.05,0.045,0.04,0.035,0.03,0.025,0.02,0.015,0.01,0.005]
oplot,t0+indgen(floor(tl-t0)+28)-50000.,Foft,color=0

; parcloud is a structure that holds information up to 50 clouds
; IF given use parameters inside, if not enter defaults

maxind=3                        ;main cloud, MC -79, magic number!
incloudn=n_elements(inpar.cname)

pass=0
IF keyword_set(parcloud) THEN BEGIN
   checkres=isa(parcloud,'struct')
   IF checkres EQ 1 THEN BEGIN
      sizeres=size(parcloud)
      IF sizeres[2] EQ 8 THEN pass=1
   ENDIF
ENDIF

IF pass EQ 0 THEN BEGIN
   IF NOT silent THEN  print,'cloudpar is not given, or not a structure with the right parameters'
   ncloud=50
   parcloud=create_struct('distance',0D,'x',fltarr(ncloud),$
                          'alpha',replicate(3.,ncloud),$
                          'beta',replicate(3.,ncloud), $
                          'meanE',replicate(3.,ncloud),'thcl', fltarr(ncloud),$
                          'weight',fltarr(ncloud),'norm',0D,'ncloud',ncloud,$
                         'absf',fltarr(ncloud),'sbpout',dblarr(ncloud,10000L),$
                         'contweight',0D,'ndcloud',4)
 
   ;free parameters
   ;place MC 79 to first index:
   parcloud.x[0]=0.985          ; distance to cloud/D, default 0.95?
   parcloud.thcl[0]=0.08        ;in kpc! default 0.07?   
   parcloud.weight[0]=7.9 ;weight of 0th cloud, not sure why we have this together with normalization, mybe weight [0] should be 1 and overall normalization should be adjusted by norm.
   parcloud.contweight=13.5D
   parcloud.norm=2D-8
                                ;adjusted parameters, place other clouds
   parcloud.distance=inpar.estdist[0,maxind]/parcloud.x[0]

    ;Populate other clouds
    ;determine clouds to use

   relx=inpar.estdist[0,*]/inpar.estdist[0,maxind] ;relative distance to the main cloud
;frclds=where(relx LT 1.) ;MC 66 may be a problem that would be handled later.

   frclds=where(inpar.estdist[0,*]/parcloud.distance LT 1.)
   frclds=frclds(where(frclds NE maxind)) ;remove main cloud
   parcloud.ndcloud=n_elements(frclds+1)
   
;frclds=[0,1,2,4,5,8,9,10] ;1630-47 is possible in MC -66 [5] handled later, 6,7 probably far away

   FOR i=1, n_elements(frclds) DO BEGIN
      parcloud.x[i]=parcloud.x[0]*relx[frclds[i-1]] 
;   thcli=thcl*inpar.csize[0,frclds[i-1]]/inpar.csize[0,maxind]>0.01
      parcloud.thcl[i]= inpar.csize[0,frclds[i-1]]/1000. ;to convert to kpc
      parcloud.weight[i]=parcloud.weight[0]*inpar.gint[0,frclds[i-1]]/inpar.gint[0,maxind] ;,0.015,0.016,0.017,0.018,0.019,0.02] eight is proportional to integrated velocity
      IF NOT SILENT THEN print,strtrim(string(i),1)+' '+$
                    inpar.cname[frclds[i-1]]+' at '+$
         string(inpar.estdist[0,frclds[i-1]])+$
         ' with thickness'+string(parcloud.thcl[i])
   ENDFOR


   ;populate continuum

   contstx=0.99                 ; starts from inside cloud+cloud thickness
   contendx=0.1 ; ends relatively close to observer, though emission should go down significantly close to observer? check
   cloudnc=parcloud.ncloud-parcloud.ndcloud ; number of continuum clouds to simulate
   totweight=parcloud.contweight
   k0=parcloud.ndcloud      ;or 2 depending on Mc -66 check
   FOR k=k0,parcloud.ncloud-1 DO BEGIN
      parcloud.x[k]=contendx+((contstx-contendx)/cloudnc)*(k-k0+1)
      parcloud.thcl[k]=((contstx-contendx)*parcloud.distance/cloudnc)
      parcloud.weight[k]=totweight/cloudnc
   ENDFOR
   
ENDIF


  
;calculate the main cloud parameter first
;there are some magic numbers, but this is a wrapper

x=parcloud.x[0]
D=parcloud.distance
IF NOT SILENT THEN print, 'Distance to source is: '+string(strtrim(D,1))+' kpc'
alpha=parcloud.alpha[0]
beta=parcloud.beta[0]
meanE=parcloud.meanE[0]
thcl=parcloud.thcl[0]

;simulate profile from the first cloud, total emission
sim_prof, x, D, alpha, beta, meanE, Foft, chand-t0, totprof1, clth=thcl,$
          /noplot,silent=silent

sbpout=dblarr(parcloud.ncloud,n_elements(totprof1)) ;create an output profile array
sbpout[0,*]=totprof1


;normalization
;exp(-1.0792) total extinction
totw=total(parcloud.weight) ;ok this includes continuum
perwfac=-1.2/totw ;per weight unit extinction (was 1.079)
xx=where(parcloud.x*D LE parcloud.x[0]*D) ;start with main cloud, 
totw0=total(parcloud.weight[xx])
parcloud.absf[0]=exp(totw0*perwfac)
sbptot=totprof1*parcloud.weight[0]*parcloud.norm*parcloud.absf[0] ;keep track of total profile


FOR i=1, parcloud.ncloud-1 DO BEGIN
   
   x=parcloud.x[i]
   alpha=parcloud.alpha[i]
   beta=parcloud.beta[i]
   meanE=parcloud.meanE[i]
   thcl=parcloud.thcl[i]
   sim_prof, x, D, alpha, beta, meanE, Foft, chand-t0, totprofi, clth=thcl, $
             /noplot,silent=silent
   sbpout[i,*]=totprofi
                                ;print,parcloud.x[k]*D,absfac
   xx=where(parcloud.x*D LE parcloud.x[i]*D)
   totwi=total(parcloud.weight[xx])
   parcloud.absf[i]=exp(totwi*perwfac)
   sbptot=sbptot+parcloud.norm*totprofi*parcloud.weight[i]*parcloud.absf[i]
ENDFOR

;   norm2=.015
;   norm3=0.7
;x=0.8
;thcl=.05
;norm4=0.3

restore,'../CHIPvsBACKG5keV/prof_rgbc67mrad5_deflare_FLAT.sav'

ploterror,rad_im, NPROF_IM2C67FLATB[0,*],$
          NPROF_IM2C67FLATB[1,*]+NPROF_IM2C67FLATB[0,*]*errsys,$
          yr=[1e-10,0.7e-8], xr=[0,600],/nohat, xtitle='Radius (arcsec)', $
          ytitle='ph cm!E-2!N s!E-1!N asec!E-2!N',charsize=1.3,ylog=logy

;yr=[1e-10,1e-8]. xr=[0,600]
dopsf=0
IF dopsf THEN BEGIN ; we can add chart simulations at different energied and do this as well
   getsimprofile, outsim, radm=5
   s_radius=outsim[1].radius
   s_nprof=outsim[1].nprof
   s_areas=outsim[1].areas
   ronorm=NPROF_IM2C67FLATB[0,0]/outsim[1].nprof[0,0]
   psf_norm=s_nprof[0,*]*ronorm
   psf_norm_int=interpol(psf_norm,s_radius,findgen(600))
   oplot, findgen(600),sbptot+psf_norm_int, color=100, line=2, thick=3
ENDIF
   
;oplot,s_radius,psf_norm,color=0,line=2,psym=10

oplot, findgen(600),sbptot, color=50, line=2, thick=4

FOR i=0,parcloud.ndcloud-1 DO oplot,findgen(600),sbpout[i,*]*parcloud.norm*parcloud.weight[i],color=100,line=1,thick=2

;continds=[parcloud.ndcloud:49]
contspb=fltarr(10000L)

FOR i=parcloud.ndcloud, parcloud.ncloud-1 DO BEGIN
   contspb=contspb+sbpout[i,*]*parcloud.norm*parcloud.weight[i]
ENDFOR
   
oplot,findgen(600),contspb[0:599],color=100,line=3,thick=2

;chi2 between 110 and 250
yy=where((rad_im GE 110.) AND (rad_im LE 250.))
sbmatch=interpol(sbptot, findgen(10000L),rad_im(yy))
dof=float(n_elements(yy)-3)
toterr=NPROF_IM2C67FLATB[1,yy]+NPROF_TOTC67FLATB[0,yy]*errsys
chi=((NPROF_IM2C67FLATB[0,yy]-sbmatch)/toterr)
redchi=total(((NPROF_IM2C67FLATB[0,yy]-sbmatch)/toterr)^2.)/dof
IF NOT SILENT THEN print, 'Reduced chi2: ',redchi
chi2tot=total(((NPROF_IM2C67FLATB[0,yy]-sbmatch)/toterr)^2.)

IF ps THEN BEGIN
   device,/close
   set_plot,'x'
ENDIF


END
