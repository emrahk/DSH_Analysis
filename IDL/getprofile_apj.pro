pro getprofile, fluxim, radius, prof, areas, $
                ps=ps, fname=namef, $
                spos=poss, annum=num_an, radm=mrad, $
                corb=corb, backpos=posback, backr=rback, brate=rateb, $
                simrad=s_radius, simareas=s_areas, simprof=s_prof,$
                yrn=nyr, xrn=nxr, dopsf=dopsf, elimps=elimps,cprof=profc,$
                expthresh=threshexp, diag=diag, chip=uchip, expo_prof=prof_expo
                  

;This program reads a Chandra fluxed image and plots the radial
;profile of the scattering halo around a source

;INPUT
;
;fluxim: name of the fluxed image file
;
;OUTPUT
;
;radius : to ease recreating plots radii are provided in arcseconds 
;prof: profile, a structure with radii, count rates and in the future
;fluxes
;areas: areas in profile to obtain total counts in each annuli
;
;OPTIONAL INPUTS
;
;ps : IF set, postscript output
;fname: name of the postscript file for profile
;spos : source position in sky pixels, best to obtain from chandra
;reposses, using image coordinates simplifies things
;  
;radm : If set, multiples of this radius in arcsec will be used to
;           create profile, default=2''
;annum : number of annulis with given radm difference
;dopsf : if set use chart simulations to oplot psf
;
;corb: IF set do background correction  
;backpos : center position to obtain background/arcsec^2, in physical units!
;rpos: side of square region (to simplify calculations)
;brate: if given use this background counts per asec2
;
;elimps: if set eliminate point sources
;
;yrn : new yrange for the profile plot
;xrn : new xrange for the profile plot
;
;expthresh: exposuremap threshold
;
;diag: if set put more information about diagnostics
;
;chip: IF set, only use the given chip, only single chip is possible
;expo_prof: If given, provide profiles of the exposure map
;
;OPTIONAL OUTPUTS
;
;simrad : simulation radius
;simareas : simulation areas
;simprof : simulated profiles
;cprof   : count profile to compare with previous work
;
; LOGS
;
; Created, EK, Aug 2017
;
; adding capacity to use different chips
;
  

IF NOT keyword_set(ps) THEN ps=0
IF NOT keyword_set(namef) THEN namef='profile.eps'
;IF NOT keyword_set(poss) THEN poss=[638.177,1905.846] ;obtained from ds9
IF NOT keyword_set(poss) THEN poss=[638,1906] ;obtained from ds9, integer
IF NOT keyword_set(mrad) THEN mrad=2.
IF NOT keyword_set(posback) THEN posback=[930,1285]
IF NOT keyword_set(rback) THEN rback=80.
IF NOT keyword_set(num_an) THEN num_an=50
IF NOT keyword_set(corb) THEN corb=0
IF NOT keyword_set(dopsf) THEN dopsf=0
IF NOT keyword_set(elimps) THEN elimps=0
IF NOT keyword_set(threshexp) THEN threshexp=0.02
IF NOT keyword_set(diag) THEN diag=0
IF NOT keyword_set(uchip) THEN uchip=0

;note the unit for fluxed image:
;photon/cm2/s/pixel

;sky coordinates to arcseconds (ok for small angles...)
pix=0.492 ;arcsec

loadct,5
IF ps THEN begin
   set_plot,'ps'
   device,/color
   device,/encapsulated
   device,filename=namef
   device,ysize=10
   device,xsize=15
   !p.font=0
   device,/times
ENDIF else device,decomposed=0

;READ image

im=readfits(fluxim, hd)
xmax=sxpar(hd,'NAXIS1')
ymax=sxpar(hd,'NAXIS2')
ltv1=sxpar(hd,'LTV1')
ltv2=sxpar(hd,'LTV2')

;set chip number

IF ((uchip GT 8) OR (uchip LT 5)) THEN BEGIN
   uchip=0
   print,'chip id must be between 5-8, the source is on 7, setting to 0'
ENDIF ELSE BEGIN
     xc=loadcol('chips.reg','x')
     yc=loadcol('chips.reg','y')
     ccd_id=loadcol('chips.reg','CCD_ID')
     chipin=where(ccd_id EQ uchip)
     xfin=where(finite(xc[*,chipin[0]]) EQ 1)
     yfin=where(finite(yc[*,chipin[0]]) EQ 1)
     xcl=xc[xfin,chipin[0]]+ltv1
     ycl=yc[xfin,chipin[0]]+ltv2
     chipobj=Obj_New('IDLanROI', xcl,ycl)
     chipobjps=Obj_New('IDLanROI', xcl-ltv1,ycl-ltv2)
ENDELSE

;READ exposure map for correctly cutting out out of borders

imgpos=strpos(fluxim,'.img')
basname=strmid(fluxim,0,imgpos-4)
expomap=readfits(basname+'thresh.expmap')
IF uchip EQ 0 THEN expoc=getareacorrection(expomap,poss,threshexp, $
                                           uchip, incexpo_prof=incprof_expo) $
ELSE expoc=getareacorrection(expomap,poss,threshexp, uchip, $
                             objchip=chipobj,incexpo_prof=incprof_expo) 

;xx=where(im GT 0.) ;use only existing events in the image
xx=where((im GT 0.) AND (expomap GT (max(expomap)*threshexp))) ;use existing events with exposure map threshold larger than the given value

y=xx/xmax
x=xx-y*xmax

;chip cut
IF uchip NE 0 THEN BEGIN
   cpoints=chipobj->ContainsPoints(x, y)
   cpin=where(cpoints EQ 1)
   x=x[cpin]
   y=y[cpin]
ENDIF


;get the fluxes and exposures

flux=im[x,y]                    ;only use non-zero points in the image
exponz=expomap[x,y]             ;non zero exposure map values to get counts and errors

;recenter images xr, yr recentered positions

xr=x-poss[0]
yr=y-poss[1]

;determine background

;make poss integer

IF corb THEN BEGIN
   
   ;check if the background region is within bounds
   IF ((posback[0]-rback) LT 0) THEN BEGIN
      print, 'background region x lower limit is out of bounds'
      x0back=0
   ENDIF ELSE x0back=posback[0]-rback

   IF ((posback[0]+rback) GE xmax) THEN BEGIN
      print, 'background region x upper limit is out of bounds'
      x1back=xmax-1
   ENDIF ELSE x1back=posback[0]+rback

      IF ((posback[1]-rback) LT 0) THEN BEGIN
      print, 'background region y lower limit is out of bounds'
      y0back=0
   ENDIF ELSE y0back=posback[1]-rback

   IF ((posback[1]+rback) GE ymax) THEN BEGIN
      print, 'background region y upper limit is out of bounds'
      y1back=ymax-1
   ENDIF ELSE y1back=posback[1]+rback

      
   totb=total(im[x0back:x1back,y0back:y1back])
   backarea=(x1back-x0back)*(y1back-y0back)*pix^2.
   normb=totb/backarea
   normb_err=sqrt(totb)/backarea

  IF keyword_set(rateb) THEN BEGIN
      normb=rateb
      normb_err=sqrt(normb)*0. ;fix this later
   ENDIF

   print, 'normalized background per arcsec^2: ', normb, normb_err
  
ENDIF ELSE BEGIN
   normb=0.
   normb_err=0.
ENDELSE
   
;convert to arcsec

ax=xr*pix
ay=yr*pix

;full profile with no cut on pile up

incpr=dblarr(2,num_an)  ;total counts at each radii
prof=dblarr(2,num_an)   ;total counts at each annulus (except first element is circle)
areas=fltarr(num_an)  ;area of each annulus (and the first circle)

radius=findgen(num_an)*mrad+(mrad/2.) ;radius refers to mid of annulus

profc=dblarr(2,num_an)
incprofc=dblarr(2,num_an)
prof_expo=dblarr(num_an)


;first element (usually highly piled up)

areas[0]=!pi*mrad^2.
xx=where((ax^2+ay^2) LE (mrad^2.),nonzero) ;number of counts inside first circle

IF nonzero NE 0 THEN BEGIN
   incprofc[0,0]=total(flux[xx]*exponz[xx])  ;get counts
   errfrac=sqrt(incprofc[0,0])/incprofc[0,0]
   incprofc[1,0]=sqrt(incprofc[0,0])
   profc[*,0]=incprofc[*,0]
   incpr[0,0]=total(flux[xx])
   incpr[1,0]=incpr[0,0]*errfrac
   prof_expo[0]=incprof_expo[floor((mrad/pix)+0.5)]
   IF diag THEN BEGIN
      print, 'innermost circle counts', nonzero, errfrac
      plot, histogram(flux[xx], min=0, binsize=1e-8), psym=10
      stop
   ENDIF
   
ENDIF


areas[0]=areas[0]*expoc[floor((mrad/pix)+0.5)]

prof[0,0]=incpr[0,0]-(normb*areas[0]) ;background corrected rate
prof[1,0]=incpr[1,0]+(normb_err*areas[0]) ;background corrected rate

;annuli
;areas1=areas
FOR i=1,num_an-1 DO BEGIN
   rad1=mrad*i
   rad2=mrad*(i+1)              ; outer radius
   xx=where((ax^2+ay^2) LE (rad2^2.),nonzero)
   IF nonzero NE 0 THEN BEGIN
      incprofc[0,i]=total(flux[xx]*exponz[xx])
      incprofc[1,i]=sqrt(incprofc[0,i]) 
      errfracall=incprofc[1,i]/incprofc[0,i]
      profc[0,i]=incprofc[0,i]-incprofc[0,i-1]
      profc[1,i]=sqrt(profc[0,i])
      errfracring=sqrt(profc[0,i])/profc[0,i]
      incpr[0,i]=total(flux[xx])
      incpr[1,i]=incpr[0,0]*errfracall
   ENDIF ELSE BEGIN
      errfracring=0.
      errfracall=0.
   ENDELSE
   
;   areas[i]=(!PI*(rad2^2.)*expoc[floor((rad2/pix)+0.5)])-$
;            (!PI*(rad1^2.)*expoc[floor((rad1/pix)+0.5)])
   areas[i]=(!PI*(rad2^2.)*interpolate(expoc,rad2/pix))-$
            (!PI*(rad1^2.)*interpolate(expoc,rad1/pix))
   prof[0,i]=(incpr[0,i]-incpr[0,i-1])-(normb*areas[i]) 
   prof[1,i]=prof[0,i]*errfracring+(normb_err*areas[i])

   prof_expo[i]=interpolate(incprof_expo,rad2/pix)-interpolate(incprof_expo,rad1/pix)
   
   IF diag THEN BEGIN
      print, i, 'radius and area: ', rad1,areas[i]
      print, i, 'th circle counts: ', nonzero, profc[0,i], incpr[0,i], errfracall,errfracring
      print, i, 'exposure profile: ', prof[1,i]
      IF xx[0] NE -1 THEN plot, histogram(flux[xx], min=0, binsize=1e-8),psym=10,/ylog,yr=[0.1,1e5]
      stop
   ENDIF
      
ENDFOR

;eliminate point sources

IF elimps THEN BEGIN
   inpfile='full_broad_outfile_2.fits'
   detposs=fltarr(2)
   detposs[0]=poss[0]-ltv1
   detposs[1]=poss[1]-ltv2
   elimpointsrc, inpfile, xr, yr, detposs, mrad, num_an, psc, $
                 areadif, diag=diag, psfr=5, sthresh=4.,inpflux=flux, $
                 outflux=fluxsrc, chip=uchip, objchip=chipobjps
   profc[0,*]=profc[0,*]-psc
   prof[0,*]=prof[0,*]-fluxsrc
   areas=areas-areadif
ENDIF

;normalized profile
nprof=fltarr(2,num_an)
nprof[0,*]=prof[0,*]/areas
nprof[1,*]=prof[1,*]/areas

;psf part

;method 1, match non-piled up part to expected counts
; calculate average energy

;destroy object
IF uchip NE 0 THEN BEGIN
   Obj_Destroy,chipobj
   Obj_Destroy,chipobjps
ENDIF

IF dopsf THEN BEGIN
   aven=2.3
   avcens=[1.,3.5,7.5]               ;available chart energies, should do 2.3
   mindelta=min(abs(avcens-aven), mind) ;mind: index of the desired energy
   getsimprofile, outsim                ; get the profiles

   print, 'index found :',mind

   s_radius=outsim[mind].radius
   s_prof=outsim[mind].prof
   s_areas=outsim[mind].areas

;rough normalization

;ronorm=total(prof[0,0:1])/total(outsim[mind].prof[0,0:1])

   ronorm=fltarr(3)
   ronorm[mind]=total(prof[0,0])/total(outsim[mind].prof[0,0])
;be careful abount energy range 5-10
   ronorm[2]=total(prof[0,0:1])/total(outsim[mind].prof[0,0:1])

ENDIF

dopucorr=0

pucorr=[1.,1.,1.]

IF dopucorr THEN BEGIN
   
;pileupcorr

;2-5 keV
;heg model predict : 0.216 after pile
;1.95 / 0.21 - pile up correction 87%
;1.38 pimms predict without pile up 

;5-10 keV
;heg model predicts: 0.05 after pile

pucorr[1]=1.38/0.216
pucorr[2]=2. ;???

ENDIF

;plot results

IF NOT keyword_set(nxr) THEN nxr=[0,max(radius)]
multiplot, [1,2], mxtitle='Radius (arcsec)',mxtitsize=1.2
IF NOT keyword_set(nyr) THEN nyr=[min(prof[0,*]-prof[1,*])*0.8,max(prof[0,*]+prof[1,*])*1.1]
ploterror,radius,prof[0,*],prof[1,*],psym=10,/ylog,ytitle='cts',/nohat,$
          yr=nyr*100.,/ystyle, xr=nxr
IF dopsf THEN oplot,outsim[mind].radius, outsim[mind].prof[0,*]*ronorm[mind]*pucorr[mind], psym=10, line=2, color=0
;oploterror, radius, normb*areas,normb_err*areas,/nohat 
multiplot
ploterror,radius,nprof[0,*],nprof[1,*],psym=10,/ylog,ytitle='cts arcsec!E-2!N',$
          /nohat, chars=1.2, yr=nyr, xr=nxr
IF dopsf THEN oplot,outsim[mind].radius, outsim[mind].nprof[0,*]*ronorm[mind]*pucorr[mind],psym=10, line=2, color=0
oplot,!x.crange,[normb,normb],line=1
multiplot,/default



IF ps THEN BEGIN
   device,/close
ENDIF


END

function getareacorrection, expomap, poss, threshexp, uchip, objchip=chipobj, incexpo_prof=incprof_expo

;
; expomap: exposure map image
; poss: position of the source
; threshexp: threshold to accept data
; uchip: chip nimber to use, if 0, use entire detector
; objchip: chip boundaries if uchip is different than zero
; incexpo_prof: incremantal exposure map profile (grows with area)

  
;  out of chip boundary and exposure area 0 pixels should not
;  contribute to the area
  yy=where(expomap GT ((max(expomap)*threshexp)))
  
  sz=size(expomap)
  yexp=(yy/sz[1]);-poss[1]
  xexp=(yy-yexp*sz[1]);-poss[0]

  ;get chip correction  
  IF uchip NE 0 THEN BEGIN
     conpoint=chipobj->ContainsPoints(xexp, yexp)
     zz=where(conpoint EQ 1)
     xexp=xexp[zz]
     yexp=yexp[zz]
  ENDIF
  
  yexpc=yexp-poss[1]
  xexpc=xexp-poss[0]
  rexp=sqrt(xexpc^2+yexpc^2)
  expoarfrac=dblarr(1500)       ;returned, active area/pi r^2
  incprof_expo=dblarr(1500); incremental exposure map profile for the given radius
  ;assume first 100 is always ok
  ;expoarfrac[0:99]=1. ;then the pir2 counting ratio is very close to 1
  FOR i=0, 1499 DO BEGIN
     rlim=where(rexp LT i,ctrlim) ;find all elements less than given radius
     
     IF (ctrlim NE 0) THEN BEGIN
        IF (i GE 180) THEN expoarfrac[i]=ctrlim/(!PI*double(i)^2.) $
        ELSE expoarfrac[i]=1  ; below 180 exposure area is always > expthresh 
        incprof_expo[i]=total(expomap[xexp[rlim],yexp[rlim]]/max(expomap))
        diagx=0
        IF (diagx AND (i GT 800)) THEN begin
           cix=findgen(3600)*2*!PI/3600.
           plot,i*cos(cix),i*sin(cix)
           oplot,xexpc[rlim],yexpc[rlim],psym=3
           stop
        ENDIF
        
     ENDIF
     
  ENDFOR
  return,expoarfrac
END

;rade=mrad/pix
                                ; this is really stupid
;radeint=floor(rade)+1
;nonzero=0                       ;for error calculation
;posexpo=0                       ;to obtain positive exposure map
;allpixs=0
;FOR k=-radeint, radeint DO BEGIN
;    FOR l=-radeint, radeint DO BEGIN
;       IF (k^2+l^2) LT rade^2 THEN BEGIN
;          allpixs=allpixs+1L
;          IF ((poss[0]+k GE 0) AND (poss[0]+k LT xmax) AND (poss[1]+l GE 0) AND (poss[1]+l LT ymax)) THEN BEGIN
;             IF im[poss[0]+k,poss[1]+l] GT 0. THEN BEGIN
;                incpr[0,0]=incpr[0,0]+im[poss[0]+k,poss[1]+l]
;                nonzero=nonzero+1L
;                posexpo=posexpo+1L
;             ENDIF ELSE BEGIN
;                IF expomap[poss[0]+k,poss[1]+l] GT 0. THEN posexpo=posexpo+1L
;             ENDELSE
;          ENDIF
;       ENDIF
;    ENDFOR
; ENDFOR

;print,'allpixs, nonzero, posexpo:',allpixs, nonzero, posexpo
;incpr[1,0]=incpr[0,0]*sqrt(nonzero)/nonzero

;areafac=posexpo/allpixs
;areas[0]=areas[0]*areafac
