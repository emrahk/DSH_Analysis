pro elimpointsrc, inpfile, x, y, sposs, mrad, numann, psc, areadif, $
                  sthresh=threshs, spix=pixs, psfr=rpsf, method=xmethod, $
                  removeall=removeall, diag=diag, inpflux=fluxinp, $
                  outflux=fluxout, chip=uchip, objchip=chipobj

;This program removes circular regions centered on detected point
;sources in images and calculates new areas in each annulus

;INPUTS
;  inpfile: input fits file with detected source information
;  x: x positions in image, physical
;  y: y positions in image, physical
;  sposs: position of central source in physical pixels
;  mrad: size of each annulus
;  numann: number of annuli
;
;OUTPUTS
;  psc: source counts in each annulus
;  areadif: area to subtract in each annulus
;
;OPTIONAL INPUTS
;
; sthresh: threshold for source significance, default 3
; spix: pixel size default Chandra
; psfr: cut-out radius in pixels, default 5 for Chandra
; removeall: IF set remove all detected sources, including the central source
;  method: TBD
; diag: IF set do diagnostics with plots, must set to 0 inside profile programs
; inpflux: if given, use image to calculate counts
; chip: if not zero, make a cut on chip
; objchip: if chip is given, this object is required to include
; sources in the chip
; 
; USED BY
;
; getprofile, getswiftprofile
;
;USES
;
; NONE
;
; LOGS
;
; Created by EK, July 2017
;
;

;set optional parameters
;

  IF NOT keyword_set(threshs) THEN threshs=3.
  IF NOT keyword_set(pixs) THEN pixs=0.492
  IF NOT keyword_set(removeall) THEN removeall=0
  IF NOT keyword_set(rpsf) THEN rpsf=5
  IF NOT keyword_set(diag) THEN diag=0

  
;get detected source data

  sxo=loadcol(inpfile,'X')
  syo=loadcol(inpfile,'Y')
  sig=loadcol(inpfile,'SRC_SIGNIFICANCE')
  psfsize=loadcol(inpfile,'PSF_SIZE')
  
; turn into radius

  sx=sxo-sposs[0]
  sy=syo-sposs[1]
  rad=sqrt(sx^2.+sy^2.)

; create results

  psc=fltarr(numann)
  areadif=fltarr(numann)
  IF keyword_set(fluxinp) THEN BEGIN
     useim=1
     fluxout=fltarr(numann)
  ENDIF
  

;go through all detected sources

;make a cut on chip

IF uchip NE 0 THEN BEGIN
   cpoints=chipobj->ContainsPoints(sxo, syo)
   cpin=where(cpoints EQ 1)
   sig=sig[cpin]
   rad=rad[cpin]
   sx=sx[cpin]
   sy=sy[cpin]
   psfsize=psfsize[cpin]
ENDIF

;make a cut on source significance and radii
  xx=where((sig GT threshs) AND (rad LE mrad*numann/pixs),nums)

  IF nums NE 0 THEN BEGIN
     
     IF removeall THEN st_ind=0 ELSE st_ind=1

     IF diag THEN print, strtrim(string(nums-st_ind),1)+ $
                         ' point sources will be removed'
     ;Make a sort on point source significance
     
     rad=rad[xx[reverse(sort(sig[xx]))]]
     sx=sx[xx[reverse(sort(sig[xx]))]]
     sy=sy[xx[reverse(sort(sig[xx]))]]
     psfsize=psfsize[xx[reverse(sort(sig[xx]))]]/pixs ;in pixels

     FOR i=st_ind, nums-1 DO BEGIN
         ;get the index of first annulus before the center
        j=floor((rad[i]*pixs/mrad)) ; center is between jth and j+1th circle
                                ;determine how many annulus the circle spreads over
        usepsf=rpsf>psfsize[i]
        jmax=floor(((rad[i]+usepsf)*pixs/mrad))+1
        jmin=floor(((rad[i]-usepsf)*pixs/mrad))
        ctind=where((((x-sx[i])^2.+(y-sy[i])^2.) LE usepsf^2.),ctsall) ;all pts
IF ctsall GT 0 THEN BEGIN
        IF (useim AND (ctsall GT 0)) THEN fluxall=fluxinp[ctind]
        angs=findgen(360)*!PI/180.
        radpl1=jmin*mrad/pixs
        
        IF diag THEN BEGIN
           print, 'position: ',sx[i]+sposs[0],sy[i]+sposs[1]
           print, 'psfsize wavdetect: ',psfsize[i], 'used: ',usepsf
           plot, x, y, $
                 xr=[min(x[ctind])-30,max(x[ctind])+30],$
                 yr=[min(y[ctind])-30,max(y[ctind])+30],psym=1
           
           IF ctsall GT 0 THEN oplot,x[ctind],y[ctind],psym=4
           oplot,radpl1*cos(angs),radpl1*sin(angs)
           plind=findgen(360)*!PI/180.
           oplot,sx[i]+usepsf*cos(plind),sy[i]+usepsf*sin(plind)
        ENDIF

        totarea=!PI*(usepsf^2.)*(pixs^2.)
        IF ((jmin EQ j) AND (jmax EQ j+1)) THEN BEGIN
           ;all psf in one annulus
           psc[j]=psc[j]+ctsall
           IF useim THEN fluxout[j]=fluxout[j]+total(fluxall)
           areadif[j]=areadif[j]+totarea

           IF diag THEN BEGIN
              print, i, ctsall, ' counts will be removed '+$
                  'from annulus ',strtrim(string(j),1)+$
                     ' in area ', totarea
              IF useim THEN print, i, total(fluxall), ' flux will be removed'
           ENDIF
        ENDIF ELSE BEGIN
           npiec=jmax-jmin      ;number of pieces=jmax-jmin
           FOR k=jmax-1,jmin+1,-1 DO BEGIN
              ctindp=where((((x[ctind]^2.)+(y[ctind]^2.) LE (k*mrad/pixs)^2.)),cts)
              IF useim THEN BEGIN
                 IF cts EQ 0 THEN fluxp=0. ELSE fluxp=fluxall[ctindp]
              ENDIF
              IF (k LT numann) THEN BEGIN
                 psc[k]=psc[k]+(ctsall-cts)
                 IF useim THEN fluxout[k]=fluxout[k]+total(fluxall)-total(fluxp)
                 midarea=getmidarea(usepsf*pixs,k*mrad,rad[i]*pixs)
              ENDIF
              
               IF (k LT numann) THEN areadif[k]=areadif[k]+(totarea-midarea)  

               IF diag THEN  BEGIN
                  print, i, (ctsall-cts),' counts will be removed '+$
                      'from annulus ',strtrim(string(k),1)+$
                         ' in area ', (totarea-midarea)
                  IF useim THEN print, i, total(fluxall)-total(fluxp),$
                                       ' flux will be removed '
               ENDIF
               ctsall=cts
               fluxall=fluxp
               totarea=midarea
               radplk=k*mrad/pixs
               IF diag THEN   oplot,radplk*cos(angs),radplk*sin(angs)        
            ENDFOR ;for k
              psc[jmin]=psc[jmin]+cts
              IF useim THEN fluxout[jmin]=fluxout[jmin]+total(fluxp)
            areadif[jmin]=areadif[jmin]+midarea
            IF diag THEN BEGIN
               print, i, cts, ' counts will be removed '+$
                   'from annulus ',strtrim(string(jmin),1)+$
                      ' in area ', midarea
               IF useim THEN print, i, total(fluxp),$
                                    ' flux will be removed '
            ENDIF
            
         ENDELSE
        
         radpl2=jmax*mrad/pixs           
IF diag THEN oplot,radpl2*cos(angs),radpl2*sin(angs)
IF diag THEN stop

ENDIF ;ctsall

    ENDFOR
     
         
  ENDIF ELSE print,'No point source found'
  

END

  
function getmidarea, psfr, annumr, dist 
;area in closer circle

;r: psfr
;R: annum radius
;d:rad
     
area=(psfr^2.)*acos((dist^2.+psfr^2.-annumr^2.)/(2*dist*psfr)) + $
(annumr^2.)*acos((dist^2+annumr^2-psfr^2)/(2*dist*annumr)) - $
     0.5*sqrt((-dist+psfr+annumr)*(dist+psfr-annumr)*(dist-psfr+annumr)*(dist+psfr+annumr))

return, area

END
