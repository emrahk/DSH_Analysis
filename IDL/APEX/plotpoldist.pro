pro plotpoldist, infile, vrange, poldist, ps=ps, fname=namef, chkcode=chkcode
  
;This program reads FITS images and already calculated polar
;distribution of brightess, and for the given velocity range plots the
;image, overplots the polar regions as well as the distribution.
;
;
; INPUTS
;
; infile: FITS file with the apex images
; vrange: range of velocities
; poldist: already calculated brightness values
;
; OPTIONAL INPUTS
;
; ps: IF set, postscript output
; fname: IF set, use namef for the postscript file
; chkcode: IF set allows checking some parts of the code
;
; OUTPUTS
;
; image
;
; USES
;
; idlastro library
; plotimage
;
; USED BY
;
; APEX and CHANDRA analysis polar distribution codes
;
; CREATED by Emrah Kalemci, Jan 2023
;
; 7 Feb 2023 Fixed radial velocity indexing
;

IF NOT keyword_set(ps) THEN ps=0
IF NOT keyword_set(namef) THEN namef='imganddist.eps' 
IF NOT keyword_set(chkcode) THEN chkcode=0

loadct,5
if ps then begin
   set_plot, 'ps'
   device,/color
;   loadct,5
   device,/encapsulated
   device, filename = namef
   device, yoffset = 2
   device, ysize = 20.
   device, xsize = 20.
   !p.font=0
   device,/times
   backc=0
   axc=0
endif

IF NOT ps THEN BEGIN
   window, 1, retain=2, xsize=600, ysize=600
   backc=0
   axc=255
   device,decomposed=0
ENDIF

cs=1.3

;multiplot, [1,3], mxtitle='Time (MJD-57000 days)', mxtitsize=1.2
;start with the image

images=readfits(infile, hdr)

;position of the source

pos=[248.5067, -47.393]

;find the velocity scale properly

CRVAL3  = sxpar(hdr,'CRVAL3')                                     
CDELT3  = sxpar(hdr,'CDELT3')
CRPIX3  = sxpar(hdr,'CRPIX3')

vstart=CRVAL3-(CDELT3*CRPIX3)
vstart=vstart/1000.             ;turn into km/s
delv=cdelt3/1000.

vindex=lonarr(2)
vindex[0]=floor((vrange[0]-vstart)/delv)-1
vindex[1]=floor((vrange[1]-vstart)/delv)-1 ;think about these ranges again

;create composite image and surfacec brightness by adding images in the range

sz=size(images)

img2pl=fltarr(sz[1],sz[2])

img2pl=images(*,*,vindex[0])
sbr2plt=poldist.sbr[vindex[0],*]*poldist.areas

nim=vindex[1]-vindex[0]+1

IF vindex[0] LT vindex[1] THEN BEGIN
   FOR k=vindex[0]+1, vindex[1] DO BEGIN
      img2pl=img2pl+images(*,*,k)
      sbr2plt=sbr2plt+poldist.sbr[k,*]*poldist.areas
   ENDFOR
ENDIF

sbr2pl=sbr2plt/poldist.areas
sbr2ple=sqrt(sbr2plt)/poldist.areas

; plot the image ahh magic numbers fix later, this problem can be
; solved using the header later, there could also be a shift

imrange=[min(img2pl,/nan),max(img2pl,/nan)]
plotimage, img2pl,imgxrange=[248.64055,248.3720],$
           imgyrange=[-47.49124,-47.2942946],range=imrange,$
           background=backc,axiscolor=axc

;oplot the source

plotsym,0,/fill

oplot, [pos[0],pos[0]],[pos[1],pos[1]], psym=8

;plot wedges

noa=poldist.noa
limrad=poldist.radlim

;this is to check the procedure, can be turned off

If chkcode THEN BEGIN
   FOR i=0, noa-1 DO BEGIN
      anglr=360./noa
      angles=[i*anglr,(i+1)*anglr]
      xx=where((poldist.mapr GE limrad[0]) AND (poldist.mapr LT limrad[1]) AND $
            (poldist.mapth GE angles[0]) AND (poldist.mapth LT angles[1]))
      xyad,hdr,(xx mod sz[1]), (xx/sz[1]), ra, dec ;CHECK THIS LATER
      oplot,ra,dec,psym=8, color=15+(i mod 2)*220
   ENDFOR
ENDIF

;oplot limiting radii and pies - learn later

;angvar=findgen(60)*6.
;oplot,pos[0]+limrad[0]*cos(angvar*!PI/180.)/3600.,$
;      pos[1]+limrad[0]*sin(angvar*!PI/180.)/3600.,color=255,thick=3
;
;oplot,pos[0]+limrad[1]*cos(angvar*!PI/180.)/3600.,$
;      pos[1]+limrad[1]*sin(angvar*!PI/180.)/3600.,color=255,thick=3

;plot the surface brightness in another plot

IF NOT ps THEN window, 2, retain=2, xsize=600, ysize=300

angvar=findgen(noa)*360./noa
ploterror, angvar, sbr2pl, sbr2ple,psym=10, xrange=[0.,360.],/xstyle,$
      xtitle='Polar angle',ytitle='Brightness'


END



