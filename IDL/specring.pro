pro specring, infile, specv, vel, radlim=limrad, ps=ps, fname=namef, allspec=specvall
  
;This program reads a FITS header and images, and find the spectrum in
;the given radial limits
;
;
; INPUTS
;
; infile: FITS file with the apex images
;
; OUTPUTS
;
; specv: spectral array
; vel: velocities for spectral array
;
; OPTIONAL INPUTS
;
; radlim: limiting radii for the ring, default: 100-200 ''
; ps: IF set provide postscript figure
; fname: IF set use this name for the postscript output
;
; OPTIONAL OUTPUTS
;
; allspec: spectrum for the entire image
;
; USES
;
; idlastro library
; mappolcoord
;
; USED BY
;
; APEX and CHANDRA analysis polar distribution codes
;
; CREATED by Emrah Kalemci, Jan 2023
;
;

  IF NOT keyword_set(limrad) THEN limrad=[100., 200.]
  IF NOT keyword_set(ps) THEN ps=0
  IF NOT keyword_set(namef) THEN namef='specring.eps'

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
   window, 1, retain=2, xsize=600, ysize=330
   backc=0
   axc=255
   device,decomposed=0
ENDIF
   
;Read the fits file

images=readfits(infile, hdr)

;map the pixel coordinates to polar coordinates for the given central
;pos

pos=[248.5067, -47.393]
mappolcoord, hdr, pos, mapr, mapth

;number of velocities/frames

nvel = sxpar(hdr,'NAXIS3')

imgarr=intarr(2)
imgarr[0]=sxpar(hdr, 'NAXIS1')
imgarr[1]=sxpar(hdr, 'NAXIS2')


specv=fltarr(nvel)    ;spectrum of the given ring
specall=fltarr(nvel) ;spectrum of the entire area

;find the velocity scale properly

CRVAL3  = sxpar(hdr,'CRVAL3')
CDELT3  = sxpar(hdr,'CDELT3')
CRPIX3  = sxpar(hdr,'CRPIX3')

vstart=CRVAL3-(CDELT3*CRPIX3)
vstart=vstart/1000.             ;turn into km/s
delv=cdelt3/1000.

vel=findgen(nvel)*delv+vstart

;integrate temperatures in the ring and overall

FOR i=0L, nvel-1L DO BEGIN
   specall[i]=total(images[*,*,i],/NAN)
   xx=where((mapr GE limrad[0]) AND (mapr LT limrad[1]))
   specv[i]=total(images[xx+(i*imgarr[0]*imgarr[1])],/NAN)
ENDFOR

plot, vel, specall, psym=10, background=backc, xtitle='Velocity (km/s)',ytitle='T (keV)',yrange=[-100, max(specall)*1.1],/ystyle, xrange=[-150.,20.],/xstyle

oplot, vel, specv, psym=10, color=150, line=2


IF ps THEN BEGIN
   device,/close
   set_plot,'x'
ENDIF



END
