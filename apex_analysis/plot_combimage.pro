pro plot_combimage, inpfits, vrange, outim, ps=ps, fname=namef, $
                    useavg=useavg, maxscl=sclmax
  
; This program reads an image cube fits file and then plots a combined
; image with the given velocity range
;
; INPUTS
;
; inpfits: input fits file
; vrange: velocity range (fltarr(2)) to be used to merge data
;
; OUTPUTS
;
; outim: combined image
;
; OPTIONAL INPUTS
;
; ps: if set output a postscript file
; fname: if set, change the default output postscripto to namef
; useavg: if set calculates the average brightness, and uses fixed
;         scale (0-8) so that all maps scale the same. default: use sum of
;         brightness and scale to max.
; maxscl: If useavg is set, maximum scale can also be chose. Default:5
;
; USES
;
; IDL image routines and readfits
;
; USED BY
;
; NONE
;
; LOGS
;
; Created by E. Kalemci, Jan 2022
;
; Magic numbers exist, must allow size option
;

IF NOT keyword_set(ps) THEN ps=0
IF NOT keyword_set(namef) THEN namef='combimage.eps'
IF NOT keyword_set(useavg) THEN useavg=0
IF NOT keyword_set(sclmax) THEN sclmax=5.

IF ps THEN BEGIN
   set_plot, 'ps'
   device,/color
   loadct,3
   device,/encapsulated
   device, filename = namef
   device, yoffset = 2
   device, ysize = 16.
   device, xsize = 16.0
   !p.font=0
   device,/times
ENDIF ELSE BEGIN
   device,decomposed=0
   loadct,3
   window, 1, retain=2, xsize=280, ysize=280
ENDELSE


;Read parameters from the FITS file

imgs=readfits(inpfits, NanValue=0)
CDELT1  = 3600.*0.3787457038519E-02

;get a velocity array
v=-175.+findgen(500)*0.5

;merge images

indices=where((v GE vrange[0]) AND (v LE vrange[1]))
combimg=fltarr(28,28)


FOR i=0, n_elements(indices)-1 DO combimg=combimg+imgs[*,*,indices[i]]
IF useavg THEN combimg=combimg/n_elements(indices)

;first plot window to overplot annulus

!x.margin=[0.,0.]
!y.margin=[0.,0.]
plot,[0,0],[0,0],xr=[0,280],yr=[0,280],/xstyle,/ystyle

IF useavg THEN imgscl=bytscl(combimg,min=0, max=sclmax) ELSE imgscl=bytscl(combimg)
imgsclrbn=rebin(imgscl, 280, 280) ; BE CAREFUL, THIS PROVIDES SMOOTHING
tv,imgsclrbn

outim=imgsclrbn

;oplot annulus

angs=findgen(360)
r1=90./(cdelt1/10.)
r2=170./(cdelt1/10.)

oplot,140.+r1*cos(angs*2.*!PI/360),140.+r1*sin(angs*2.*!PI/360),color=255, thick=1.5
oplot,140.+r2*cos(angs*2.*!PI/360),140.+r2*sin(angs*2.*!PI/360),color=255, thick=1.5



IF ps THEN BEGIN
   device,/close
   set_plot,'x'
ENDIF

END
