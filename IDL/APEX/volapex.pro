pro volapex, infile, regimagesws, neard=neard, fard=fard, minmax=mimax, drange=pranged

;This program reads CO map file as XYZ, and using a galactic model
;converts velocities to depths and plots a volume
;
; INPUTS
;
; infile: CO MAP
;
; OUTPUTS
;
; 3D volume plot of given fard or neard
; regimagews: regular gridded image
;
; OPTIONAL INPUTS
;
; neard: If set assume all clouds and the source are in near distance
; fard: If set assume all clouds and the source in the far distance
; minmax: If set, use the given scale in determining pixel colors
; drange: If set, use the given distange range
;
; USES
;
; volume
;
; USED BY
;
; NONE
;
; LOGS
;
;created by EK, May 2023
;currently only uses negative velocities, positive velocities are not
;of interest
;

IF NOT keyword_set(neard) THEN neard=0
IF NOT keyword_set(fard) THEN fard=0
IF NOT keyword_set(mimax) THEN mimax=[0.5,15.]
IF NOT keyword_set(pranged) THEN BEGIN   ;this is only for plotting purposes
   pranged=[1.,15.]
   IF fard THEN pranged=[9.,15.]
   IF neard THEN pranged=[1.,7.]
ENDIF


;basic info

;source galactic coordinates
b=0.250
l=336.911

;pixel size
pixsz = 13.634845D


;Read the fits file
images=readfits(infile, hdr)

nimages=sxpar(hdr, 'NAXIS3')
naxis1=sxpar(hdr, 'NAXIS1')
naxis2=sxpar(hdr, 'NAXIS2')


;find the velocity scale properly

CRVAL3  = sxpar(hdr,'CRVAL3')
CDELT3  = sxpar(hdr,'CDELT3')
CRPIX3  = sxpar(hdr,'CRPIX3')

vstart=CRVAL3-(CDELT3*CRPIX3)
vstart=vstart/1000.             ;turn into km/s
delv=cdelt3/1000.

vel=findgen(nimages)*delv+vstart

;Calculate distances for given velocities CHECK THIS

lu=(360.-l)*!PI/180.  ; in radians
Ro=8.5D
Vo=220D
tp=Ro*cos(lu)
Vc=Vo
R=Vc/((abs(vel))/(Ro*sin(lu)) + (Vo/Ro))
dR=sqrt(R^2.-(Ro*sin(lu))^2.)
neardist=tp-dR
fardist=tp+dR

;Not sure about positive velocities, drop them for now

useneg=where(vel LT -5)

;volume likes regularly spaced grids
;interpolate for all available pixels
;some pixels are missing, turn them into 0

nregd=600 ;each slice is 10pc 6 kpc for both far and near distances
regimagesn=fltarr(NAXIS1,NAXIS2,nregd)
regimagesf=fltarr(NAXIS1,NAXIS2,nregd)


;far the entire range one needs to calculate cases for far and near
;distances separately

;for far distances

FOR i=0, naxis1-1 DO BEGIN
   FOR j=0, naxis2-1 DO BEGIN
      f=transpose(images[i,j,*])
      xx=where(finite(f, /NAN))
      f[xx]=0.
      newf=interpol(f(useneg),fardist(useneg),9.+findgen(600)/100.,/NAN,/spline)
      regimagesf[i,j,*]=newf
   ENDFOR
ENDFOR

;for near distances, remember the distances decreases with index now


FOR i=0, naxis1-1 DO BEGIN
   FOR j=0, naxis2-1 DO BEGIN
      f=transpose(images[i,j,*])
      xx=where(finite(f, /NAN))
      f[xx]=0.
      newf=interpol(f(useneg),neardist(useneg),7.-findgen(600)/100.,/NAN,/spline)
      regimagesn[i,j,*]=newf
   ENDFOR
ENDFOR

;add point source as a cube?

xs=24
ys=26
fards=[11.5, 0.3] ;apj results
neards=[4.7,0.3] ;apj results

fdsi=floor(((fards[0]-9.)-[fards[1]/2.,-fards[1]/2.])*100)
ndsi=floor(((7.-neards[0])-[neards[1]/2.,-neards[1]/2.])*100)


regimageswsn=regimagesn
regimageswsf=regimagesf

regimageswsf[xs-1:xs+1,ys-1:ys+1,fdsi[0]:fdsi[1]]=15.
regimageswsn[xs-1:xs+1,ys-1:ys+1,ndsi[0]:ndsi[1]]=15.

;volume location and dimensions
;let's try to do this properly

;crval1=sxpar(hdr, 'CRVAL1')

;magic number method FIX THIS BEFORE PAPER
imgxrange=[248.64055,248.3720]
imgyrange=[-47.49124,-47.2942946]
imgzrange=pranged


vloc=[imgxrange[0],imgyrange[0],imgzrange[0]]
vdim=[(imgxrange[1]-imgxrange[0]),(imgyrange[1]-imgyrange[0]),(imgzrange[1]-imgzrange[0])]

IF neard THEN BEGIN
   ;we need a cut on give pranged
   ndspi=floor([(7.-pranged[1]),(7.-pranged[0])])*100
   IF ndspi[0] LT 0 THEN ndspi[0]=0
   IF ndspi[1] GE 600 THEN ndspi[1]=599
   regimagesws=regimageswsn[*,*,ndspi[0]:ndspi[1]]
   vloc=[imgxrange[0],imgyrange[0],imgzrange[1]]
   vdim=[(imgxrange[1]-imgxrange[0]),(imgyrange[1]-imgyrange[0]),(imgzrange[0]-imgzrange[1])]
ENDIF

IF fard THEN BEGIN
   fdspi=floor([(pranged[0]-9.),(pranged[1]-9.)])*100
   IF fdspi[0] LT 0 THEN fdspi[0]=0
   IF fdspi[1] GE 600 THEN fdspi[1]=599
   regimagesws=regimageswsf[*,*,fdspi[0]:fdspi[1]]
ENDIF



v = VOLUME(regimagesws, RENDER_EXTENTS=0, $
   HINTS = 3, /AUTO_RENDER, $
   RGB_TABLE0=15, AXIS_STYLE=2, $
   RENDER_QUALITY=2, BACKGROUND_COLOR='gray', $
   DEPTH_CUE=[0, 2], /PERSPECTIVE, $
   VOLUME_LOCATION=vloc, $
   VOLUME_DIMENSIONS=vdim, min_value=mimax[0],max_value=mimax[1])


END
