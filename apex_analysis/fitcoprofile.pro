pro fitcoprofile, infile, b, l, outpar, $
                  ps=ps, fname=namef

;This program fits and plots CO J=1-0 map with multiple gausses saves
;important results like integrated velocities and near and far distances
;
; INPUTS
;
; infile: input FITS file, DHT36_Quad4_interp.fits for example
; l: galactic l
; b: galactic b
;
; OUTPUTS
;
; outpar: structure that holds information on gauss profiles,
;integrated profiles, and near and far distances
;
; OPTIONAL INPUTS
;
; ps: IF set postscript plot
; fname: base name for postscript plots
;
; USES
;
; IDL FITS routines
; mpcurvefit
; gaussfits
; 
; USED BY
;
; None
;
; Logs
;
; Created by EK Dec 2017
;

IF NOT keyword_set(ps) THEN ps=0
IF NOT keyword_set(namef) THEN namef='CO_gaussfits.eps'


;Read parameters from the FITS file

  dat=readfits(infile) ;datacube
  hdr=headfits(infile)
  v0=sxpar(hdr,'CRVAL1') ; v: velocity
  Dv=sxpar(hdr,'CDELT1')
  vROWS=sxpar(hdr,'NAXIS1')
  l0=sxpar(hdr,'CRVAL2') ; l: Galactic longitude
  Dl=-sxpar(hdr,'CDELT2') ;decreasing
  REFb=sxpar(hdr,'CRPIX3') ; b: galactic latitude
  Db=sxpar(hdr,'CDELT3')
  b0=-(Refb-1)*Db

;get a velocity array
  v=v0+findgen(vROWS)*Dv

;I think these are obsolete
;4U1630 336.91, 0.250

lval4u=floor(-((336.91)-l0)/Dl+0.5)
bval4u=floor(((0.25)-b0)/Db+0.5)

;SGR:336.98,-0.110

lvalsgr=floor(-((336.98)-l0)/Dl+0.5)
bvalsgr=floor(((-0.11)-b0)/Db+0.5)

IF ps THEN BEGIN
   set_plot, 'ps'
   device,/color
   loadct,5
   device,/encapsulated
   device, filename = namef
   device, yoffset = 2
   device, ysize = 16.
   device, xsize = 25.0
   !p.font=0
   device,/times
ENDIF ELSE BEGIN
   device,decomposed=0
   loadct,0
   window, 2, retain=2, xsize=600, ysize=400
ENDELSE

bval=floor((b-b0)/Db+0.5)
lval=floor(-(l-l0)/Dl+0.5)

cs=1.5

IF NOT ps THEN plot,v,dat[*,lval,bval],xr=[-150,0],psym=10, yr=[-0.2,7.8],$
     xtitle='v!Dlsr!N (km s!E-1!N)', ytitle='T!Dmb!N',charsize=cs,/ystyle,$
     color=255,/xstyle

;Make negative temperatures 0

tofit=dat[*,lval,bval]
xx=where(tofit LT 0.)
tofit[xx]=0.

;create model

;the program I use take FWHM as inputs, but the plotting uses sigma

 ; region 1

a11 = [.4,-113.,4.]
a12 = [1.5,-105.,1.4]
a13 = [1.4,-98.,1.6]


 ; region 2 (main region)

a21 = [4.6,-80.,4.]
a22 = [2.,-73,3.]
a23 = [1.1,-66,2.]
a24 = [.9,-56,1.8]

                                ; region 3

a31 = [.8,-47.,.8]
a32 = [1.8,-39.5,1.5]
a33 = [1.1,-34.,1.1]

;region 4

a41 = [1.6,-19.,1.2]
   
a = [a11,a12,a13,a21,a22,a23,a24,a31,a32,a33,a41]
afit=a
siginds=indgen(11)*3+2 ;gaussian uses sigma, not width
afit[siginds]=afit[siginds];*2.355

;a=[a11,a12,a13, a21,a22,a23,a24]

ngauss=11

IF NOT ps THEN BEGIN
   datai=dblarr(ngauss+1,n_elements(v))

   FOR i=0,ngauss-1 DO BEGIN
      datai[i,*] = gaussian(v,a[i*3:i*3+2])
      datai[ngauss,*]=datai[ngauss,*]+datai[i,*]
      oplot,v,datai[i,*],line=2
   ENDFOR
   oplot, v, datai[ngauss,*]
ENDIF

 ;   stop

    ; specify 4 regions
nregions = 4
wregions4=where((v GE -25.) AND (v LE -10.))
wregions3=where((v GE -50.) AND (v LT -25.))
wregions2=where((v GE -95.) AND (v LT -50.))
wregions1=where((v GE -125.) AND (v LT -95.))

regions=[[min(wregions1),max(wregions1)],[min(wregions2),max(wregions2)],$
            [min(wregions3),max(wregions3)],[min(wregions4),max(wregions4)]]

;    regions = [[-125,-95],[-95,-50],[-50,-25],[-25,-10]]
;inits = [[a11],[a12],[a13],[a21],[a22],[a23],[a24],[a31],[a32],[a33],[a41]]

inits = [[afit[0:2]],[afit[3:5]],[afit[6:8]],[afit[9:11]],[afit[12:14]],$
         [afit[15:17]],[afit[18:20]],[afit[21:23]],[afit[24:26]],$
         [afit[27:29]],[afit[30:32]]]

max_iters = 500
    
p = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0,0]}, 33) 
                       ; 33 = 11 gauss * 3 parameter per guass

;limit some unphysical parameters
    p[*].value = afit
    p[2].limited=[0,1]
    p[2].limits=[0.,9.]
    p[23].limited=[0,1]
    p[23].limits=[0.,4.]
    ;p[20].limited=[0,1]
    ;p[20].limits=[0.,8.]
    
    ; hold the first gaussians height fixed
  ;  p[0].fixed = 1

                                ; find all the fits at once
 
   yfit = gauss_fits(v,tofit,nregions,regions,inits,ngauss,max_iters,coefficients,errors,parinfo=p,quiet=1)

;    print,coefficients
;    print,errors
    ; unwrap the results and plot them

plot,v,dat[*,lval,bval],xr=[-130,-10],psym=10, yr=[-0.2,5.5],$
     xtitle='v!Dlsr!N (km s!E-1!N)', ytitle='T!Dmb!N',charsize=cs,/ystyle,$
     color=0,/xstyle


;prepare the output file
;cname: cloud name
;gpar: gauss fit parameters
;gparerr: gauss fit errors
;gint: integrated profile and its error
;vlsr: to be able to plot from the output
;l and b galactic l and b values

outpar=create_struct('cname',strarr(ngauss),'gpar',dblarr(3,ngauss),$
                     'gparerr',dblarr(3,ngauss),'gint',dblarr(2,ngauss),$
                      'Tmb',dat[*,lval,bval],'vlsr',v,'l',l,'b',b,$
                     'neardist',dblarr(2, ngauss),'fardist',dblarr(2,ngauss))


;plot and populate outpar
    
an=coefficients ; fit results
siginds=indgen(11)*3+2 ;gaussian uses sigma, not width
an[siginds]=an[siginds]/2.355 ;convert to sigma
errors[siginds]=errors[siginds]/2.355
dataf=dblarr(ngauss+1,n_elements(v))

;for far and near distances
lu=(360.-l)*!PI/180.  ; in radians
Ro=8.5D
Vo=220D
tp=Ro*cos(lu)
Vc=Vo
;Vr=Ro*sin(lu) * ((Vc/R) - (Vo/Ro))
;(Vr/(Ro*sin(lu))) + (Vo/Ro) = Vc/R
;R=Vc/((Vr/(Ro*sin(lu))) + (Vo/Ro))
;dR=sqrt(R^2.-(Ro*sin(lu))^2.)


FOR i=0,ngauss-1 DO BEGIN
    dataf[i,*] = gaussian(v,an[i*3:i*3+2])
    dataf[ngauss,*]=dataf[ngauss,*]+dataf[i,*]
    oplot,v,dataf[i,*],line=2
    txtp='MC'+strtrim(string(floor(an[i*3+1]+0.5)),1)
    xyouts,an[i*3+1]-3,an[i*3]+0.2,txtp,charsize=0.9
    outpar.cname[i]=txtp
    outpar.gpar[*,i]=an[i*3:i*3+2]
    outpar.gparerr[*,i]=errors[i*3:i*3+2]
    outpar.gint[0,i]=an[i*3]*sqrt(!PI)*sqrt(2)*an[i*3+2]
    errf=((errors[i*3]/an[i*3])+(errors[i*3+2]/an[i*3+2]))/sqrt(2.)
    outpar.gint[1,i]=outpar.gint[0,i]*errf
    R=Vc/((abs(an[i*3+1]))/(Ro*sin(lu)) + (Vo/Ro))
    dR=sqrt(R^2.-(Ro*sin(lu))^2.)
    outpar.neardist[0,i]=tp-dR
    outpar.fardist[0,i]=tp+dR
ENDFOR
oplot, v, dataf[ngauss,*]


IF ps THEN BEGIN
   device,/close
   set_plot,'x'
ENDIF



END
