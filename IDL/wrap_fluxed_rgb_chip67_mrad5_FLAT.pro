restore,'prof_rgbc6mrad5_deflare_flb.sav'
restore,'prof_rgbc7mrad5_deflare_flb.sav'
restore,'prof_rgbc7mrad5_flat.sav'
restore,'prof_rgbc6mrad5_flat.sav'

prof_im1c67=fltarr(2,num_an)
prof_im1c67[0,*]=prof_im1c7[0,*]+prof_im1c6[0,*] ;add profiles
errorfrac=sqrt(profc_im1c7[0,*]+profc_im1c6[0,*])/float(profc_im1c7[0,*]+profc_im1c6[0,*])
prof_im1c67[1,*]=prof_im1c67[0,*]*errorfrac
areas_im1c67=areas_im1c6+areas_im1c7

prof_im1c67flatb=fltarr(2,num_an)

prof_im1c6flatb=prof_im1c6
prof_im1c6flatb[0,*]=prof_im1c6[0,*]-prof_im1c6flat[0,*] ;areas may not be equal.
prof_im1c7flatb=prof_im1c7
prof_im1c7flatb[0,*]=prof_im1c7[0,*]-prof_im1c7flat[0,*]

;make a cut on chip 7 at 450? and 120 at chip6

yy=where(rad_im LT 120.)
zz=where(rad_im GE 400.)

prof_im1c67flatb[0,*]=prof_im1c7flatb[0,*]+prof_im1c6flatb[0,*]
prof_im1c67flatb[0,yy]=prof_im1c7flatb[0,yy]
prof_im1c67flatb[0,zz]=prof_im1c6flatb[0,zz]

prof_im1c67flatb[1,*]=prof_im1c67flatb[0,*]*errorfrac
prof_im1c67flatb[1,yy]=prof_im1c7flatb[0,yy]*errorfrac[yy]
prof_im1c67flatb[1,zz]=prof_im1c6flatb[0,zz]*errorfrac[zz]

areas_im1c67flatb=areas_im1c6+areas_im1c7
areas_im1c67flatb[yy]=areas_im1c7[yy]
areas_im1c67flatb[zz]=areas_im1c6[zz]

nprof_im1c67flatb=fltarr(2,num_an)
nprof_im1c67flatb[0,*]=prof_im1c67flatb[0,*]/areas_im1c67flatb
nprof_im1c67flatb[1,*]=prof_im1c67flatb[1,*]/areas_im1c67flatb

;------

prof_im2c67=fltarr(2,num_an)
prof_im2c67[0,*]=prof_im2c7[0,*]+prof_im2c6[0,*]
errorfrac=sqrt(profc_im2c7[0,*]+profc_im2c6[0,*])/float(profc_im2c7[0,*]+profc_im2c6[0,*])
prof_im2c67[1,*]=prof_im2c67[0,*]*errorfrac
areas_im2c67=areas_im2c6+areas_im2c7

prof_im2c67flatb=fltarr(2,num_an)

prof_im2c6flatb=prof_im2c6
prof_im2c6flatb[0,*]=prof_im2c6[0,*]-prof_im2c6flat[0,*]
prof_im2c7flatb=prof_im2c7
prof_im2c7flatb[0,*]=prof_im2c7[0,*]-prof_im2c7flat[0,*]


prof_im2c67flatb[0,*]=prof_im2c7flatb[0,*]+prof_im2c6flatb[0,*]
prof_im2c67flatb[0,yy]=prof_im2c7flatb[0,yy]
prof_im2c67flatb[0,zz]=prof_im2c6flatb[0,zz]

prof_im2c67flatb[1,*]=prof_im2c67flatb[0,*]*errorfrac
prof_im2c67flatb[1,yy]=prof_im2c7flatb[0,yy]*errorfrac[yy]
prof_im2c67flatb[1,zz]=prof_im2c6flatb[0,zz]*errorfrac[zz]


areas_im2c67flatb=areas_im2c6+areas_im2c7
areas_im2c67flatb[yy]=areas_im2c7[yy]
areas_im2c67flatb[zz]=areas_im2c6[zz]

nprof_im2c67flatb=fltarr(2,num_an)
nprof_im2c67flatb[0,*]=prof_im2c67flatb[0,*]/areas_im2c67flatb
nprof_im2c67flatb[1,*]=prof_im2c67flatb[1,*]/areas_im2c67flatb

;===

prof_im3c67=fltarr(2,num_an)
prof_im3c67[0,*]=prof_im3c7[0,*]+prof_im3c6[0,*]
errorfrac=sqrt(profc_im3c7[0,*]+profc_im3c6[0,*])/float(profc_im3c7[0,*]+profc_im3c6[0,*])
prof_im3c67[1,*]=prof_im3c67[0,*]*errorfrac
areas_im3c67=areas_im3c6+areas_im3c7

prof_im3c67flatb=fltarr(2,num_an)

prof_im3c6flatb=prof_im3c6
prof_im3c6flatb[0,*]=prof_im3c6[0,*]-prof_im3c6flat[0,*]
prof_im3c7flatb=prof_im3c7
prof_im3c7flatb[0,*]=prof_im3c7[0,*]-prof_im3c7flat[0,*]

prof_im3c67flatb[0,*]=prof_im3c7flatb[0,*]+prof_im3c6flatb[0,*]
prof_im3c67flatb[0,yy]=prof_im3c7flatb[0,yy]
prof_im3c67flatb[0,zz]=prof_im3c6flatb[0,zz]

prof_im3c67flatb[1,*]=prof_im3c67flatb[0,*]*errorfrac
prof_im3c67flatb[1,yy]=prof_im3c7flatb[0,yy]*errorfrac[yy]
prof_im3c67flatb[1,zz]=prof_im3c6flatb[0,zz]*errorfrac[zz]

areas_im3c67flatb=areas_im3c6+areas_im3c7
areas_im3c67flatb[yy]=areas_im3c7[yy]
areas_im3c67flatb[zz]=areas_im3c6[zz]

nprof_im3c67flatb=fltarr(2,num_an)
nprof_im3c67flatb[0,*]=prof_im3c67flatb[0,*]/areas_im3c67flatb
nprof_im3c67flatb[1,*]=prof_im3c67flatb[1,*]/areas_im3c67flatb


doplot=1

IF doplot THEN BEGIN

ps=1

IF ps THEN BEGIN
   set_plot,'ps'
;   device,/color
;   loadct,5
   device,/encapsulated
   device, filename = 'rgb_profc67mrad5_deflare_FLAT.eps'
   device, yoffset = 2
   device, ysize = 18.
   ;device, xsize = 12.0
   ;!p.font=0
   ;device,/times
ENDIF

multiplot,/default

multiplot, [1,4], mxtitle='Radius (arcsec)', $
           mytitle='10!E-9!N Photons cm!E-2!N s!E-1!N asec!E-2!N',$
           mxtitsize=1.2, mytitsize=1.2

;nyr=[2e-10,4.2e-9]*1E9

nyr=[0,4.2e-9]*1e9

ploterror, rad_im, nprof_im1c67flatb[0,*]*1E9,nprof_im1c67flatb[1,*]*1E9,$
           yr=nyr,/nohat,/ystyle, chars=1.;,/ylog

xyouts,!x.crange[1]*0.7,!y.crange[1]*0.8,'1.5-2.25 keV'

;oplot,[80.,80.],nyr,color=0
;oplot,[260.,260.],nyr,color=0,line=1
;oplot,[140.,140.],nyr,color=0,line=2


multiplot

ploterror, rad_im, nprof_im2c67flatb[0,*]*1E9,nprof_im2c67flatb[1,*]*1E9,$
           yr=nyr,/nohat,/ystyle, chars=1.;,/ylog

xyouts,!x.crange[1]*0.7,!y.crange[1]*0.8,'2.25-3.15 keV'

;oplot,[80.,80.],nyr,color=0
;oplot,[260.,260.],nyr,color=0,line=1
;oplot,[140.,140.],nyr,color=0,line=2


multiplot

ploterror, rad_im, nprof_im3c67flatb[0,*]*1E9, nprof_im3c67flatb[1,*]*1E9,$
           yr=nyr,/nohat,/ystyle, chars=1.;,/ylog

xyouts,!x.crange[1]*0.7,!y.crange[1]*0.8,'3.15-5.0 keV'

;oplot,[80.,80.],nyr,color=0
;oplot,[260.,260.],nyr,color=0,line=1
;oplot,[140.,140.],nyr,color=0,line=2

multiplot

nprof_totc67flatb=fltarr(2,num_an)
nprof_totc67flatb[0,*]=nprof_im1c67flatb[0,*]+nprof_im2c67flatb[0,*]+nprof_im3c67flatb[0,*]

;get error
;yy=where(rad_im LT 120.)
;zz=where(rad_im GE 400.)

totctsc67=transpose(profc_im1c7[0,*]+profc_im2c7[0,*]+profc_im3c7[0,*])+$
	transpose(profc_im1c6[0,*]+profc_im2c6[0,*]+profc_im3c6[0,*])
totctsc67[yy]=totctsc67[yy]-transpose(profc_im1c6[0,yy]+profc_im2c6[0,yy]+profc_im3c6[0,yy])
totctsc67[zz]=totctsc67[zz]-transpose(profc_im1c7[0,zz]+profc_im2c7[0,zz]+profc_im3c7[0,zz])

errfracc67=sqrt(totctsc67)/totctsc67
nprof_totc67flatb[1,*]=nprof_totc67flatb[0,*]*errfracc67

;nyrt=[8e-10,9e-9]*1E9
nyrt=[0,10.5e-9]*1e9
ploterror, rad_im, nprof_totc67flatb[0,*]*1E9, nprof_totc67flatb[1,*]*1E9,$
           yr=nyrt,/nohat,/ystyle, chars=1.

xyouts,!x.crange[1]*0.7,!y.crange[1]*0.8,'1.5-5.0 keV'

;oplot,[80.,80.],nyrt,color=0
;oplot,[260.,260.],nyrt,color=0,line=1
;oplot,[140.,140.],nyrt,color=0,line=2

multiplot,/default

IF ps THEN BEGIN
   device,/close
   set_plot,'x'
ENDIF

ENDIF


savepar=0

IF savepar THEN save, rad_im, prof_im1c67, prof_im2c67, prof_im3c67, $
   areas_im1c67, areas_im2c67, areas_im3c67, mrad, num_an, $
   prof_im1c67flatb, prof_im2c67flatb, prof_im3c67flatb, $
   areas_im1c67flatb, areas_im2c67flatb, areas_im3c67flatb, $
   nprof_im1c67flatb, nprof_im2c67flatb, nprof_im3c67flatb, nprof_totc67flatb, $
   filename='prof_rgbc67mrad5_deflare_FLAT.sav'

;fitline

fitline=0

IF fitline THEN BEGIN

;multiplot,/default
nyrt=[0.15,12.]
ploterror, rad_im, nprof_totc67flb[0,*]*1E9, nprof_totc67flb[1,*]*1E9,$
            yr=nyrt,/nohat,/ystyle, chars=1.3,/xlog,/ylog,xr=[60.,500.],$
	    xtitle='Radius (arcsec)', $
            ytitle='10!E-9!N Photons cm!E-2!N s!E-1!N asec!E-2!N'

radstart=135.
radend=300.

xx=where((rad_im GT radstart) AND (rad_im LT radend) AND (nprof_totc67flb[0,*] GT 0.))
nprofit=nprof_totc67flb[0,xx]
newrad=rad_im[xx]

lognprof=transpose(alog(nprofit))
logCi=12.
lograd=alog(newrad)
ai=-3.
errs=transpose(nprof_totc67flb[1,xx]/nprof_totc67flb[0,xx])

;plot, lograd,lognprof
;yy=findgen(60.)/10.+1.
;oplot,yy,logci+ai*yy

res=linfit(lograd,lognprof,chi=nchi,/double,measure_errors=errs,sigma=nsig,yfit=fitres)

oplot,exp(lograd),exp(fitres)*1E9,color=0

nrads=findgen(470)+130
oplot,nrads,exp(res[0])*nrads^res[1]*1E9,color=0

oplot,[60.,500.],[1.65,1.65],color=0
oplot,[275,275],10^(!y.crange),color=0

ENDIF


END

