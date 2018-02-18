; this is a wrapper for chandra profile for rgb equal count case

;fluxim='~/CHANDRA/4U1630_expoc/fluxall/1.5-6.0_flux.img'
fluxim1='~/CHANDRA/4U1630_backgr/rgb_back_c6/band1_back_c6_fluxed.img'
poss=[515,-173]                 ;use image and integer!
;poss=[160,477]
posback=[914.,1252.]
rback=50.
;rateb=0.075
num_an=120
corb=0
;nyr=[4e-3,5.]
;nxr=[0.,300.]
dopsf=0
;ps=0
elimps=0
mrad=5
diag=0
threshexp=0.15
uchip=6
ps=0
namef='psc.eps'

getprofiles=1

IF getprofiles THEN BEGIN

getprofile, fluxim1, rad_im, prof_im1c6flat, areas_im1c6flat, $
                ps=ps, fname=namef, $
                spos=poss, annum=num_an, radm=mrad, $
                corb=corb, backpos=posback, backr=rback, brate=rateb, $
                simrad=s_radius, simareas=s_areas, simprof=s_prof,$
                yrn=nyr, xrn=nxr, dopsf=dopsf, elimps=elimps,$
		cprof=profc_im1c6flat,$
                expthresh=threshexp, diag=diag,chip=uchip,$
                expo_prof=prof_expo_im1c6flat


nprof_im1c6flat=fltarr(2,num_an)
nprof_im1c6flat[0,*]=prof_im1c6flat[0,*]/areas_im1c6flat
nprof_im1c6flat[1,*]=prof_im1c6flat[1,*]/areas_im1c6flat
nprof_expo_im1c6flat=prof_expo_im1c6flat/areas_im1c6flat

fluxim2='~/CHANDRA/4U1630_backgr/rgb_back_c6/band2_back_c6_fluxed.img'

getprofile, fluxim2, rad_im, prof_im2c6flat, areas_im2c6flat, $
                ps=ps, fname=namef, $
                spos=poss, annum=num_an, radm=mrad, $
                corb=corb, backpos=posback, backr=rback, brate=rateb, $
                simrad=s_radius, simareas=s_areas, simprof=s_prof,$
                yrn=nyr, xrn=nxr, dopsf=dopsf, elimps=elimps,$
		cprof=profc_im2c6flat,$
                expthresh=threshexp, diag=diag,chip=uchip,$
                expo_prof=prof_expo_im2c6flat


nprof_im2c6flat=fltarr(2,num_an)
nprof_im2c6flat[0,*]=prof_im2c6flat[0,*]/areas_im2c6flat
nprof_im2c6flat[1,*]=prof_im2c6flat[1,*]/areas_im2c6flat
nprof_expo_im2c6flat=prof_expo_im2c6flat/areas_im2c6flat

fluxim3='~/CHANDRA/4U1630_backgr/rgb_back_c6/band3_back_c6_fluxed.img'

getprofile, fluxim3, rad_im, prof_im3c6flat, areas_im3c6flat, $
                ps=ps, fname=namef, $
                spos=poss, annum=num_an, radm=mrad, $
                corb=corb, backpos=posback, backr=rback, brate=rateb, $
                simrad=s_radius, simareas=s_areas, simprof=s_prof,$
                yrn=nyr, xrn=nxr, dopsf=dopsf, elimps=elimps,$
		cprof=profc_im3c6flat,$
                expthresh=threshexp, diag=diag,chip=uchip,$
                expo_prof=prof_expo_im3c6flat


nprof_im3c6flat=fltarr(2,num_an)
nprof_im3c6flat[0,*]=prof_im3c6flat[0,*]/areas_im3c6flat
nprof_im3c6flat[1,*]=prof_im3c6flat[1,*]/areas_im3c6flat
nprof_expo_im3c6flat=prof_expo_im3c6flat/areas_im3c6flat


ENDIF

savepar=1

IF savepar THEN BEGIN
   save, rad_im, prof_im1c6flat, prof_im2c6flat, prof_im3c6flat, $
   areas_im1c6flat, areas_im2c6flat, areas_im3c6flat, mrad, num_an, $
   profc_im1c6flat, profc_im2c6flat, profc_im3c6flat, $
   nprof_im1c6flat, nprof_im2c6flat, nprof_im3c6flat, $
   nprof_expo_im3c6flat, nprof_expo_im2c6flat, nprof_expo_im1c6flat, $
   filename='prof_rgbc6mrad5_flat.sav'
ENDIF
   

END