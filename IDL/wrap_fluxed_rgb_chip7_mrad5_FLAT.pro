; this is a wrapper for chandra profile for rgb equal count case

;fluxim='~/CHANDRA/4U1630_expoc/fluxall/1.5-6.0_flux.img'
fluxim1='~/CHANDRA/4U1630_backgr/rgb_back_c7/band1_back_c7_fluxed.img'
poss=[576,868]                 ;use image and integer!
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
uchip=7
ps=0
namef='psc.eps'

getprofiles=1

IF getprofiles THEN BEGIN

getprofile, fluxim1, rad_im, prof_im1c7flat, areas_im1c7flat, $
                ps=ps, fname=namef, $
                spos=poss, annum=num_an, radm=mrad, $
                corb=corb, backpos=posback, backr=rback, brate=rateb, $
                simrad=s_radius, simareas=s_areas, simprof=s_prof,$
                yrn=nyr, xrn=nxr, dopsf=dopsf, elimps=elimps,$
		cprof=profc_im1c7flat,$
                expthresh=threshexp, diag=diag,chip=uchip,$
                expo_prof=prof_expo_im1c7flat


nprof_im1c7flat=fltarr(2,num_an)
nprof_im1c7flat[0,*]=prof_im1c7flat[0,*]/areas_im1c7flat
nprof_im1c7flat[1,*]=prof_im1c7flat[1,*]/areas_im1c7flat
nprof_expo_im1c7flat=prof_expo_im1c7flat/areas_im1c7flat

fluxim2='~/CHANDRA/4U1630_backgr/rgb_back_c7/band2_back_c7_fluxed.img'

getprofile, fluxim2, rad_im, prof_im2c7flat, areas_im2c7flat, $
                ps=ps, fname=namef, $
                spos=poss, annum=num_an, radm=mrad, $
                corb=corb, backpos=posback, backr=rback, brate=rateb, $
                simrad=s_radius, simareas=s_areas, simprof=s_prof,$
                yrn=nyr, xrn=nxr, dopsf=dopsf, elimps=elimps,$
		cprof=profc_im2c7flat,$
                expthresh=threshexp, diag=diag,chip=uchip,$
                expo_prof=prof_expo_im2c7flat


nprof_im2c7flat=fltarr(2,num_an)
nprof_im2c7flat[0,*]=prof_im2c7flat[0,*]/areas_im2c7flat
nprof_im2c7flat[1,*]=prof_im2c7flat[1,*]/areas_im2c7flat
nprof_expo_im2c7flat=prof_expo_im2c7flat/areas_im2c7flat

fluxim3='~/CHANDRA/4U1630_backgr/rgb_back_c7/band3_back_c7_fluxed.img'

getprofile, fluxim3, rad_im, prof_im3c7flat, areas_im3c7flat, $
                ps=ps, fname=namef, $
                spos=poss, annum=num_an, radm=mrad, $
                corb=corb, backpos=posback, backr=rback, brate=rateb, $
                simrad=s_radius, simareas=s_areas, simprof=s_prof,$
                yrn=nyr, xrn=nxr, dopsf=dopsf, elimps=elimps,$
		cprof=profc_im3c7flat,$
                expthresh=threshexp, diag=diag,chip=uchip,$
                expo_prof=prof_expo_im3c7flat


nprof_im3c7flat=fltarr(2,num_an)
nprof_im3c7flat[0,*]=prof_im3c7flat[0,*]/areas_im3c7flat
nprof_im3c7flat[1,*]=prof_im3c7flat[1,*]/areas_im3c7flat
nprof_expo_im3c7flat=prof_expo_im3c7flat/areas_im3c7flat


ENDIF

savepar=1

IF savepar THEN BEGIN
   save, rad_im, prof_im1c7flat, prof_im2c7flat, prof_im3c7flat, $
   areas_im1c7flat, areas_im2c7flat, areas_im3c7flat, mrad, num_an, $
   profc_im1c7flat, profc_im2c7flat, profc_im3c7flat, $
   nprof_im1c7flat, nprof_im2c7flat, nprof_im3c7flat, $
   nprof_expo_im3c7flat, nprof_expo_im2c7flat, nprof_expo_im1c7flat, $
   filename='prof_rgbc7mrad5_flat.sav'
ENDIF
   

END