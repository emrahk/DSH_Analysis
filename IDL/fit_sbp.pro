;PRO FIT_SBP

restore,'DATA/resparex_reid_use1.sav' ; cloud parameters
restore,'DATA/parcloud_fit.sav'
restore,'DATA/prof_rgbc67mrad5_deflare_FLAT.sav'

COMMON myVars, cl_dist, weight, thcl, time_hires, Foft_hires, tobs

  
;IF ps THEN BEGIN
;   set_plot,'ps'
;   device,/color
;   loadct,5
;   device,/encapsulated
;   device, filename = namef
;   device, yoffset = 2
;   device, ysize = 18.
;   ;device, xsize = 12.0
;   !p.font=0
;   device,/times
;ENDIF  

LOADCT, 39, /SILENT
DEVICE, decomposed = 0

!p.multi=[0,1,2]

;PLOT MAXI LIGHT CURVE 2-4 keV

sed_data=read_csv('DATA/glcbin24.0h_regbg_hv0.csv')
mjdm=sed_data.field1
m24=sed_data.field4
e24=sed_data.field5
;m410=sed_data.field6
;e410=sed_data.field7

chand=57789.4 ;date of chandra (or earlier observations)

ploterror,mjdm-50000.,m24,e24,psym=4,/nohat,xr=[7600,7800],$
          xtitle='Time (MJD-50000 days)',$
          ytitle='ph cm!E-2!N s!E-1!N',charsize=1.3,/ylog,yr=[1e-5,1.]


;swift data

sdates=[57769.845,57771.775,57773.237,57777.35]
sflux24=[0.00110, 0.000471,0.000305,0.000076]
sflux15=[0.00165,0.000772,0.000494,0.000124]

plotsym,0,/fill
oplot,sdates-50000.,sflux24,psym=8,color=250

;oploterror,mjdm,m410,e410,psym=5,/nohat
oplot,[chand,chand]-50000.,10^(!y.crange),color=250,line=2,thick=2

;arrow, chand-50000., 0.3, chand-50000., 1E-5, /data, color=0

t0=57625.  ;start of outburst
tl=57761.  ;end of dates that we use maxi data

delt=tl-t0
;oplot,[t0,t0]-50000.,!y.crange,color=0    ;show if necessary
;oplot,[tl,tl]+9-50000.,!y.crange,color=0  ;show if necessary

xx=where((mjdm GT t0) AND (mjdm LT tl))

time=indgen(floor(tl-t0)) ;get maxi part

m24[2157:2159]=0.285  ; average out some stupid data
m24[2182:2183]=0.18
F=interpol(m24[xx],mjdm[xx],time+t0) ; interpolate maxi data per day
;oplot,time+t0,F,color=0
;from MJD 57755 to 57780 add exponential decay

edecay=0.34 ;too high

dd=findgen(28)+1 ; decay almost all the way to Chandra data (2 days before)
Foft_add=F[n_elements(F)-1L]*exp(-dd*edecay) ; 
Foft=[F,Foft_add]

;Foft=[F,0.05,0.045,0.04,0.035,0.03,0.025,0.02,0.015,0.01,0.005]
oplot,t0+indgen(floor(tl-t0)+28)-50000.,Foft,color=250

;====== END OF Plotting F(t) ===========


;====== Fitting the SBP(theta) ===========
errsys = 0.04

data_x = rad_im
data_y = NPROF_IM2C67FLATB[0,*]
data_e = NPROF_IM2C67FLATB[1,*]+NPROF_IM2C67FLATB[0,*]*errsys

ploterror,data_x, data_y, data_e,$
          yr=[1e-10,0.7e-8], xr=[0,600],/nohat, xtitle='Radius (arcsec)', $
          ytitle='ph cm!E-2!N s!E-1!N asec!E-2!N',charsize=1.3,ylog=logy

par_name = ['Source Distance', 'Normalization']
units = ['kpc', '']
p0 = [11.5, 60]  ;initial values
nparam = N_ELEMENTS(par_name)

 ;--- initial parameter choice ---
  myparam = DIALOG_INPUT(NFIELDS = nparam, $
           PROMPT = [par_name[0] + ' (' + units[0] + '): ', $
                     par_name[1] + ' (' + units[1] + '): '], $
           INITIAL = [STRTRIM(STRING(p0[0],format='(F6.2)'),2), $
	    	      STRTRIM(STRING(p0[1],format='(E8.1)'),2)], $
           TITLE = 'Input the initial parameter values: ')

  IF (myparam[0] EQ '') THEN GOTO, cancel ELSE p0 = DOUBLE(myparam)
    
;--- make the parinfo structure ---
parinfo = REPLICATE({fixed:0, limited:[1,0], limits:[0.0,0.0], parname:'', STEP:0D}, nparam);, relstep:0D}, nparam)
;parinfo = REPLICATE({fixed:0, limited:[1,0], limits:[0.0,0.0], parname:'', MPMAXSTEP:0D}, nparam);, STEP:0D}, nparam);, relstep:0D}, nparam)
parinfo.parname = par_name
parinfo[0].step = 0.5D; myparam[0]/5D
;parinfo[1].relstep = 1D ;doesn't affect the fit much
;parinfo[0].MPMAXSTEP = 1.5D

  ;free or fixed?
  parinfo.fixed = 0
  FOR fr = 0, nparam-1 DO BEGIN
    free = DIALOG_CHECKLIST(['free', 'fixed'], title = '     '+par_name[fr]+':     ', init = 1)
    fixed = WHERE(free EQ 1, nf)
    IF nf EQ 0 THEN GOTO, cancel
    parinfo[fr].fixed = fixed
  ENDFOR

parinfo[0].limits = [1d, 40d]
parinfo[1].limits = [0, 1]

parinfo[0].limited = 1
parinfo[1].limited = 0


t0=57625.  ;start of outburst
tl=57761.  ;end of dates that we use maxi data
chand=57789.4 ;date of chandra (or earlier observations)
src_time = Findgen(floor(tl-t0)+28)
tobs = chand-t0
src_time[N_ELEMENTS(src_time)-1] = tobs ;set the max time to chandra time


 ;--- fitting range choice ---
  select_range:
  myrange = DIALOG_INPUT(NFIELDS = 2, $
           PROMPT = ['Minimum (arcsec):', 'Maximum (arcsec):'], $
           INITIAL = [100D, 250D], $
           ;INITIAL = [MIN(data_x), MAX(data_x)], $
           TITLE = 'Enter the fitting range: ')

  IF (myrange[0] EQ '') THEN GOTO, cancel ELSE xran = FLOAT(myrange)
  IF (MIN(myrange) LT MIN(data_x) OR MAX(myrange) GT MAX(data_x)) THEN BEGIN
    msg = ['The range you provided out of the data rage.', 'Please select a valid range!']
    res = DIALOG_MESSAGE(msg, /cancel, /center)
    IF (res EQ 'OK') THEN GOTO, select_range ELSE GOTO, cancel
  ENDIF
    
  fit_range = WHERE(data_x GE xran[0] AND data_x LE xran[1])
  fit_x = data_x(fit_range)
  fit_y = data_y(fit_range)
  fit_e = data_e(fit_range)



hires = 5000.
time_hires = FINDGEN(hires)*(MAX(src_time)-MIN(src_time))/(hires-1) + MIN(src_time)
Foft_hires = INTERPOL(foft, src_time, time_hires)

xn_arr = parcloud_fit.x
weight = parcloud_fit.weight
thcl = parcloud_fit.thcl

cl_dist = DBLARR(51)
OPENR, lun, 'cloud_distance.dat', /GET_LUN
READF, lun, cl_dist
CLOSE, lun & FREE_LUN, lun

cl_dist = cl_dist[1:*]

;STOP
;result = SBP_MODEL(fit_x, p0)


model_fit = MPFITFUN('SBP_MODEL', fit_x, fit_y, fit_e, p0, parinfo=parinfo, $
    	    	    	perror=perror, dof=dof, bestnorm=chi)

OPLOT, fit_x, SBP_MODEL(fit_x, model_fit), color = 100, thick = 2    ; Overplot the model

PRINT, ''
PRINT, 'Best-Fit Parameters:'
FOR i = 0, nparam-1 DO $
    PRINT, par_name[i] + ' = ' + STRTRIM(model_fit[i], 2) + ' +/- ' + STRTRIM(perror[i], 2) $ 
    +	     ' ' + units[i]
PRINT, 'CHI-SQUARE/DOF = ' + STRTRIM(chi, 2) + ' / ' + STRTRIM(dof, 2)
PRINT, 'Fit Range = [' + STRTRIM(xran[0], 2) +', ' + STRTRIM(xran[1], 2) + '] arcsec'

;oplot, data_x, SBP_MODEL(data_x, [15,6e-9]), color = 250

cancel:

END
