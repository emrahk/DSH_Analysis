pro plotmaxi_triggerx, drange, ps=ps, fname=namef, lowe=lowe, efold=folde, today=dtoday

IF NOT keyword_set(ps) THEN ps=0
IF NOT keyword_set(namef) THEN namef='plotmaxi_triggerx.eps'
IF NOT keyword_set(lowe) THEN lowe=0
IF NOT keyword_set(dtoday) Then dtoday=59819.5
device,decomposed=0

readcol,'J1634-473_g_lc_1day_all.dat', mjdm, m220,e220, m24,e24, m410, e410, m1020, e1020

;sed_data=read_csv('glcbin24.0h_regbg_hv0.csv')
;mjdm=sed_data.field1
;m220=sed_data.field2
;e220=sed_data.field3;
;m24=sed_data.field4
;e24=sed_data.field5
;m410=sed_data.field6
;e410=sed_data.field7
;m1020=sed_data.field8
;e1020=sed_data.field9

xx=where((mjdm gt drange[0]) AND (mjdm lt drange[1]))

IF lowe THEN BEGIN
   rate=m24[xx]
   e_rate=e24[xx]
ENDIF ELSE BEGIN
   rate=m410[xx]
   e_rate=e410[xx]
ENDELSE

dnew=mjdm[xx]-mjdm[xx[0]]

hr=m410/m24
hrx=hr[xx]

;stop
IF lowe THEN yrx=[0.01,0.5] ELSE yrx=[0.01,0.5]
derr=fltarr(n_elements(dnew))
xtit='Date (MJD-'+strtrim(string(mjdm[xx[0]]),1)+')'
IF lowe THEN ytit='MAXI 2-4 keV rate' ELSE ytit='MAXI 4-10 keV rate'
ploterror, dnew, rate, derr, e_rate, xrange=[0.,max(dnew)+40.], $
           yrange=yrx,/nohat,$
           /xstyle,/ystyle,xtitle=xtit, ytitle=ytit,$
           psym=5, chars=1,/ylog,err_col=255


swdays=[55392.5, 57769.8, 57771.8, 57773.2,57777.9]

FOR i=0,4 DO BEGIN
                                ;find count rates
   yy=where((mjdm GE swdays[i]-1.) AND (mjdm LE swdays[i]+1.),nxx)
   IF nxx GT 1 THEN BEGIN
      m24l=avg(m24[yy])
      m410l=avg(m410[yy])
   ENDIF ELSE BEGIN
      m24l=m24[yy]
      m410l=m410[yy]
   ENDELSE
 print,mjdm[yy],m24l,m410l

   
   lis=0
   IF i eq 1 THEN lis=1
   IF i eq 0 THEN lis=2
   
IF lowe THEN oplot,[!x.crange[0],!x.crange[1]],[m24l,m24l],col=0,line=lis ELSE $
   oplot,[!x.crange[0],!x.crange[1]],[m410l,m410l],line=lis,col=0
;55392.5
ENDFOR


cond1 = 'No' ;Set initial prompt response.

;IF NOT soft THEN PRINT,'Finding the hardening transition' ELSE $
;  PRINT,'Finding the softening transition'

;IF NOT SOFT THEN BEGIN
WHILE cond1 EQ 'No' do begin

  print,'Please set the range of line!'
  print,'Move the mouse to the startpoint of region and click'
  cursor, x1, y1, /down
  oplot,[x1,x1],10.^(!y.crange),line=0,col=0
  print,'Move the mouse to the endpoint of region and click'
  if !mouse.button ne 4 then cursor, x2, y2, /down
  oplot,[x2,x2],10.^(!y.crange),line=0,col=0
  wait,0.5
  rt=[x1,x2]
  cond1 = DIALOG_MESSAGE('Are you happy with the range you chose?', /Ques)

ENDWHILE

ti=where((dnew gt rt[0]) AND (dnew lt rt[1]))
it=rate[ti]
ite=e_rate[ti]
tt=dnew[ti]

hhh=findgen(max(dnew)*5)

IF NOT keyword_set(folde) THEN BEGIN
   
   aexp=[10.,-0.125]
   l1fit=exponential_fit(tt,it,guess=aexp,errors=ite,sigma=sigma_l1)
   print,l1fit
ENDIF ELSE BEGIN
                                ;some magic numbers required
                                ;ds=141.5, there at 0.3
   ;norm*e(141.*fit)=0.3/

   l1fit=[0.3/exp(64.5D*folde),folde]
   print,l1fit
ENDELSE


oplot,hhh,l1fit[0]*exp(hhh*l1fit[1]),col=0

;find where it hits target
IF lowe THEN tarr=0.001 ELSE tarr=0.0077

;oplot,[!x.crange[0],!x.crange[1]],tarr*[1.4,1.4],col=0
;oplot,[!x.crange[0],!x.crange[1]],tarr*[0.2,0.2]/0.7,col=0


;dtarget1=alog(tarr*1.4/l1fit[0])/(l1fit[1])
;dtarget2=alog(tarr*0.3/l1fit[0])/(l1fit[1])
;trigval1=l1fit[0]*exp(l1fit[1]*(dtarget1-15.))
;trigval2=l1fit[0]*exp(l1fit[1]*(dtarget2-15.))
;print,dtarget1,dtarget2
;print,trigval1,trigval2
;oplot,[!x.crange[0],!x.crange[1]],trigval1*[1.,1.],col=0,line=2
;oplot,[!x.crange[0],!x.crange[1]],trigval2*[1.,1.],col=0,line=2
;print,-1./l1fit[1]


;oplot,58416-[mjdm[xx[0]],mjdm[xx[0]]],10^(!y.crange),col=0
;print,58416-[mjdm[xx[0]],mjdm[xx[0]]],10^(!y.crange)

arrow, dtoday-mjdm[xx[0]], 0.2, dtoday-mjdm[xx[0]], 0.1, col=0,/data

;dd=findgen(20)
;oplot,54.+dd,0.1*exp(-dd/4.),col=0
;print,0.2*exp(-dd/4.)

;boxc,58416-mjdm[xx[0]], 10^(!y.crange[0]),!x.crange[1],10^(!y.crange[1]),100
;xyouts, 0.9, 0.2,'No Chandra Visibility',color=0,orientation=90.,$
;         size=2.5,/normal
;xyouts, 0.3, 0.3,'Chandra Target',color=0,$
;         size=2.5,/normal

;swift

plotsym, 0, /fill

;oplot,[58397.6,58397.6]-mjdm[xx[0]],[0.14,0.14],psym=8, color=0


if ps then begin
   set_plot, 'ps' 
   device,/color
;   loadct,5
   device,/encapsulated
   device, filename = namef
   device, yoffset = 2
   device, ysize = 12.
   ;device, xsize = 12.0
   !p.font=0
   device,/times
;   backc=0
endif

;postscript later

;1.47 e-9 = 0.14 cts/s absorbed. 2.2e-9 2-10 keV absorbed, 6.6e-9
;unabsorbed.

;multiply, 6.6e-9/0.14

IF ps THEN BEGIN
   device,/close
   set_plot,'ps'
ENDIF

; NOTE TO SELF for next proposal use Tomsick, the exponential time
; scale decay is around 3.5 - 4 days. Need a larger Maxi/swift trigger
;                             time. swift flux should have been higher
;                             either for the worst case scenario

END
