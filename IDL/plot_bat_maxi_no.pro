pro plot_bat_maxi_no, drange, ps=ps, fname=namef, shobsdates=shobsdates, lowe=lowe

;This program plots MAXI and BAT light curves


IF NOT keyword_set(ps) THEN ps=0
IF NOT keyword_set(namef) THEN namef='maxi_bat_'+strtrim(string(drange[0]),1)+$
                                     '_'+strtrim(string(drange[1]),1)+'.eps'
IF NOT keyword_set(shobsdates) THEN shobsdates=0
IF NOT keyword_set(lowe) THEN lowe=0

;read data (need to change filenames at some point)

readcol,'J1634-473_g_lc_1day_all.dat', mjdm, m220,e220, m24,e24, m410, e410, m1020, e1020c


;swiftdays
swday=[57769.8,57771.8,57773.2,57777.9,57779.9,57781.9,57783.8,$
       57785.5,57787.0,57789.0]-57000.

;chandra day
chanday=57789.4-57000.

;sed_data=read_csv('maxi_hv0.csv')
;mjdm=sed_data.field1
;m220=sed_data.field2
;e220=sed_data.field3
;m24=sed_data.field4
;e24=sed_data.field5
;m410=sed_data.field6
;e410=sed_data.field7
;m1020=sed_data.field8
;e1020=sed_data.field9

  readcol,'4U1630-472.lc.txt', mjdb, batr, bate, y,d,sterr,syserr,df,tdel,tdelc,tdeldq
  
loadct,5
if ps then begin
   set_plot, 'ps' 
   device,/color
;   loadct,5
   device,/encapsulated
   device, filename = namef
   device, yoffset = 2
   device, ysize = 21.
   ;device, xsize = 12.0
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

multiplot, [1,3], mxtitle='Time (MJD-57000 days)', mxtitsize=1.2

;limit drange

ma=where((mjdm GE drange[0]) AND (mjdm LE drange[1]))
ba=where((mjdb GE drange[0]) AND (mjdb LE drange[1]))

;character size and plot symbol
cs=1.2
plotsym,0,/fill
nxr=drange-57000.
;MAXI 2-20 keV
;determine y range of plot

IF lowe THEN BEGIN
   mpl=m24
   epl=e24
   nytit='MAXI 2-4 keV rate'
ENDIF ELSE BEGIN
   mpl=m220
   epl=e220
   nytit='MAXI 2-20 keV rate'
ENDELSE


ymax=1.1*max(mpl[ma]+epl[ma])
nyr=[0.,ymax]
IF lowe THEN nyr=[0.,0.47]

ploterror, mjdm[ma]-57000, mpl[ma], epl[ma], yr=nyr, /xstyle, /ystyle, $
           /nohat, $
           chars=cs, psym=8, xr=nxr, background=backc, $
           ytitle=nytit,axiscolor=axc

IF shobsdates THEN BEGIN
   FOR i=0,n_elements(swday)-1 DO arrow, swday[i],!y.crange[1]*0.5,swday[i],!y.crange[1]*0.4,color=100,/data
   arrow, chanday,!y.crange[1]*0.5,chanday,!y.crange[1]*0.4,color=30,$
          /data,thick=2.5
ENDIF


multiplot

;BAT 15-50 keV
;determine y range of plot

ymaxb=1.1*max(batr[ba]+bate[ba])
;nyrb=[0.,ymaxb]
nyrb=[0.,0.022]

ploterror, mjdb[ba]-57000., batr[ba], bate[ba], yr=nyrb, /xstyle, /ystyle, $
           /nohat, $
           chars=cs, psym=8, xr=nxr, background=backc, $
           ytitle='BAT 15-50 keV rate',axiscolor=axc

multiplot

;MAXI RATIO

;eliminate negative values
pos=where((m410[ma] GE 0.) AND (m24[ma] GE 0.))
marat=m410[ma[pos]]/m24[ma[pos]]
emarat=sqrt((e24[ma[pos]]/m24[ma[pos]])^2.+(e410[ma[pos]]/m410[ma[pos]])^2.)/sqrt(2.) ;assuming independent
ymaxr=1.1*max(marat+emarat)
yminr=0.9*min(marat-emarat)

nyrrm=[yminr,ymaxr]
if ymaxr GT 5. THEN nyrrm=[0.2,2.9]

;ploterror, mjdm[ma[pos]], marat, marat*emarat, yr=nyrrm, /xstyle, /ystyle, /nohat, chars=cs, xr=nxr, psym=8, ytitle='MAXI 4-10/2-4 keV'

;multiplot

;BAT / MAXI

;bring dates to the same order
mjdb2=mjdb[ba]-0.5
mjdm2=mjdm[ma]
;get indices, start with first index

cont=0
mi=0 ;maxi index
WHILE cont EQ 0 DO BEGIN
   bi=where(floor(mjdb2) eq floor(mjdm2[mi]))
   IF bi eq -1 THEN BEGIN
      mi=mi+1
      IF mi eq n_elements(mjdm2) THEN BEGIN
         cont=1
         mi=mi-1
         ENDIF
      ENDIF ELSE cont=1
ENDWHILE

FOR i=mi+1, n_elements(ma)-1 DO BEGIN
   bix=where(floor(mjdb2) eq floor(mjdm2[i]))
   IF bix NE -1 THEN BEGIN
      bi=[bi,bix]
      mi=[mi,i]
   ENDIF
ENDFOR


bmrat=batr[ba[bi]]/m24[ma[mi]]
ebmrat=sqrt((bate[ba[bi]]/batr[ba[bi]])^2.+(e24[ma[mi]]/m24[ma[mi]])^2.)/sqrt(2.) ;assuming independent
ymaxbr=1.1*max(bmrat+ebmrat)
yminbr=0.9*min(bmrat-ebmrat)

nyrbr=[yminbr,ymaxbr]
IF ymaxbr GT 0.35 THEN nyrbr=[0.,.29]

ploterror, mjdm2[mi]-57000., bmrat, bmrat*ebmrat, yr=nyrbr, /xstyle, /ystyle, /nohat, chars=cs, xr=nxr, psym=8, ytitle='BAT(15-50) / MAXI(2-4)',axiscolor=axc

multiplot,/default

IF ps THEN BEGIN
   device,/close
   set_plot,'x'
ENDIF

END


