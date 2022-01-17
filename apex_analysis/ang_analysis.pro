function ang_analysis, inpimg, cdelt, rin, rout, theta, del_theta, $
                       offset=setoff, plot=plot, ps=ps, fname=namef
  
; This program calculates the brightness of a wedge with inner radius
; rin, outer radius rout from the given center, theta-del_theta/2 to
;  theta+del_theta/2 where theta is the angle from the x axis.
;
; INPUTS
;
; inpimg: input image array
; cdelt: converts pixels to arcseconds
; rin: inner radius of the wedge in arcsecond
; rout: outer radius of the wedge in arcsecond
; theta: mid angle from X axis that passes through the center
; del_theta: delta theta to construct the wedge
;
; OUTPUTS
;
; returns sum of pixels in the wedge
;
; OPTIONAL INPUTS
;
; offset: (X, Y offset from the center of image in arcsecond
; plot: if set, plot the image with the wedge
; ps: if set output a postscript file
; fname: if set, change the default output postscripto to namef
;
; USES
;
; IDL image routines
;
; USED BY
;
; NONE
;
; LOGS
;
; Created by E. Kalemci, Jan 2022
;
; Magic numbers exist. Generalize
; Negative angles problem

IF NOT keyword_set(ps) THEN ps=0
IF NOT keyword_set(namef) THEN namef='plotwedge.eps'
IF NOT keyword_set(setoff) THEN setoff=[0.,0.]
IF NOT keyword_set(plot) THEN plot=0

rin=rin/cdelt
rout=rout/cdelt


IF plot THEN BEGIN
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

   ;first plot window to overplot annulus

   !x.margin=[0.,0.]
   !y.margin=[0.,0.]
   plot,[0,0],[0,0],xr=[0,280],yr=[0,280],/xstyle,/ystyle ;magic number!!!! fix later
   tv, inpimg

   ;oplot wedge

   angs=(theta-(del_theta/2.)+findgen(del_theta))*!PI/180.

   oplot, 140.+rin*cos(angs),140.+rin*sin(angs),$
                                  color=255, thick=1.5
   oplot, [140.+rin*cos(angs[del_theta-1]),140.+rout*cos(angs[del_theta-1])],$
          [140.+rin*sin(angs[del_theta-1]),140.+rout*sin(angs[del_theta-1])],$
          color=255, thick=1.5

   oplot, 140.+rout*cos(angs),140.+rout*sin(angs),$
                                  color=255, thick=1.5

   oplot, [140.+rin*cos(angs[0]),140.+rout*cos(angs[0])],$
          [140.+rin*sin(angs[0]),140.+rout*sin(angs[0])],$
          color=255, thick=1.5
ENDIF

;calculate sum

res=0.
sz=size(inpimg)
FOR i=0, sz[1]-1 DO BEGIN
   FOR j=0, sz[2]-1 DO BEGIN
      r=sqrt(((i-140.))^2.+((j-140.))^2.)
      tht=atan((j-140.)/(i-140.))*180./!PI
      IF (((j-140.) LT 0.) AND ((i-140.) LT 0.)) THEN tht=180.+tht
      IF (((j-140.) LT 0.) AND ((i-140.) GT 0.)) THEN tht=360.+tht
      IF (((j-140.) GT 0.) AND ((i-140.) LT 0.)) THEN tht=180.+tht
      IF ((r GE rin) AND (r LT rout) AND $
          (tht GE theta-del_theta/2.) AND (tht LT theta+del_theta/2.)) THEN $
             BEGIN
         res=res+inpimg[i,j]
         IF plot THEN oplot, [i,i]+0.5,[j,j]+0.5,psym=3, color=255
      ENDIF
   ENDFOR
ENDFOR


return,res

END

   
