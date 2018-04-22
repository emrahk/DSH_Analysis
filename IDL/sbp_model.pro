FUNCTION SBP_MODEL, theta_arr, p    ;p = [D, norm]

COMMON myVars, cl_dist, weight, thcl, time_hr, Foft_hr, t_obs

c = 2.998D10 ;cm/s
kpc = 3.086D21 ;cm
day = 86400D ;s
alpha=3.    ;parcloud.alpha[0]
beta=3.     ;parcloud.beta[0]
meanE=3.    ;parcloud.meanE[0]

p = DOUBLE(p)
weight = DOUBLE(weight)
thcl = DOUBLE(thcl)
time_hr = DOUBLE(time_hr)
Foft_hr = DOUBLE(Foft_hr)
t_obs = DOUBLE(t_obs)

n_data = N_ELEMENTS(theta_arr)
intensity = DBLARR(n_data)
int_n = DBLARR(n_data)
intensity_abs = DBLARR(n_data)
x_arr = cl_dist/P[0]
xn_arr = 0D
wt_arr = 0D

n_cl = N_ELEMENTS(cl_dist)

;----- divide all the clouds into 10pc thickness -----
FOR cl = 0, n_cl-1 DO BEGIN

  thcl_ind = FLOOR(thcl[cl]/0.01)
  xn_arr = [xn_arr, x_arr[cl]-FINDGEN(thcl_ind)*0.01/P[0]]
  wt_arr = [wt_arr, REPLICATE(weight[cl], thcl_ind)]
  
ENDFOR
  
xn_arr = xn_arr[1:*]
wt_arr = wt_arr[1:*]

FOR n = 0, N_ELEMENTS(xn_arr)-1 DO BEGIN

  xn = DOUBLE(xn_arr[n])    ;of thickness 10pc
  
  ;for 1 cloud at xn
  FOR th = 0, N_ELEMENTS(theta_arr)-1 DO BEGIN
    theta = DOUBLE(theta_arr[th])
    dt = (theta*!PI/(180D*3600D))^2 * xn*P[0]*kpc/(2*c*day*(1-xn))
    flux_t = t_obs - dt

    t_ind = MIN(WHERE(time_hr GE flux_t))

    int_n = Foft_hr[t_ind]*$
       	    	    ((theta/(1-xn)/1000.)^(-alpha))*meanE^(-beta)/(1-xn)^2.
		    
    intensity[th] += int_n
    
    intensity_abs[th] += P[1]*wt_arr[n]*int_n * 1d-10

  ENDFOR

ENDFOR	;for all clouds


RETURN, intensity_abs
END
