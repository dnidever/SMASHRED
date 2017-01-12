;+
;
; SMASHRED_MORPHCUTS
;
; This program performs star galaxy separation using
; the sharpness morphology parameter.
;
; INPUTS:
;  obj            Structure with SMASH photometry.
;  =brightthresh  The sharpness threshold to use for the bright
;                   sources.  By default this is 0.2.
;  =nsigthresh    The number of sigmas to use for the sharpness cut
;                   at the faint end.  The default is 2.0.
;  /pl            Make some plots
;  /stp           Stop at the end of the program.
;
; OUTPUTS:
;  ind            The index array of sources passing the sharpness cuts.
;  fitstr         The structure with the fitting parameters.
;
; USAGE:
;  IDL>smashred_morphcuts,obj,gdsharp
;
; By D.Nidever  Jan 2017
;-

pro smashred_morphcuts,obj,ind,fitstr,brightthresh=brightthresh0,nsigthresh=nsigthresh0,pl=pl,stp=stp

; Apply morphological star/galaxy cuts

undefine,gdsharp
undefine,fitstr
  
; Not enough inputs
if n_elements(obj) eq 0 then begin
  print,'Syntax - smashred_morphcuts,obj,ind,fitstr,brightthresh=brightthresh,nsigthresh=nsigthresh,pl=pl,stp=stp'
  return
endif

; Defaults
if n_elements(brightthresh0) gt 0 then brightthresh=brightthresh0[0]>0.01 else brightthresh=0.2
if n_elements(nsigthresh0) gt 0 then nsigthresh=nsigthresh0[0]>0.5 else nsigthresh=2.0


; make the density diagram and find the maxima at any color
; step from left to right and only take the maximum that is closest
; to the one on the left (the last one).


; Loop through magnitude bins and fit Gaussian to the sharpness distribution
undefine,dum,im
gdobj = where(obj.g lt 50 and abs(obj.sharp) lt 1 and obj.ndet gt 5,ngd)
dx = 0.2
dy = 0.05
hess,obj[gdobj].g,obj[gdobj].sharp,dum,im,dx=dx,dy=dy,xr=[14,25.6],yr=[-1,1],/noplot,xarr=xarr,yarr=yarr
nx = n_elements(xarr)
ybin = fltarr(nx)
pararr = fltarr(nx,3)+999999.
for i=0,nx-1 do begin
  line = reform(im[i,*])
  dln_maxmin,line,minind,maxind
  if i eq 0 then begin
    hiind = first_el(maxloc(line[maxind]))
    bestind = maxind[hiind]
    ybin[i] = yarr[bestind]
  endif else begin
    ymax = yarr[maxind]
    diff = abs(ymax-ybin[i-1])
    clsind = first_el(minloc(diff))
    bestind = maxind[clsind]
    ybin[i] = yarr[bestind]
  endelse
  ; fit Gaussian to get scatter
  parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},3)
  parinfo[0].limited=1 & parinfo[0].limits=[0.9,1.1]*line[bestind]
  parinfo[1].limited=1 & parinfo[1].limits=[-0.05,0.05]+ybin[i]
  parinfo[2].limited=1 & parinfo[2].limits=[0.05,max(yarr)]
  estimates = [line[bestind],ybin[i],0.1]
  yfit = mpfitpeak(yarr,line,par,nterms=3,estimates=estimates,parinfo=parinfo,status=status)
  pararr[i,*] = par
  ;stop
endfor

; this works well!!!!

if keyword_set(pl) then begin
  hess,obj[gdobj].g,obj[gdobj].sharp,dum,im,dx=0.2,dy=0.05,xr=[14,25.6],yr=[-1,1],xarr=xarr,yarr=yarr
  oplot,xarr,pararr[*,1],ps=-1
endif
  
; Polynomial fit
gdpar = where(xarr le 24 and pararr[*,1] lt 1e4,ngd)
sharpcoef = robust_poly_fitq(xarr[gdpar],reform(pararr[gdpar,1]),2)
x = scale_vector(findgen(100),0,30)
if keyword_set(pl) then oplot,x,poly(x,sharpcoef),co=250

; Loop over mags and measure sigma with MAD
sigarr = fltarr(nx)
xind = floor((obj.g-xarr[0])/dx)
for i=0,nx-1 do begin
  gmed = xarr[i]
  sharpmed = poly(gmed,sharpcoef)
  ;ind = where(obj.g ge xarr[i]-0.5*dx and obj.g lt xarr[i]+0.5*dx and abs(obj.sharp) lt 1,nind)
  match,xind,i,ind1,ind0,/sort,count=nmatch
  sharp = obj[ind1].sharp
  sharpdiff = sharp-sharpmed
  
  ; Perform sigma on negative and positive values separately
  ; and use the lowest one
  posind = where(sharpdiff ge 0,npos)
  negind = where(sharpdiff le 0,nneg)
  sigpos = mad(sharpdiff[posind],/zero)
  signeg = mad(sharpdiff[negind],/zero)
  sig0 = mad(sharpdiff,/zero)
  sig = (sigpos < signeg)

  ; Iterative sigma on all points
  ;sig0 = mad(sharpdiff,/zero)
  ;ind2 = where(abs(sharpdiff) lt 2.5*sig0,nind2)
  ;sig = mad(sharpdiff[ind2],/zero)

  sigarr[i] = sig
endfor
if keyword_set(pl) then begin
  oplot,xarr,poly(xarr,sharpcoef)+sigarr,co=150,linestyle=2
  oplot,xarr,poly(xarr,sharpcoef)-sigarr,co=150,linestyle=2
endif
  

; Apply the cut
;--------------
gd = where(obj.g lt 50,ngd)
sharpdiff = obj.sharp*0+99.99
sharpdiff[gd] = obj[gd].sharp - poly(obj[gd].g,sharpcoef)
; Use 0.2 threshold at bright end
; at faint end use 2*sigma
thresharr = nsigthresh*sigarr > brightthresh
sharpthresh = sharpdiff*0+1.0
; Interpolate
interp,xarr,thresharr,obj[gd].g,sharpthresh1
sharpthresh[gd] = sharpthresh1
; Apply the sharpness cut
ind = where(abs(sharpdiff) lt sharpthresh,nind,comp=cind,ncomp=ncind)

; Make some plots
if keyword_set(pl) then begin
  hess,obj.g,sharpdiff,dx=0.2,dy=0.05,xr=[14,26],yr=[-1,1]
  oplot,[0,30],[0,0],linestyle=2
  oplot,xarr,sigarr,linestyle=0,co=150
  oplot,xarr,-sigarr,linestyle=0,co=150
  oplot,xarr,thresharr,linestyle=2,co=250
  oplot,xarr,-thresharr,linestyle=2,co=250
endif

; The fitting parameters
fitstr = {gmag:xarr,medsharp:reform(pararr[*,1]),gpar:pararr,$
          sharpcoef:sharpcoef,sigma:sigarr,thresh:thresharr}

;hess,obj[ind].g,sharpdiff[ind],dx=0.2,dy=0.05,xr=[14,26],yr=[-1,1]
;hess,obj[ind].g-obj[ind].i,obj[ind].g,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[26,15],/log


if keyword_set(stp) then stop

end
