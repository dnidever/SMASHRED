;+
;
; SMASHRED_MORPHCUTS
;
; This program performs star galaxy separation using
; the sharpness morphology parameter.
;
; INPUTS:
;  obj      Structure with SMASH photometry.
;  /pl      Make some plots
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  gdsharp  The index array of sources passing the sharpness cuts.
;  fitstr   The structure with the fitting parameters.
;
; USAGE:
;  IDL>smashred_morphcuts,obj,gdsharp
;
; By D.Nidever  Jan 2017
;-

pro smashred_morphcuts,obj,gdsharp,fitstr,pl=pl,stp=stp

; Apply morphological star/galaxy cuts

undefine,gdsharp
undefine,fitstr
  
; Not enough inputs
if n_elements(obj) eq 0 then begin
  print,'Syntax - smashred_morphcuts,obj,gdsharp,fitstr,pl=pl,stp=stp'
  return
endif
  
;gd = where(obj.ndet gt 5 and obj.chi lt 3 and obj.sharp lt 50 and obj.g lt 50,ngd)
;obj2 = obj[gd]
;
;bindata,obj2.g,obj2.sharp,xbin,numbin,bin=0.2,/hist
;gdind = where(numbin gt 5,ngdind)
;bindata,obj2.g,obj2.sharp,xbin,medbin,bin=0.2,/med
;bindata,obj2.g,obj2.sharp,xbin,sigbin,bin=0.2,/mad
;xbin2 = xbin[gdind]
;
;sigcoef = robust_poly_fitq(xbin[gdind],alog10(sigbin[gdind]),2)
;interp,xbin2,medbin[gdind],obj2.g,medsharp
;sigsharp = 10^poly(obj2.g,sigcoef)
;gd = where(abs(obj2.sharp-medsharp) lt 3*sigsharp,ngd)
;
;hess,obj2.g,obj2.sharp,dx=0.1,dy=0.05,xr=[13,25.5],yr=[-2,2],/log
;oplot,xbin2,medbin[gdind]
;oplot,xbin2,10^poly(xbin2,sigcoef)+medbin[gdind]
;oplot,xbin2,-(10^poly(xbin2,sigcoef))+medbin[gdind]
;
;; Redo
;bindata,obj2[gd].g,obj2[gd].sharp,xbin,numbin2,bin=0.2,/hist
;gdind2 = where(numbin gt 5,ngdind2)
;bindata,obj2[gd].g,obj2[gd].sharp,xbin,medbin2,bin=0.2,/med
;bindata,obj2[gd].g,obj2[gd].sharp,xbin,sigbin2,bin=0.2,/mad
;xbin2 = xbin[gdind2]
;
;; the median still deviates upward at the faintest magnitudes
;; how to fix that?
;
;; should I try doing this fitting for each exposure separately??
;
;; MAYBE TRY B-SPLINES INSTEAD!!!!
;
;dx = 0.1
;bin = 0.20
;g0 = 14.0
;g1 = 25.5
;nbin = (g1-g0)/dx + 1
;tstr = replicate({gbin:0.0,sharp:0.0,scatter:0.0},nbin)
;gtemp = lindgen(nbin)*dx+g0+0.5*dx
;tstr.gbin = gtemp
;tags = tag_names(tstr)
;
;nord = 2 ;3 ;4
;bkspace = 2.0 ;0.3 ; 50
;gd = where(obj.g lt 50 and abs(obj.sharp) lt 1,ngd)
;xx = obj[gd].g
;yy = obj[gd].sharp
;invvar = fltarr(ngd)+1.0
;dum = bspline_iterfit(xx,yy,invvar=invvar,nord=nord,bkspace=bkspace,yfit=model)
;diff = yy-model
;sig = mad(diff,/zero)
;;if keyword_set(pl) then begin
;  plot,xx,diff,ps=3,xr=[-1,4],yr=[-4,4]*sig
;  oplot,[-10,10],[0,0],linestyle=2
;;endif

; It also deviates at the faintest mags because of the giant "blob"
; at higher sharp values.
; Maybe I should try to follow the "ridgeline".

; make the density diagram and find the maxima at any color
; step from left to right and only take the maximum that is closest
; to the one on the left (the last one).

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
  match,xind,i,ind,ind0,/sort,count=nmatch
  sharp = obj[ind].sharp
  sig0 = mad(sharp-sharpmed,/zero)
  ind2 = where(abs(sharp-sharpmed) lt 2.5*sig0,nind2)
  sig = mad(sharp[ind2]-sharpmed,/zero)
  sigarr[i] = sig
endfor
if keyword_set(pl) then begin
  oplot,xarr,poly(xarr,sharpcoef)+sigarr,co=150,linestyle=2
  oplot,xarr,poly(xarr,sharpcoef)-sigarr,co=150,linestyle=2
endif
  
; now measure scatter, the gaussian-fitting sigma
; grows too big because of the "blob", need something that grows
; more slowly at the end

; Apply the cut
;--------------
ind = where(obj.g lt 50,nind)
sharpdiff = obj.sharp*0+99.99
sharpdiff[ind] = obj[ind].sharp - poly(obj[ind].g,sharpcoef)
; Use 0.2 threshold at bright end
; at faint end use 2*sigma
thresharr = sigarr > 0.2
sharpthresh = 2*sharpdiff*0+1.0
; Interpolate
interp,xarr,thresharr,obj[ind].g,sharpthresh1
sharpthresh[ind] = sharpthresh1
; Apply the sharpness cut
gdsharp = where(abs(sharpdiff) lt sharpthresh,ngdsharp,comp=bdsharp,ncomp=nbdsharp)

; Make some plots
if keyword_set(pl) then begin
  hess,obj.g,sharpdiff,dx=0.2,dy=0.05,xr=[14,26],yr=[-1,1]
  oplot,[0,30],[0,0],linestyle=2
  oplot,xarr,sigarr,linestyle=1,co=150
  oplot,xarr,-sigarr,linestyle=1,co=150
  oplot,xarr,thresharr,linestyle=2,co=250
  oplot,xarr,-thresharr,linestyle=2,co=250
endif

; The fitting parameters
fitstr = {gmag:xarr,medsharp:reform(pararr[*,1]),gpar:pararr,$
          sharpcoef:sharpcoef,sigma:sigarr,thresh:thresharr}

;hess,obj[gdsharp].g,sharpdiff[gdsharp],dx=0.2,dy=0.05,xr=[14,26],yr=[-1,1]
;hess,obj[gdsharp].g-obj[gdsharp].i,obj[gdsharp].g,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[26,15],/log


if keyword_set(stp) then stop

end
