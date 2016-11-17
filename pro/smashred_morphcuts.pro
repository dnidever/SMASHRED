pro smashred_morphcuts,obj

; Apply morphological star/galaxy cuts

gd = where(obj.ndet gt 5 and obj.chi lt 3 and obj.sharp lt 50 and obj.g lt 50,ngd)
obj2 = obj[gd]

bindata,obj2.g,obj2.sharp,xbin,numbin,bin=0.2,/hist
gdind = where(numbin gt 5,ngdind)
bindata,obj2.g,obj2.sharp,xbin,medbin,bin=0.2,/med
bindata,obj2.g,obj2.sharp,xbin,sigbin,bin=0.2,/mad
xbin2 = xbin[gdind]

sigcoef = robust_poly_fitq(xbin[gdind],alog10(sigbin[gdind]),2)
interp,xbin2,medbin[gdind],obj2.g,medsharp
sigsharp = 10^poly(obj2.g,sigcoef)
gd = where(abs(obj2.sharp-medsharp) lt 3*sigsharp,ngd)

hess,obj2.g,obj2.sharp,dx=0.1,dy=0.05,xr=[13,25.5],yr=[-2,2],/log
oplot,xbin2,medbin[gdind]
oplot,xbin2,10^poly(xbin2,sigcoef)+medbin[gdind]
oplot,xbin2,-(10^poly(xbin2,sigcoef))+medbin[gdind]

; Redo
bindata,obj2[gd].g,obj2[gd].sharp,xbin,numbin2,bin=0.2,/hist
gdind2 = where(numbin gt 5,ngdind2)
bindata,obj2[gd].g,obj2[gd].sharp,xbin,medbin2,bin=0.2,/med
bindata,obj2[gd].g,obj2[gd].sharp,xbin,sigbin2,bin=0.2,/mad
xbin2 = xbin[gdind2]

; the median still deviates upward at the faintest magnitudes
; how to fix that?

; should I try doing this fitting for each exposure separately??

; MAYBE TRY B-SPLINES INSTEAD!!!!


dx = 0.1
bin = 0.20
g0 = 14.0
g1 = 25.5
nbin = (g1-g0)/dx + 1
tstr = replicate({gbin:0.0,sharp:0.0,scatter:0.0},nbin)
gtemp = lindgen(nbin)*dx+g0+0.5*dx
tstr.gbin = gtemp
tags = tag_names(tstr)

nord = 2 ;3 ;4
bkspace = 2.0 ;0.3 ; 50
gd = where(obj.g lt 50 and abs(obj.sharp) lt 1,ngd)
xx = obj[gd].g
yy = obj[gd].sharp
invvar = fltarr(ngd)+1.0
dum = bspline_iterfit(xx,yy,invvar=invvar,nord=nord,bkspace=bkspace,yfit=model)
diff = yy-model
sig = mad(diff,/zero)
;if keyword_set(pl) then begin
  plot,xx,diff,ps=3,xr=[-1,4],yr=[-4,4]*sig
  oplot,[-10,10],[0,0],linestyle=2
;endif

; It also deviates at the faintest mags because of the giant "blob"
; at higher sharp values.
; Maybe I should try to follow the "ridgeline".

; make the density diagram and find the maxima at any color
; step from left to right and only take the maximum that is closest
; to the one on the left (the last one).

undefine,dum,im
hess,obj[gd].g,obj[gd].sharp,dum,im,dx=0.2,dy=0.05,xr=[14,25.6],yr=[-1,1],/noplot,xarr=xarr,yarr=yarr
nx = n_elements(xarr)
ybin = fltarr(nx)
pararr = fltarr(nx,3)
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

; now measure scatter, the gaussian-fitting sigma
; grows too big because of the "blob", need something that grows
; more slowly at the end

stop

; Measure sigs and scatters
si = sort(xx)
nbinsigs = 100 ;200
nsigs = floor(ngd/nbinsigs)
sigs = fltarr(nsigs)
scatters = sigs
mngi = sigs
mederr = sigs
for j=0L,nsigs-1 do begin
  lo = nbinsigs*j  &  hi=(lo+nbinsigs-1) < (ngd-1)
  sigs[j] = mad(diff[si[lo:hi]],/zero)
  scatters[j] = sqrt( sqrt(median(diff[si[lo:hi]]^2))^2 - sqrt(median(colerr[si[lo:hi]]^2))^2 )
  mngi[j] = mean(xx[si[lo:hi]])
  mederr[j] = median(colerr[si[lo:hi]])
endfor

; Smooth
smsigs = gsmooth(sigs,3)
smscatters = gsmooth(scatters,3)
smmederr = gsmooth(mederr,3)
  
; Make the templates
interp,xx[si],model[si],gitemp,coltemp
interp,mngi,smsigs,gitemp,sigstemp
interp,mngi,smscatters,gitemp,scatterstemp
interp,mngi,smmederr,gitemp,mederrtemp

tind1 = where(tags eq strupcase(bands[mind[i]])+'I',ntind1)
tstr.(tind1[0]) = coltemp
tind2 = where(tags eq strupcase(bands[mind[i]])+'ISIG',ntind2)
tstr.(tind2[0]) = sigstemp
tind3 = where(tags eq strupcase(bands[mind[i]])+'ISCATTER',ntind3)
tstr.(tind3[0]) = scatterstemp

tstr.cols[mind[i]] = coltemp
tstr.sigs[mind[i]] = sigstemp
tstr.scatters[mind[i]] = scatterstemp
tstr.mederr[mind[i]] = mederrtemp
tstr.n = ngd

if keyword_set(stp) then stop

end
