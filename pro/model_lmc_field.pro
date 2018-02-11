function lmc_hess_model,x,par,mwmodel=mwmodel,lmcsynth=lmcsynth,compim=compim,dx=dx,dy=dy,xr=xr,yr=yr,$
                        gcoef=gcoef,icoef=icoef,smlen=smlen,mask=mask

; Construct the model
; thindisk_scale, thickdisk_scale, halo_scale lmcdist, lmcdistspread, lmcscale
mwim = mwmodel[*,*,0] * par[0]  ; thin disk
mwim += mwmodel[*,*,1] * par[1] ; thick disk
mwim += mwmodel[*,*,2] * par[2] ; halo

; Add mean distance and distance spread
distmod = par[3]+randomn(seed,n_elements(lmcsynth))*par[4]
gobs = lmcsynth.gobs+distmod
iobs = lmcsynth.iobs+distmod
; Add uncertainties
gobs += randomn(seed,n_elements(lmcsynth))*10.0^poly(gobs,gcoef)
iobs += randomn(seed,n_elements(lmcsynth))*10.0^poly(iobs,icoef)

undefine,dum,lmcim
hess,gobs-iobs,gobs,dum,lmcim,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
lmcim *= par[5]  ; LMC scaling

; Apply completeness
im0 = (mwim + lmcim) * compim
im = gsmooth(im0,smlen)   ; smooth

; Mask the image
if n_elements(mask) gt 0 then im*=mask

return,im

end

;-----------

pro model_lmc_field,field

; Model the hess diagram of an outer LMC field

rootdir = smashred_rootdir()+'cp/red/photred/'
if n_elements(version) eq 0 then version='v6'
catdir = rootdir+'catalogs/final/'+version+'/'

; Not enough inputs
if n_elements(field) eq 0 then begin
  print,'Syntax - model_lmc_field,field'
  return
endif

; Restore the stars file
stars = mrdfits(catdir+'stars1/'+field+'_allobj_stars.fits.gz',1)
;stars = mrdfits(catdir+field+'_combined_allobj_deep.fits.gz',1)
; Make star/galaxy cuts
gd = where(stars.g lt 50 and stars.i lt 50,ngd)
;gd = where(stars.g lt 50 and stars.i lt 50 and stars.r lt 50 and abs(stars.sharp) lt 1 and stars.chi lt 3 and stars.prob gt 0.7,ngd)
stars = stars[gd]
; completeness
;ast = mrdfits(catdir+'stars1/'+field+'_complete_stars.fits.gz',1)
ast = mrdfits(catdir+'ast/'+field+'_complete.fits.gz',1)
; the "stars" ASTs look horrible



; MAYBE USE "DEEP" STARS CATALOG??

; Restore Thomas' model
model = importascii(rootdir+'mwmodels/dartmouth/CMD_'+strlowcase(field),/header)
nmodel = n_elements(model)

; Restore the LMC synthetic isochrone
;synth = mrdfits(rootdir+'lmchalo/parsec_10gyr_-1.00_synth.fits.gz',1)
;synth = mrdfits(rootdir+'lmchalo/parsec_7gyr_-1.00_synth.fits.gz',1)
synth = mrdfits(rootdir+'lmchalo/parsec_7gyr_-0.50_synth.fits.gz',1)
; no distance modulus
synth.gobs -= synth.dmod
synth.iobs -= synth.dmod
synth.dmod = 0

;DECam_u, R_u = 3.9631
;DECam_g, R_g = 3.1863
;DECam_r, R_r = 2.1401
;DECam_i, R_i = 1.5690
;DECam_z, R_z = 1.1957
;DECam_Y, R_Y = 1.0476
;For example:
;A_u = R_u * EBV_SFD98
ext = {u:3.9631,g:3.1863,r:2.1401,i:1.5690,z:1.1957,y:1.0476}

; Make the hess diagrams
dx = 0.05
dy = 0.1
xr = [-0.4, 1.6]
yr = [23.0,19.5]
;yr = [24.0,18.0]
smlen = 3
;nx = range(xr)/dx+1
;ny = range(yr)/dy+1

; Observed, dereddened
col = stars.g-stars.ebv*ext.g - (stars.i-stars.ebv*ext.i)
mag = stars.g-stars.ebv*ext.g
undefine,dum,obsim
hess,col,mag,dum,obsim0,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
obsim = gsmooth(obsim0,smlen)  ; smoothing
sz = size(obsim)
nx = sz[1]
ny = sz[2]

; Completeness
cdx = 0.20
cdy = 0.40
undefine,dum,cim
grec = where(ast.g gt 0 and ast.i gt 0 and ast.g lt 50 and ast.i lt 50 and ast.recovered eq 1,ngrec)
hess,ast.inp_g-ast.inp_i,ast.inp_g,dum,imall,dx=cdx,dy=cdy,xr=xr,yr=yr,/noplot
hess,ast[grec].inp_g-ast[grec].inp_i,ast[grec].inp_g,dum,imrec,dx=cdx,dy=cdy,xr=xr,yr=yr,/noplot
compim0 = imrec/(imall>1.0)
compim = congrid(compim0,nx,ny,/interp)

; Errors in g and i
gdg = where(stars.g lt 50,ng)
gcoef = robust_poly_fit(stars[gdg].g,alog10(stars[gdg].gerr),2)
gdi = where(stars.i lt 50,ni)
icoef = robust_poly_fit(stars[gdi].i,alog10(stars[gdi].ierr),2)
xx = scale_vector(findgen(100),14,26)

; Model, each component separately
; add uncertainties
gobs = model.g+randomn(seed,nmodel)*10.0^poly(model.g,gcoef)
iobs = model.i+randomn(seed,nmodel)*10.0^poly(model.i,icoef)
g = where(model.popid le 6,ng)  ; thin disk, 2-6
undefine,dum,mim1
hess,gobs[g]-iobs[g],gobs[g],dum,mim1,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
;hess,model[g].g-model[g].i,model[g].g,dum,mim1,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
g = where(model.popid eq 8,ng)  ; thick disk, 8
undefine,dum,mim2
hess,gobs[g]-iobs[g],gobs[g],dum,mim2,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
;hess,model[g].g-model[g].i,model[g].g,dum,mim2,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
g = where(model.popid eq 9,ng)  ; halo, 9
undefine,dum,mim3
hess,gobs[g]-iobs[g],gobs[g],dum,mim3,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
;hess,model[g].g-model[g].i,model[g].g,dum,mim3,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
sz = size(mim1)
mwmodel = fltarr(sz[1],sz[2],3)
mwmodel[*,*,0] = mim1
mwmodel[*,*,1] = mim2
mwmodel[*,*,2] = mim3

; 2D color and magnitude arrays
colim = xarr#replicate(1,sz[2])
magim = replicate(1,sz[1])#yarr


; Now fit the multi-component model
; 0 - thindisk_scale
; 1 - thickdisk_scale
; 2 - halo_scale
; 3 - lmcdist
; 4 - lmcspread
; 5 - lmcscale

parinfo = replicate({limited:[1,0],limits:[0.0,0.0],fixed:0},6)
initpars = [1.0,1.0,1.0,18.5,0.03,1.0]
xx = findgen(sz[1])#replicate(1,sz[2])
yy = obsim
err = sqrt(yy)>1.0

; Fit the MW components first
initpars[5] = 0
parinfo[3:5].fixed = 1
mask = float(colim gt 0.6 or magim lt 21.3)
yy = obsim*mask
fa = {mwmodel:mwmodel,lmcsynth:synth,compim:compim,dx:dx,dy:dy,xr:xr,yr:yr,gcoef:gcoef,icoef:icoef,smlen:smlen,mask:mask}
mvpars = MPFITFUN('lmc_hess_model',xx,yy,err,initpars,functargs=fa,parinfo=parinfo,dof=dof,$
                bestnorm=chisq,status=status,yfit=yfit)  ;,/quiet)

zmax = max([yy,yfit])
wset,0
displayc,yy,xarr,yarr,/yflip,min=0,max=zmax
wset,1
displayc,yfit,xarr,yarr,/yflip,min=0,max=zmax

im = lmc_hess_model(xx,mwpars,_extra=fa)

stop

end
