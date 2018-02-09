function lmc_hess_model,x,par,mwmodel=mwmodel,lmcsynth=lmcsynth,compim=compim,dx=dx,dy=dy,xr=xr,yr=yr

; Construct the model
; thindisk_scale, thickdisk_scale, halo_scale lmcdist, lmcscale
mwim = mwmodel[*,*,0] * par[0]  ; thin disk
mwim += mwmodel[*,*,1] * par[1] ; thick disk
mwim += mwmodel[*,*,2] * par[2] ; halo

undefine,dum,lmcim
hess,lmcsynth.gobs-lmcsynth.iobs,lmcsynth.gobs+par[3],dum,lmcim,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
lmcim *= par[4]

im = (mwim + lmcim) * compim

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
gd = where(stars.g lt 50 and stars.i lt 50,ngd)
stars = stars[gd]
; completeness
ast = mrdfits(catdir+'stars1/'+field+'_complete_stars.fits.gz',1)

; Restore Thomas' model
model = importascii(rootdir+'mwmodels/dartmouth/CMD_'+strlowcase(field),/header)

; Restore the LMC synthetic isochrone
synth = mrdfits(rootdir+'lmchalo/parsec_10gyr_-1.00_synth.fits.gz',1)
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
dx = 0.02
dy = 0.05
xr = [-1.0, 3.5]
yr = [25.0,14.0]

; Observed, dereddened
col = stars.g-stars.ebv*ext.g - (stars.i-stars.ebv*ext.i)
mag = stars.g-stars.ebv*ext.g
undefine,dum,obsim
hess,col,mag,dum,obsim,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot

; Model, each component separately
g = where(model.popid le 6,ng)  ; thin disk, 2-6
undefine,dum,mim1
hess,model[g].g-model[g].i,model[g].g,dum,mim1,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
g = where(model.popid eq 8,ng)  ; thick disk, 8
undefine,dum,mim2
hess,model[g].g-model[g].i,model[g].g,dum,mim2,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
g = where(model.popid eq 9,ng)  ; halo, 9
undefine,dum,mim3
hess,model[g].g-model[g].i,model[g].g,dum,mim3,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
sz = size(mim1)
mwmodel = fltarr(sz[1],sz[2],3)
mwmodel[*,*,0] = mim1
mwmodel[*,*,1] = mim2
mwmodel[*,*,2] = mim3

; Completeness
undefine,dum,cim
grec = where(ast.g lt 50 and ast.i lt 50,ngrec)
hess,ast.inp_g-ast.inp_i,ast.inp_g,dum,imall,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
hess,ast[grec].inp_g-ast[grec].inp_i,ast[grec].inp_g,dum,imrec,dx=dx,dy=dy,xr=xr,yr=yr,/noplot
compim = imrec/(imall>1.0)


; Now fit the multi-component model
; thindisk_scale, thickdisk_scale, halo_scale lmcdist, lmcscale

parinfo = replicate({limited:[1,0],limits:[0.0,0.0],fixed:0},5)
initpars = [1.0,1.0,1.0,18.5,1.0]
xx = findgen(sz[1])#replicate(1,sz[2])
yy = obsim
err = sqrt(yy)>1.0
fa = {mwmodel:mwmodel,lmcsynth:synth,compim:compim,dx:dx,dy:dy,xr:xr,yr:yr}
pars = MPFITFUN('lmc_hess_model',xx,yy,err,initpars,functargs=fa,parinfo=parinfo,dof=dof,$
                bestnorm=chisq,status=status,yfit=yfit)  ;,/quiet)

; NEED TO ADD BETTER ERRORS!!


stop

end
