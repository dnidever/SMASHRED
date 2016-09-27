pro measure_gaiasmash_colorterms,save=save

; Get calibrated fields with full ugriz
; and not in the LMC/SMC main bodies
restore,'/data/smash/cp/red/photred/catalogs/pro/check_calibrated.dat'
gd = where(info.calflag eq 2 and info.nuchips gt 0 and long(strmid(info.field,5)) gt 60,ngd)
fields = info[gd].field
nfields = n_elements(fields)

catdir = '/data/smash/cp/red/photred/catalogs/final/'
gaiadir = '/data/smash/cp/red/photred/gaia/'
;fields = ['Field101','Field110','Field134']
;nfields = n_elements(fields)

restore,'/data/smash/cp/red/photred/gaia/gaiasmash_matchedcats.dat'
goto,fitcoef

undefine,mobj,mgaia
fieldstr = replicate({field:'',fieldid:0,nmatch:0L,ucoef:dblarr(5),urms:99.99,ubin:fltarr(31),gcoef:dblarr(4),grms:99.99,$
                      gbin:fltarr(31),rcoef:dblarr(4),rrms:99.99,rbin:fltarr(31),icoef:dblarr(4),irms:99.99,ibin:fltarr(31),$
                      zcoef:dblarr(5),zrms:99.99,zbin:fltarr(31)},nfields)
for i=0,nfields-1 do begin
  print,strtrim(i+1,2),'/',strtrim(nfields,2),' ',fields[i]
  obj = mrdfits(catdir+fields[i]+'_combined_allobj.fits.gz',1)
  gaia = mrdfits(gaiadir+fields[i]+'_gaia.fits',1)

  srcmatch,obj.ra,obj.dec,gaia.ra_icrs,gaia.de_icrs,5.0,ind1,ind2,/sph,count=nmatch
  print,strtrim(nmatch,2),' matches'
  obj1 = obj[ind1]
  gaia1 = gaia[ind2]

  ; Convert to common obj schema without indices
  ;  gaia too
  schema_obj = {id:'',fieldid:0,ra:0.d0,dec:0.0d0,u:0.0,uerr:0.0,g:0.0,gerr:0.0,r:0.0,rerr:0.0,$
            i:0.0,ierr:0.0,z:0.0,zerr:0.0,chi:0.0,sharp:0.0,flag:0,prob:0.0,ebv:0.0}
  newobj1 = replicate(schema_obj,nmatch)
  struct_assign,obj1,newobj1
  newobj1.fieldid = long(strmid(fields[i],5))
  schema_gaia = {source:0L,ra_icrs:0.0d0,e_ra_icrs:0.0,de_icrs:0.0,e_de_icrs:0.0,_fg_:0.0d0,$
                 e__fg_:0.0d0,_gmag_:0.0}
  newgaia1 = replicate(schema_gaia,nmatch)
  struct_assign,gaia1,newgaia1

  plotc,obj1.g-obj1.i,obj1.g-gaia1._gmag_,gaia1._gmag_,ps=3,xr=[0,4],yr=[-1,3],xs=1,ys=1,xtit='g-i',ytit='g-G',tit=fields[i]

  fieldstr[i].field = fields[i]
  fieldstr[i].fieldid = long(strmid(fields[i],5))
  fieldstr[i].nmatch = nmatch
  ucoef = robust_poly_fitq(obj1.g-obj1.i,obj1.u-gaia1._gmag_,4)
  urms = mad(obj1.u-gaia1._gmag_-poly(obj1.g-obj1.i,ucoef))
  bindata,obj1.g-obj1.i,obj1.u-gaia1._gmag_,xbin,ybin,bin=0.1,/med,min=0.0,max=3.0
  fieldstr[i].ucoef = ucoef
  fieldstr[i].urms = urms
  fieldstr[i].ubin = ybin
  gcoef = robust_poly_fitq(obj1.g-obj1.i,obj1.g-gaia1._gmag_,3)
  grms = mad(obj1.g-gaia1._gmag_-poly(obj1.g-obj1.i,gcoef))
  bindata,obj1.g-obj1.i,obj1.g-gaia1._gmag_,xbin,ybin,bin=0.1,/med,min=0.0,max=3.0
  fieldstr[i].gcoef = gcoef
  fieldstr[i].grms = grms
  fieldstr[i].gbin = ybin
  rcoef = robust_poly_fitq(obj1.g-obj1.i,obj1.r-gaia1._gmag_,3)
  rrms = mad(obj1.r-gaia1._gmag_-poly(obj1.g-obj1.i,rcoef))
  bindata,obj1.g-obj1.i,obj1.r-gaia1._gmag_,xbin,ybin,bin=0.1,/med,min=0.0,max=3.0
  fieldstr[i].rcoef = rcoef
  fieldstr[i].rrms = rrms
  fieldstr[i].rbin = ybin
  icoef = robust_poly_fitq(obj1.g-obj1.i,obj1.i-gaia1._gmag_,3)
  irms = mad(obj1.i-gaia1._gmag_-poly(obj1.g-obj1.i,icoef))
  bindata,obj1.g-obj1.i,obj1.i-gaia1._gmag_,xbin,ybin,bin=0.1,/med,min=0.0,max=3.0
  fieldstr[i].icoef = icoef
  fieldstr[i].irms = irms
  fieldstr[i].ibin = ybin
  zcoef = robust_poly_fitq(obj1.g-obj1.i,obj1.z-gaia1._gmag_,4)
  zrms = mad(obj1.z-gaia1._gmag_-poly(obj1.g-obj1.i,zcoef))
  bindata,obj1.g-obj1.i,obj1.z-gaia1._gmag_,xbin,ybin,bin=0.1,/med,min=0.0,max=3.0
  fieldstr[i].zcoef = zcoef
  fieldstr[i].zrms = zrms
  fieldstr[i].zbin = ybin

  push,mobj,newobj1
  push,mgaia,newgaia1

  ;stop
endfor

; RMS in the binned values
urms = mad(fieldstr.ubin,dim=2)  ; ~5%
grms = mad(fieldstr.gbin,dim=2)  ; ~1%
rrms = mad(fieldstr.rbin,dim=2)  ; ~1%
irms = mad(fieldstr.ibin,dim=2)  ; ~1%
zrms = mad(fieldstr.zbin,dim=2)  ; ~1.5%

;save,fieldstr,mobj,mgaia,file='/data/smash/cp/red/photred/gaia/gaiasmash_matchedcats.dat'
stop

FITCOEF:

; Field132 shows higher scatter and multiple "sequences"
;  gaia or smash issue?
; Field171 as well
setdisp
!p.font = 0
nobj = n_elements(mobj)
otags = tag_names(mobj)
ftags = tag_names(fieldstr)
filters = ['u','g','r','i','z']
colr = mobj.g-mobj.i
x = scale_vector(findgen(100),0.0,3.0)

gaiastr = replicate({filter:'',bestcoef:dblarr(4),bestgirange:fltarr(2),decentcoef:dblarr(5),decentgirange:fltarr(2)},5)
gaiastr.filter = filters

; Filter loop
for i=0,4 do begin
  filt = filters[i]

  if keyword_set(save) then begin
    file = '/data/smash/cp/red/photred/gaia/plots/gaiasmash_colorterms_'+filt
    ps_open,file,/color,thick=4,/encap
    device,/inches,xsize=8.5,ysize=9.5
  endif

  ; Get the decent/best color ranges
  case filt of
  'u': begin
    ; Decent fit for 0.30 < g-i < 2.7
    gaiastr[i].decentgirange = [0.30,2.7]
    ; Best fit for 0.60< g-i < 2.2
    gaiastr[i].bestgirange = [0.60,2.2]
    yr = [-1.0,6.0]
  end
  'g': begin
    ; Becent fit for 0.0 < g-i < 2.8
    gaiastr[i].decentgirange = [0.0,2.8]
    ; Best fit for 0.50 < g-i < 1.7
    gaiastr[i].bestgirange = [0.5,1.7]
    yr = [-0.5,3.0]
  end
  'r': begin
    ; Decent fit for 0.25 < g-i < 2.60
    gaiastr[i].decentgirange = [0.25,2.60]
    ; Best fit for 0.50 < g-i < 2.0
    gaiastr[i].bestgirange = [0.5,2.0]
    yr = [-1.0,1.5]
  end
  'i': begin
    ; Decent fit for 0.0 < g-i < 2.8
    gaiastr[i].decentgirange = [0.0,2.8]
    ; Best fit for  0.40 < g-i < 1.9 
    gaiastr[i].bestgirange = [0.4,1.9]
    yr = [-1.5,1.0]
  end
  'z': begin
    ; Decent fit for 0.15 < g-i < 2.9
    gaiastr[i].decentgirange = [0.15,2.9]
    ; Best fit for 0.50 < g-i < 2.3
    gaiastr[i].bestgirange = [0.5,2.3]
    yr = [-1.5,0.5]
  end
  else: stop,' not an option'
  endcase

  ; Get good fields and their stars
  rmsind = where(ftags eq strupcase(filt)+'RMS',nrmsind)
  gdfields = where(fieldstr.(rmsind) lt 0.5,ngdfields)
  goodmask = intarr(nobj)
  for j=0,ngdfields-1 do begin
    MATCH,mobj.fieldid,fieldstr[gdfields[j]].fieldid,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then goodmask[ind1]=1
  endfor
  magind = where(otags eq strupcase(filt),nmagind)
  ; Get good stars to use
  gdobj = where(goodmask eq 1 and mobj.g lt 21 and mobj.i lt 21 and mgaia._gmag_ lt 21 and mobj.(magind) lt 21,ngdobj)

  ; Hess diagram
  hess,mobj[gdobj].g-mobj[gdobj].i,mobj[gdobj].(magind)-mgaia[gdobj]._gmag_,dx=0.01,dy=0.01,xr=[0,3.2],yr=yr,$
       xtit='g-i',ytit=filt+'-G',tit=filt+'-G',/log

  ; Overplot the individual fits
  coefind = where(ftags eq strupcase(filt)+'COEF',ncoefind)
  ;for j=0,ngdfields-1 do oplot,x,poly(x,fieldstr[gdfields[j]].(coefind)),co=255
  ; Get values for individual fits and calculate RMS
  fit = fltarr(100,ngdfields)
  for j=0,ngdfields-1 do fit[*,j]=poly(x,fieldstr[gdfields[j]].(coefind))
  ;plot,x,mad(fit,dim=2),ps=-1

  ; Poly fit to all stars
  coef_all = robust_poly_fitq(colr[gdobj],mobj[gdobj].(magind)-mgaia[gdobj]._gmag_,3)
  ;oplot,x,poly(x,coef_all),co=250

  ; Bin the data and poly fitting
  bindata,colr[gdobj],mobj[gdobj].(magind)-mgaia[gdobj]._gmag_,xbin,ybin,binsize=0.1,/med,min=0.0,max=3.0,gdind=gdind
  ;oplot,xbin,ybin,ps=-1,sym=3,co=80
  coef_bin_all = robust_poly_fitq(xbin[gdind],ybin[gdind],4)
  ;oplot,x,poly(x,coef_bin_all),co=150

  ; Get decent fit
  decentind = where(xbin ge gaiastr[i].decentgirange[0] and xbin le gaiastr[i].decentgirange[1] and $
                    finite(ybin) eq 1,ndecentind)
  coef_bin_decent = robust_poly_fitq(xbin[decentind],ybin[decentind],4)
  gaiastr[i].decentcoef = coef_bin_decent
  oplot,xbin[decentind],poly(xbin[decentind],coef_bin_decent),co=255

  ; Get best fit
  bestind = where(xbin ge gaiastr[i].bestgirange[0] and xbin le gaiastr[i].bestgirange[1] and $
                  finite(ybin) eq 1,nbestind)
  coef_bin_best = robust_poly_fitq(xbin[bestind],ybin[bestind],3)
  gaiastr[i].bestcoef[0:3] = coef_bin_best
  oplot,xbin[bestind],poly(xbin[bestind],coef_bin_best),co=250,linestyle=2

  ; Legend
  ;al_legend,['Binned','Decent','Best'],textcolor=[80,200,250],psym=[1,0,0],charsize=1.5,/top,/left
  al_legend,['Decent','Best'],textcolor=[255,250],charsize=1.5,/top,/left

  if keyword_set(save) then begin
    ps_close
    ps2png,file+'.eps',/eps
  endif

  ;stop

endfor

; MWRFITS,gaiastr,gaiadir+'gaiasmash_colorterms.fits',/create

stop

end
