pro measure_gaiasmash_colorterms

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

undefine,mobj,mgaia
fieldstr = replicate({field:'',fieldid:0,nmatch:0L,ucoef:dblarr(5),urms:99.99,gcoef:dblarr(4),grms:99.99,$
                      rcoef:dblarr(4),rrms:99.99,icoef:dblarr(4),irms:99.99,zcoef:dblarr(5),zrms:99.99},nfields)
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
  fieldstr[i].ucoef = ucoef
  fieldstr[i].urms = urms
  gcoef = robust_poly_fitq(obj1.g-obj1.i,obj1.g-gaia1._gmag_,3)
  grms = mad(obj1.g-gaia1._gmag_-poly(obj1.g-obj1.i,gcoef))
  fieldstr[i].gcoef = gcoef
  fieldstr[i].grms = grms
  rcoef = robust_poly_fitq(obj1.g-obj1.i,obj1.r-gaia1._gmag_,3)
  rrms = mad(obj1.r-gaia1._gmag_-poly(obj1.g-obj1.i,rcoef))
  fieldstr[i].rcoef = rcoef
  fieldstr[i].rrms = rrms
  icoef = robust_poly_fitq(obj1.g-obj1.i,obj1.i-gaia1._gmag_,3)
  irms = mad(obj1.i-gaia1._gmag_-poly(obj1.g-obj1.i,icoef))
  fieldstr[i].icoef = icoef
  fieldstr[i].irms = irms
  zcoef = robust_poly_fitq(obj1.g-obj1.i,obj1.z-gaia1._gmag_,4)
  zrms = mad(obj1.z-gaia1._gmag_-poly(obj1.g-obj1.i,zcoef))
  fieldstr[i].zcoef = zcoef
  fieldstr[i].zrms = zrms

  push,mobj,newobj1
  push,mgaia,newgaia1

  ;stop
endfor

;save,fieldstr,mobj,mgaia,file='/data/smash/cp/red/photred/gaia/gaiasmash_matchedcats.dat'
stop

; Field132 shows higher scatter and multiple "sequences"
;  gaia or smash issue?
; Field171 as well

x = scale_vector(findgen(100),-1,4)

;plot,mobj.g-mobj.i,mobj.u-mgaia._gmag_,ps=3,xr=[0,4],yr=[-2,6],xs=1,ys=1,xtit='g-i',ytit='u-G',tit='u-mag'
colr = mobj.g-mobj.i
gdu = where(mobj.g lt 21 and mobj.i lt 21 and mgaia._gmag_ lt 21 and mobj.u lt 21,ngdu)
hess,mobj[gdu].g-mobj[gdu].i,mobj[gdu].u-mgaia[gdu]._gmag_,dx=0.01,dy=0.01,xr=[0,4],yr=[-2,6],xtit='g-i',ytit='u-G',tit='u-G',/log
bindata,colr[gdu],mobj[gdu].u-mgaia[gdu]._gmag_,xbin,ybin,binsize=0.1,/med,min=0.2,max=3.0
ucoef = robust_poly_fitq(colr[gdu],mobj[gdu].u-mgaia[gdu]._gmag_,3)
oplot,x,poly(x,ucoef),co=250

;plot,mobj.g-mobj.i,mobj.g-mgaia._gmag_,ps=3,xr=[0,4],yr=[-1,3],xs=1,ys=1,xtit='g-i',ytit='g-G',tit='g-mag'
gdg = where(mobj.g lt 21 and mobj.i lt 21 and mgaia._gmag_ lt 21,ngdg)
colr = mobj.g-mobj.i
hess,mobj[gdg].g-mobj[gdg].i,mobj[gdg].g-mgaia[gdg]._gmag_,dx=0.01,dy=0.01,xr=[0,4],yr=[-1,3],xtit='g-i',ytit='g-G',tit='g-G',/log
bindata,colr[gdg],mobj[gdg].g-mgaia[gdg]._gmag_,xbin,ybin,binsize=0.1,/med,min=0.2,max=3.0
gcoef = robust_poly_fitq(colr[gdg],mobj[gdg].g-mgaia[gdg]._gmag_,3)
oplot,x,poly(x,gcoef),co=250

;plot,mobj.g-mobj.i,mobj.r-mgaia._gmag_,ps=3,xr=[0,4],yr=[-2,2],xs=1,ys=1,xtit='g-i',ytit='r-G',tit='r-mag'
colr = mobj.g-mobj.i
gdr = where(mobj.g lt 21 and mobj.i lt 21 and mgaia._gmag_ lt 21 and mobj.r lt 21,ngdr)
hess,mobj[gdr].g-mobj[gdr].i,mobj[gdr].r-mgaia[gdr]._gmag_,dx=0.01,dy=0.01,xr=[0,4],yr=[-1,2],xtit='g-i',ytit='r-G',tit='r-G',/log
bindata,colr[gdr],mobj[gdr].r-mgaia[gdr]._gmag_,xbin,ybin,binsize=0.1,/med,min=0.2,max=3.0
rcoef = robust_poly_fitq(colr[gdr],mobj[gdr].r-mgaia[gdr]._gmag_,3)
oplot,x,poly(x,rcoef),co=250

;plot,mobj.g-mobj.i,mobj.i-mgaia._gmag_,ps=3,xr=[0,4],yr=[-2,2],xs=1,ys=1,xtit='g-i',ytit='i-G',tit='i-mag'
colr = mobj.g-mobj.i
gdi = where(mobj.g lt 21 and mobj.i lt 21 and mgaia._gmag_ lt 21,ngdi)
hess,mobj[gdi].g-mobj[gdi].i,mobj[gdi].i-mgaia[gdi]._gmag_,dx=0.01,dy=0.01,xr=[0,4],yr=[-1.5,1],xtit='g-i',ytit='i-G',tit='i-G',/log
bindata,colr[gdi],mobj[gdi].i-mgaia[gdi]._gmag_,xbin,ybin,binsize=0.1,/med,min=0.2,max=3.0
icoef = robust_poly_fitq(colr[gdi],mobj[gdi].i-mgaia[gdi]._gmag_,3)
oplot,x,poly(x,icoef),co=250

;plot,mobj.g-mobj.i,mobj.z-mgaia._gmag_,ps=3,xr=[0,4],yr=[-2,2],xs=1,ys=1,xtit='g-i',ytit='z-G',tit='z-mag'
colr = mobj.g-mobj.i
gdz = where(mobj.g lt 21 and mobj.i lt 21 and mgaia._gmag_ lt 21 and mobj.z lt 21,ngdz)
hess,mobj[gdz].g-mobj[gdz].i,mobj[gdz].z-mgaia[gdz]._gmag_,dx=0.01,dy=0.01,xr=[0,4],yr=[-2,1],xtit='g-i',ytit='z-G',tit='z-G',/log
bindata,colr[gdz],mobj[gdz].z-mgaia[gdz]._gmag_,xbin,ybin,binsize=0.1,/med,min=0.2,max=3.0
zcoef = robust_poly_fitq(colr[gdz],mobj[gdz].z-mgaia[gdz]._gmag_,3)
oplot,x,poly(x,zcoef),co=250

stop

end
