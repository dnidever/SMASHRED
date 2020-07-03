pro dr2_depth

; Measure 5sigma depth

files = file_search('/net/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/v6/4*combined_allobj.fits.gz',count=nfiles)
fstr = replicate({field:'',ra:0.0d0,dec:0.0d0,nobj:0L,depth:fltarr(5)+99.99},nfiles)

for i=0,nfiles-1 do begin
  fstr[i].field = file_basename(files[i],'_object.fits')
  obj = mrdfits(files[i],1,/silent)
  nobj = n_elements(obj)
  fstr[i].nobj = nobj
  if nobj lt 100 then goto,BOMB
  g = where(obj.depthflag eq 2 and obj.uerr gt 0.19 and obj.uerr lt 0.21,ng)
  if ng eq 0 then g = where(obj.uerr gt 0.19 and obj.uerr lt 0.21,ng)
  if ng gt 5 then fstr[i].depth[0] = median(obj[g].u)

  g = where(obj.depthflag eq 2 and obj.gerr gt 0.19 and obj.gerr lt 0.21,ng)
  if ng eq 0 then g = where(obj.gerr gt 0.19 and obj.gerr lt 0.21,ng)
  if ng gt 5 then fstr[i].depth[1] = median(obj[g].g)

  g = where(obj.depthflag eq 2 and obj.rerr gt 0.19 and obj.rerr lt 0.21,ng)
  if ng eq 0 then g = where(obj.rerr gt 0.19 and obj.rerr lt 0.21,ng)
  if ng gt 5 then fstr[i].depth[2] = median(obj[g].r)

  g = where(obj.depthflag eq 2 and obj.ierr gt 0.19 and obj.ierr lt 0.21,ng)
  if ng eq 0 then g = where(obj.ierr gt 0.19 and obj.ierr lt 0.21,ng)
  if ng gt 5 then fstr[i].depth[3] = median(obj[g].i)

  g = where(obj.depthflag eq 2 and obj.zerr gt 0.19 and obj.zerr lt 0.21,ng)
  if ng eq 0 then g = where(obj.zerr gt 0.19 and obj.zerr lt 0.21,ng)
  if ng gt 5 then fstr[i].depth[4] = median(obj[g].z)

  print,i+1,' ',fstr[i].field,' ',fstr[i].depth[0],' ',fstr[i].depth[1],' ',fstr[i].depth[2],' ',fstr[i].depth[3],' ',fstr[i].depth[4]
  BOMB:
  ;stop
endfor

radeg = 180.0d0 / !dpi
pix = long(strmid(fstr.field,0,5))
PIX2ANG_RING,64,pix,theta,phi
ra = phi*radeg
dec = 90-theta*radeg
fstr.ra = ra
fstr.dec = dec

;;MWRFITS,fstr,'/net/dl1/users/dnidever/smash/cp/red/photred/catalogs/final/v6/dr2_depth.fits',/create

stop

end
