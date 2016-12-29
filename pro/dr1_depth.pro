pro dr1_depth

; Measure 5sigma depth

files = file_search('/data/smash/cp/red/photred/catalogs/final/v4/db/dr1/Field*object.fits',count=nfiles)
fstr = replicate({field:'',depth:fltarr(5)},nfiles)

for i=0,nfiles-1 do begin
  fstr[i].field = file_basename(files[i],'_object.fits')
  obj = mrdfits(files[i],1,/silent)
  g = where(obj.depthflag eq 2 and obj.uerr gt 0.19 and obj.uerr lt 0.21,ng)
  fstr[i].depth[0] = median(obj[g].umag)
  g = where(obj.depthflag eq 2 and obj.gerr gt 0.19 and obj.gerr lt 0.21,ng)
  fstr[i].depth[1] = median(obj[g].gmag)
  g = where(obj.depthflag eq 2 and obj.rerr gt 0.19 and obj.rerr lt 0.21,ng)
  fstr[i].depth[2] = median(obj[g].rmag)
  g = where(obj.depthflag eq 2 and obj.ierr gt 0.19 and obj.ierr lt 0.21,ng)
  fstr[i].depth[3] = median(obj[g].imag)
  g = where(obj.depthflag eq 2 and obj.zerr gt 0.19 and obj.zerr lt 0.21,ng)
  fstr[i].depth[4] = median(obj[g].zmag)
  print,i+1,' ',fstr[i].field,' ',fstr[i].depth[0],' ',fstr[i].depth[1],' ',fstr[i].depth[2],' ',fstr[i].depth[3],' ',fstr[i].depth[4]
  ;stop
endfor

stop

end
