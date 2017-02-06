pro stargalaxy_separation_all,field,version

;; Make "star" catalogs.

rootdir = smashred_rootdir()
if n_elements(version) eq 0 then version='v5'
catdir = rootdir+'cp/red/photred/catalogs/final/'+version+'/'

redfields = mrdfits(rootdir+'cp/red/photred/catalogs/pro/check_calibrated_v4.fits',1)
redfields.field = strtrim(redfields.field,2)
fields = redfields.field
nfields = n_elements(fields)

;; Loop through the fields
for i=0,nfields-1 do begin
  ifield = fields[i]
  print,strtrim(i+1,2),' ',ifield
   
  STARGALAXY_SEPARATION,ifield,version
  STARGALAXY_SEPARATION,ifield,version,/deep

  ;stop
  
endfor

stop

end
