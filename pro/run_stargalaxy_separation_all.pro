pro run_stargalaxy_separation_all,field,version

;; Make "star" catalogs.

rootdir = smashred_rootdir()
if n_elements(version) eq 0 then version='v5'
catdir = rootdir+'cp/red/photred/catalogs/final/'+version+'/'

redfields = mrdfits(rootdir+'cp/red/photred/catalogs/pro/check_calibrated_v4.fits',1)
redfields.field = strtrim(redfields.field,2)
fields = redfields.field
nfields = n_elements(fields)

; Field118 crashed in make_stellar_locus.pro
; Field183 crashed because not all bands observed
; Field184 crashed as well
; Field185, 186, same
; AHH, they have only depthflag=1

;; Loop through the fields
for i=74,nfields-1 do begin
;for i=0,nfields-1 do begin
  ifield = fields[i]
  print,strtrim(i+1,2),' ',ifield
   
  RUN_STARGALAXY_SEPARATION_SINGLE,ifield,version
  RUN_STARGALAXY_SEPARATION_SINGLE,ifield,version,/deep,/doast

  ;stop
  
endfor

stop

end
