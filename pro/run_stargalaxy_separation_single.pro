pro run_stargalaxy_separation_single,field,version,sversion=sversion,deep=deep,doast=doast,redo=redo

;; Make "star" catalogs.

if n_elements(field) eq 0 then begin
  print,'Syntax - run_stargalaxy_separation_single,field,version,sversion=sversion,deep=deep,doast=doast,redo=redo'
  return
endif

rootdir = smashred_rootdir()
if n_elements(version) eq 0 then version='v6'
if n_elements(sversion) eq 0 then sversion='1'
catdir = rootdir+'cp/red/photred/catalogs/final/'+version+'/'
if keyword_set(deep) then tag='_deep' else tag=''

; Does the output file already exist
outfile = catdir+'stars'+sversion+'/'+field+'_allobj'+tag+'_stars.fits'
if ((file_info(outfile)).exists eq 1 or (file_info(outfile+'.gz')).exists eq 1) and not keyword_set(redo) then begin
  print,outfile,' exists and /redo not set'
  return
endif

print,'Performing star/galaxy separation for field ',field
if keyword_set(deep) then print,'Using "deep-only" catalog'

;; Load the allobj file
file = catdir+field+'_combined_allobj'+tag+'.fits.gz'
if file_test(file) eq 0 then begin
  print,file+' NOT FOUND'
  return
endif
obj = mrdfits(file,1,/silent)

; Doing ASTs
if keyword_set(deep) and keyword_set(doast) then begin
  astfile = catdir+'ast/'+field+'_complete.fits.gz'
  if file_test(astfile) eq 1 then begin
    print,'Running AST catalog through star/galaxy separation code'
    astobj = MRDFITS(astfile,1,/silent)
  endif else begin
    print,astfile,' NOT FOUND'
    undefine,astobj
  endelse
endif

STARGALAXY_SEPARATION,obj,ind,astobj=astobj,astind=astind
obj2 = obj[ind]

; Write out pruned catalogs
print,'Writing out pruned catalogs'
if file_test(catdir+'stars'+sversion,/directory) eq 0 then file_mkdir,catdir+'stars'+sversion
file_delete,[outfile,outfile+'.gz'],/allow
MWRFITS,obj2,outfile,/create
spawn,['gzip',outfile],out,errout,/noshell   ; compress

; Write out pruned AST catalog
if keyword_set(deep) and keyword_set(doast) and n_elements(astobj) gt 0 then begin
  ;astobj2 = astobj[astind]
  ; The ones that passed all of the cuts are the final "recovered"
  ; artificial stars, leave the rest in there.
  recovered = lonarr(n_elements(astobj))
  recovered[astind] = 1
  astobj2 = astobj
  astobj2.recovered = recovered
  outastfile = catdir+'stars'+sversion+'/'+field+'_complete_stars.fits'
  MWRFITS,astobj2,outastfile,/create
  if file_test(outastfile+'.gz') eq 1 then file_delete,outastfile+'.gz',/allow
  spawn,['gzip',outastfile],out,errout,/noshell   ; compress
endif

;stop

end
