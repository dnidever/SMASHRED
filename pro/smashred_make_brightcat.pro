;+
;
; SMASHRED_MAKE_BRIGHTCAT
;
; Make small allobj catalogs of bright stars
; for help in crossmatching across overlapping
; fields.
;
; INPUTS:
;  field    The name of the field, e.g. Field130.
;  =dir     The directory where the catalogs live.
;  =maglim  The magnitude limit to impose.  The
;             default is 22.
;  /redo    Redo this field.
;
; OUTPUTS:
;  A small version of the allobj catalog with
;  similar columns and only bright stars is
;  saved the name FIELD_combined_allobj_bright.fits
;  Stars are selected with detections in all
;  observed bands, good sharp/chi values, and
;  are brighter than MAGLIM in the reference band
;  (first band in g/r/i/z with some detections).
;
; USAGE:
;  IDL>smashred_make_brightcat,'Field130'
;
; By D.Nidever  Oct 2016
;-

pro smashred_make_brightcat,field,redo=redo,dir=dir,maglim=maglim

; Not enough inputs
if n_elements(field) eq 0 then begin
  print,'Syntax smashred_make_brightcat,field,redo=redo,dir=dir,maglim=maglim'
  return
endif

; Defaults
if n_elements(dir) eq 0 then dir = '/data/smash/cp/red/photred/catalogs/gaiacal/'
if n_elements(maglim) eq 0 then maglim=22.0

; Check that the allobj catalog exists
allobjfile = dir+field+'_combined_allobj.fits.gz'
if file_test(allobjfile) eq 0 then begin
  print,allobjfile,' NOT FOUND'
  return
endif

outfile = dir+field+'_combined_allobj_bright.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,outfile,' EXISTS and /redo NOT set.'
  return
endif

; Load the file
print,'Loading ',allobjfile
allobj = mrdfits(allobjfile,1)
nallobj = n_elements(allobj)
tags = tag_names(allobj)

; Make a mask to filter out the "good" bright ones
;  want detections in all bands that were observed
mask = intarr(nallobj)+1
filters = ['u','g','r','i','z']
nfilters = n_elements(filters)
numgood = lonarr(nfilters)
for i=0,nfilters-1 do begin
  magind = where(tags eq strupcase(filters[i]),nmagind)
  gdmag = where(allobj.(magind) lt 50,ngdmag)
  numgood[i] = ngdmag
  if ngdmag gt 0 then mask=mask and (allobj.(magind) lt 50)
endfor
mask = mask and (abs(allobj.sharp) lt 1 and allobj.chi lt 4)
; Magnitude limit
;  g,r,i,z, first one with some good measurements
reffiltind = first_el(where(filters ne 'u' and numgood gt 0))
refmagind = where(tags eq strupcase(filters[reffiltind]),nrefmagind)
gd = where(mask eq 1 and (allobj.(refmagind) le maglim),ngd)
print,strtrim(ngd,2),' bright stars'

; Make new output structure
schema_new = {id:'',ra:0.0d0,dec:0.0d0,u:0.0,uerr:0.0,g:0.0,gerr:0.0,r:0.0,rerr:0.0,$
              i:0.0,ierr:0.0,z:0.0,zerr:0.0,chi:0.0,sharp:0.0,prob:0.0,ebv:0.0}
newobj = replicate(schema_new,ngd)
STRUCT_ASSIGN,allobj[gd],newobj

; Write out the bright file
print,'Writing bright allobj catalog to ',outfile
MWRFITS,newobj,outfile,/create

;stop

end
