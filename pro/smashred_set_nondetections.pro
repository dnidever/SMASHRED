;+
;
; SMASHRED_SET_NONDETECTIONS
;
; Set the magnitudes for non-detections to 99.99 and to NAN
; for bad/missing/no data.
;
; INPUTS:
;  field   The SMASH field name, e.g. "Field100".
;  allobj  The structure with the average values per object.
;  =dir    The directory for the exposure map.
;
; OUTPUTS:
;  allobj  The magnitude columns are updated for non-detections.
;
; USAGE:
;  IDL>smashred_set_nondetections,field,allobj
;
; By D. Nidever  Sep. 2016
;-

pro smashred_set_nondetections,field,allobj,dir=dir

if n_elements(dir) eq 0 then dir = '/data/smash/cp/red/photred/catalogs/final/'

; Not enough inputs
if n_elements(field) eq 0 or n_elements(allobj) eq 0 then begin
  print,'Syntax - smashred_set_nondetections,field,allobj'
  return
endif

; Exposure map file
expmapfile = dir+field+'_combined_expmap.fits'
if file_test(expmapfile) eq 0 and file_test(expmapfile+'.gz') eq 0 then begin
  print,expmapfile,' NOT FOUND'
  return
endif
if file_test(expmapfile) eq 0 then expmapfile+='.gz'

tags = tag_names(allobj)

filters = ['u','g','r','i','z']
nfilters = n_elements(filters)

; Filter loop
FOR f=0,nfilters-1 do begin

  ; Load the exposure map file
  ;FITS_READ,expmapfile,expmap,exphead,exten=f
  ; fits_read has problems reading large gzip files
  expmap = MRDFITS(expmapfile,f,exphead)

  ; Deal with "bad" mags
  magind = where(tags eq strupcase(filters[f]),nmagind)
  bdobj = where(allobj.(magind) gt 50,nbdobj)
  if nbdobj gt 0 then begin
    ; Convert to exposure map coordinates
    HEAD_ADXY,exphead,allobj[bdobj].ra,allobj[bdobj].dec,xx,yy,/deg
    expval = expmap[round(xx),round(yy)]
    
    ; Set all bad to start
    allobj[bdobj].(magind) = !values.f_nan

    ; Set non-detections
    nondet = where(expval gt 0,n_nondet)
    if n_nondet gt 0 then allobj[bdobj[nondet]].(magind) = 99.99
    print,filters[f],' ',strtrim(n_nondet,2),' non-detections','  ',strtrim(nbdobj-n_nondet,2),' bad/no data'
  endif
ENDFOR

end
