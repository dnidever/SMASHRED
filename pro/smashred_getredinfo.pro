;+
;
; SMASHRED_GETREDINFO
;
; combine all the catalogs for a given field (short/long and multiple
; nights), use the overlap to do "ubercal", use APASS to do
; calibration of ~g-band, and use stellar locus regression (SLR) to
; fix the colors.
;
;-match and concatenate all of the long/short catalogs for all nights for a given field
;-use the overlapping short exposures to find a common zeropoint for each band
;-use SLR to calibrate the colors, maybe use APASS to set zeropoint for g-band??
;-use the "summary" files to figure out which field number it actually is
;
; INPUTS:
;  =dir       The directory with the catalogs.  By default,
;               dir='/data/smash/cp/red/photred/'
;  =sumfiles  Input list of summary files.  Otherwise
;               dir+'/20??????/*summaryf.its' is searched.
;  /silent    Don't print anything to the screen.
;
; OUTPUTS:
;  info     The structures with information for each SMASH PHOTRED catalog.
;
; USAGE:
;  IDL>smashred_getredinfo,info
;
; By D.Nidever Feb 2016
;-

pro smashred_getredinfo,info,dir=dir,sumfiles=sumfiles,silent=silent

undefine,info

rootdir = SMASHRED_ROOTDIR()
if n_elements(dir) eq 0 then dir=rootdir+'cp/red/photred/'
;if n_elements(dir) eq 0 then dir='/data/smash/cp/red/photred/'
;dir='/data/smash/cp/red/photred/catalogs/inst/'

; Directory not found
if file_test(dir,/directory) eq 0 then begin
  print,dir,' NOT FOUND'
  return
endif

; Load the SMASH fields file
;smash = importascii('/data/smash/cp/red/photred/catalogs/pro/smash_fields_final.txt',/header,/silent)
smash = importascii(rootdir+'cp/red/photred/catalogs/pro/smash_fields_final.txt',/header,/silent)

; Get all of the summary files
;sumfiles = file_search(dir+'inst/20*/F*summary.fits',count=nsumfiles)
;sumfiles = file_search(dir+'inst/20*/*summary.fits',count=nsumfiles)
;sumfiles = file_search(dir+'/20*/*summary.fits',count=nsumfiles)
if n_elements(sumfiles) eq 0 then begin
  sumfiles = file_search(dir+'/20??????/*summary.fits',count=nsumfiles)
  if not keyword_set(silent) then print,'Found ',strtrim(nsumfiles,2),' PHOTRED summary files in ',dir,'/20??????/'
endif else begin
  if not keyword_set(silent) then print,'Using ',strtrim(nsumfiles,2),' input PHOTRED summary files'
  nsumfiles = n_elements(sumfiles)
endelse
if nsumfiles eq 0 then begin
  print,'No PHOTRED summary files'
  return
endif

; Create the output structure
info = replicate({file:'',object:'',field:'',sh:-1,night:'',nexp:0L,bands:'',fstr:ptr_new()},nsumfiles)

; Load the summary file information
if not keyword_set(silent) then print,'   NUM     CATNAME SMASHFIELD     NIGHT   NEXP  FILTERS'
for i=0,nsumfiles-1 do begin
  fstr = mrdfits(sumfiles[i],1,/silent)
  base = file_basename(sumfiles[i],'_summary.fits')
  info[i].file = sumfiles[i]
  if strmid(base,1,/reverse) eq 'sh' then begin
    info[i].sh = 1
    len = strlen(base)
    info[i].object = strmid(base,0,len-2)
  endif else begin
    info[i].sh = 0
    info[i].object = base
  endelse
  info[i].night = strmid(file_basename(file_dirname(sumfiles[i])),0,8)  ; 20130320
  ; /datalab/users/dnidever/smash/cp/red/photred/deep/Field110/Field110_summary.fits'
  if stregex(sumfiles[i],'deep',/boolean) then info[i].night='deep'
  info[i].nexp = n_elements(fstr)
  ui = uniq(fstr.filter,sort(fstr.filter))
  info[i].bands = strjoin(fstr[ui].filter)

  ; Get the "standard" SMASH field number
  mnra = median([fstr.ra])
  mndec = median([fstr.dec])
  dist = sphdist(mnra,mndec,smash.radeg,smash.dedeg,/deg)
  bestind = first_el(minloc(dist))
  if dist[bestind] lt 1.2 then begin
    smashfield = 'Field'+strtrim(smash[bestind].num,2)
    info[i].field = smashfield
    ;print,'This is SMASH field: ',smashfield
  endif else begin
    ;print,'NO SMASH field MATCH. Using observed field name: ',ofield
    ;smashfield = ofield
  endelse

  info[i].fstr = ptr_new(fstr)
  if not keyword_set(silent) then $
    print,i+1,base,info[i].field,info[i].night,info[i].nexp,info[i].bands,format='(I5,A13,A10,A12,I5,A8)'
endfor

end
