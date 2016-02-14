pro smashred_getredinfo,info

; combine all the catalogs for a given field (short/long and multiple
; nights), use the overlap to do "ubercal", use APASS to do
; calibration of ~g-band, and use stellar locus regression (SLR) to
; fix the colors.
;
;-match and concatenate all of the long/short catalogs for all nights for a given field
;-use the overlapping short exposures to find a common zeropoint for each band
;-use SLR to calibrate the colors, maybe use APASS to set zeropoint for g-band??
;-use the "summary" files to figure out which field number it actually is

dir = '/data/smash/cp/red/photred/catalogs/'

; Load the SMASH fields file
smash = importascii(dir+'pro/smash_fields_final.txt',/header)

; Get all of the summary files
sumfiles = file_search(dir+'inst/20*/F*summary.fits',count=nsumfiles)
; Load the summary file information
info = replicate({file:'',object:'',field:'',sh:-1,night:'',nexp:0L,bands:'',fstr:ptr_new()},nsumfiles)
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
  info[i].night = file_basename(file_dirname(sumfiles[i]))
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
  print,i+1,base,info[i].field,info[i].night,info[i].nexp,info[i].bands,format='(I5,A13,A10,A12,I5,A8)'
endfor

end
