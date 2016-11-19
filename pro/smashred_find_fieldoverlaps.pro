pro smashred_find_fieldoverlaps,version

; Find field overlaps using their chips catalogs

if n_elements(version) eq 0 then begin
  print,'Please enter the version number, e.g. "v2"'
  return
endif
dir = '/data/smash/cp/red/photred/catalogs/final/'+version+'/'

filters = ['u','g','r','i','z']
nfilters = n_elements(filters)

chipfiles = file_search(dir+'Field*_combined_chips.fits.gz',count=nchipfiles)

nfields = nchipfiles
fieldstr = replicate({field:'',cenra:0.0d0,cendec:0.0d0,rar:dblarr(2),decr:dblarr(2),$
                      ndet:intarr(5),calib:intarr(5),noverlap:0,overlapfield:strarr(10),$
                      overlapindex:intarr(10)-1},nfields)

; Load all of the chip information
print,'Loading the chips files'
for i=0,nchipfiles-1 do begin
  chstr = mrdfits(chipfiles[i],1,/silent)
  field = file_basename(chipfiles[i],'_combined_chips.fits.gz')
  fieldstr[i].field = field
  fieldstr[i].rar = minmax(chstr.vertices_ra)
  fieldstr[i].decr = minmax(chstr.vertices_dec)
  fieldstr[i].cenra = mean(fieldstr[i].rar)
  fieldstr[i].cendec = mean(fieldstr[i].decr)
  ui = uniq(chstr.expnum,sort(chstr.expnum))
  uexpnum = chstr[ui].expnum
  ufilters = chstr[ui].filter
  ndet = intarr(5)
  for j=0,nfilters-1 do begin
    dum = where(ufilters eq filters[j],ndet1)
    ndet[j] = ndet1
    dum = where(chstr[ui].calibrated eq 1 and chstr[ui].photometric eq 1 and chstr[ui].badsoln eq 0 and chstr[ui].filter eq filters[j],ncalib)
    if ncalib gt 0 then fieldstr[i].calib[j]=1
  endfor
  fieldstr[i].ndet = ndet
endfor


; Now figure out the overlaps
for i=0,nfields-1 do begin
  print,strtrim(i+1,2),' ',fieldstr[i].field
  vertices_ra = [fieldstr[i].rar[0],fieldstr[i].rar[1],fieldstr[i].rar[1],fieldstr[i].rar[0]]
  vertices_dec = [fieldstr[i].decr[0],fieldstr[i].decr[0],fieldstr[i].decr[1],fieldstr[i].decr[1]]

  ; Find potential overlaps within the bounding rectangle
  ;  loop over the rest of the fields
  overlap = intarr(nfields)
  for j=0,nfields-1 do begin
    vertices_ra2 = [fieldstr[j].rar[0],fieldstr[j].rar[1],fieldstr[j].rar[1],fieldstr[j].rar[0]]
    vertices_dec2 = [fieldstr[j].decr[0],fieldstr[j].decr[0],fieldstr[j].decr[1],fieldstr[j].decr[1]]
    if j ne i then overlap[j] = dopolygonsoverlap(vertices_ra,vertices_dec,vertices_ra2,vertices_dec2)
  endfor

  ; Some overlaps, check further
  overlapind = where(overlap eq 1,noverlapind)
  print,strtrim(noverlapind,2),' potential overlap(s) found'
  if noverlapind gt 0 then begin
    brightfile1 = dir+fieldstr[i].field+'_combined_allobj_bright.fits'
    bright1 = mrdfits(brightfile1,1,/silent)
  endif
  for j=0,noverlapind-1 do begin
    ; Restore the bright file
    brightfile2 = dir+fieldstr[overlapind[j]].field+'_combined_allobj_bright.fits'
    bright2 = mrdfits(brightfile2,1,/silent)
    ; Crossmatch
    SRCMATCH,bright1.ra,bright1.dec,bright2.ra,bright2.dec,1.0,ind1,ind2,/sph,count=nmatch
    print,strtrim(nmatch,2),' cross-matches with ',fieldstr[overlapind[j]].field
    if nmatch gt 0 then begin
      fieldstr[i].overlapfield[fieldstr[i].noverlap] = fieldstr[overlapind[j]].field
      fieldstr[i].overlapindex[fieldstr[i].noverlap] = overlapind[j]
      fieldstr[i].noverlap += 1
    endif
  endfor

  ;stop

endfor  ; field loop

; Save the final structure
;MWRFITS,fieldstr,dir+'smash_fieldoverlaps.fits',/create

stop

end
