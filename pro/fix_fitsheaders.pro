pro fix_fitsheaders,files,night=night

; the "BEGIN " lines are causing problems for daophot
; add a "COMMENT " before that

if n_elements(files) eq 0 and n_elements(night) eq 0 then begin
  print,'Syntax - fix_fitsheaders,files,night=night'
  return
endif

dir = '/data/smash/cp/red/photred/'

if n_elements(files) eq 0 then begin
  files = file_search(dir+'/'+night[0]+'F*/F*-*_??.fits',count=nfiles)
  if nfiles eq 0 then begin
    print,'NO fits files found in ',dir,'/',night[0]
    return
  endif
endif
nfiles = n_elements(files)
print,strtrim(nfiles,2),' files to fix.'

for i=0,nfiles-1 do begin
  head = headfits(files[i])
  bd1 = where(strmid(head,0,6) eq 'BEGIN ',nbd1)
  bd2 = where(strmid(head,0,6) eq 'EXTVER',nbd2)
  ; There's something to fix
  if nbd1 gt 0 or nbd2 gt 0 then begin
    FITS_READ,files[i],im,head
    if nbd1 gt 0 then head[bd1] = 'COMMENT '+head[bd1]  ; fix BEGIN line
    if nbd2 gt 1 then REMOVE,bd2[1:*],head  ; remove extra EXTVER
    MWRFITS,im,files[i],head,/create
    print,strtrim(i+1,2),' ',files[i],' fixed'
  endif else print,strtrim(i+1,2),' ',files[i],' nothing to fix'
endfor

;stop

end
