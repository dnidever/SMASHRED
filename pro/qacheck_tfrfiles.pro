pro qacheck_tfrfiles

; Check that we have the right TFR files and that allframe
; wasn't accidentally run when it shouldn't have.

smashred_getredinfo,info,/silent
bd = where(strtrim(info.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,info
ninfo = n_elements(info)
print,strtrim(ninfo,2),' photred files to check'

; Loop through the files
for i=0,ninfo-1 do begin
  print,strtrim(i+1,2),' ',info[i].file
  chstr = mrdfits(info[i].file,2,/silent)
  if total(chstr.alf_nsources) gt 0 then alf=1 else alf=0
  base = file_basename(info[i].file)
  photredfield = first_el(strsplit(base,'_',/extract))
  ;len = strlen(photredfield)
  ;if strmid(base,len-2,2) eq 'sh' then short=1 else short=0

  if info[i].sh eq 1 and alf eq 1 then stop,'short field but ALF=1'

  ; ALS, are there tfr.orig files in this directory
  if alf eq 0 then begin
    nightdir = file_dirname(info[i].file)
    tfrorigfiles = file_search(nightdir+'/'+chstr[0].field+'/'+chstr[0].field+'-????????_??.tfr.orig',count=ntfrorigfiles)
    if ntfrorigfiles gt 0 then stop,'TFR.ORIG files found'
    ;stop
  endif

endfor

stop

end
