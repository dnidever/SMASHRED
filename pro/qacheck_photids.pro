pro qacheck_photids,night,verbose=verbose,allfields=allfields

; Make sure that none of the IDs coming phot photcalib.pro are *****  

if n_elements(night) gt 0 then begin
  dirs = file_search('/data/smash/cp/red/photred/'+night,/test_directory,count=ndirs)
endif else begin
  dirs = file_search('/data/smash/cp/red/photred/20??????',/test_directory,count=ndirs)
endelse
print,strtrim(ndirs,2),' nights'

for i=0,ndirs-1 do begin
  print,strtrim(i+1,2),'  '+dirs[i]
  cd,dirs[i]
  readcol,'fields',shortname,fieldname,format='A,A',/silent
  nfields = n_elements(fieldname)
  ; Loop through the fields
  for j=0,nfields-1 do begin
    field1 = shortname[j]
    if keyword_set(verbose) then print,' ',strtrim(j+1,2),'  '+field1
    cd,dirs[i]+'/'+field1
    photfiles = file_search(field1+'-*_??.phot',count=nphotfiles)
    if nphotfiles eq 0 then goto,bomb
    bad = 0
    for k=0,nphotfiles-1 do begin
      spawn,['grep','"\*"',photfiles[k]],out,errout,/noshell
      nout = n_elements(out)
      if nout le 1 then begin
        cmt = 'okay'
      endif else begin
        cmt = 'bad *****'
        bad = 1
      endelse
      if keyword_set(verbose) then print,'  ',strtrim(k+1,2),'  ',file_basename(photfiles[k]),' ',cmt
    endfor
    if bad eq 0 then cmt='okay' else cmt='bad *****'
    print,strtrim(j+1,2),'  ',dirs[i]+'/'+field1,'  ',cmt

    bomb:

    ;stop
  endfor ; field loop
  print,''
  ;stop
endfor  ; directory loop


stop

end
