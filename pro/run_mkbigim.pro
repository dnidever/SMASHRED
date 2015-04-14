pro run_mkbigim

; Run mkbigim to make binned versions of the stacked images
; for all the osi/osj image files

files1 = file_search('*osi*.fits',count=nfiles1)
;files2 = file_search('*osj*.fits',count=nfiles2)
if nfiles1 gt 0 then push,files,files1
;if nfiles2 gt 0 then push,files,files2
bd = where(stregex(files,'_bin.fits',/boolean) eq 1,nbd)
if nbd gt 0 then remove,bd,files
nfiles = n_elements(files)

print,strtrim(nfiles,2),' files'

for i=0,nfiles-1 do begin
  print,strtrim(i+1,2),'/',strtrim(nfiles,2),' ',files[i]
  mkbigim,files[i]
  print,''
endfor

stop

end
