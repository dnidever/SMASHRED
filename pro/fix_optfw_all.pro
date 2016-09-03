pro fix_optfw_all,pl=pl,verbose=verbose

; Check that all of the FWHM values in the daophot opt
; files are okay and consistent
; Use the summary files

setdisp,/silent

dirs = file_search('/data/smash/cp/red/photred/20??????',/test_directory,count=ndirs)
print,strtrim(ndirs,2),' nights'

undefine,allchstr
for i=0,ndirs-1 do begin
  night = file_basename(dirs[i])
  undefine,chstr
  FIX_OPTFW,night,chstr
  push,allchstr,chstr
endfor  ; directory loop

;save,allchstr,file='fix_optfw_all.dat'

stop

end
