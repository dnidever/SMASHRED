pro run_imfwhm,night

dir = '/data/smash/cp/red/photred/'

; Journal
journal,dir+night+'/run_imfwhm.log'


print,'Running run_imfwhm on ',night

files = file_search(dir+night+'/F*/F*-*_??.fits',count=nfiles)
print,strtrim(nfiles,2),' FITS files found'

str = replicate({night:'',file:'',fwhm:-1.0,ellip:-1.0},nfiles)
str.night = night
str.file = files

for i=0,nfiles-1 do begin
  undefine,fwhm,ellip,error
  IMFWHM,files[i],fwhm,ellip,error=error
  if n_elements(error) eq 0 then begin
    str[i].fwhm = fwhm
    str[i].ellip = ellip
  endif
endfor

mwrfits,str,dir+night+'/run_imfwhm.fits',/create

journal  ; end journal file

;stop

end
