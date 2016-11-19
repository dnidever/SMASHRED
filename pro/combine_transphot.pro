pro combine_transphot,files,outfile

; This combines the output transformation equations for all the filters

files = file_search('transphot_?_eqns.fits',count=nfiles)
;nfiles = n_elements(files)
undefine,fitstr,chipstr,ntstr
allfilters = ''
for i=0,nfiles-1 do begin
  print,'Loading ',files[i]
  fitstr1 = mrdfits(files[i],1)
  push,fitstr,fitstr1
  chipstr1 = mrdfits(files[i],2)
  push,chipstr,chipstr1
  ntstr1 = mrdfits(files[i],3)
  push,ntstr,ntstr1
  allfilters += strtrim(fitstr1[0].filter,2)
endfor
head = headfits(files[0],ext=0)
sxaddpar,head,'filter',allfilters
if n_elements(outfile) eq 0 then outfile = 'smashred_transphot_eqns.fits'
FITS_WRITE,outfile,0,head
MWRFITS,fitstr,outfile,/silent
MWRFITS,chipstr,outfile,/silent
MWRFITS,ntstr,outfile,/silent

stop

end
