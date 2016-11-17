pro combine_transphot_201610

; This combines the output transformation equations for all the filters
; for 201610 and with the older ones

files = file_search('/data/smash/cp/red/photred/stdred_201610/transphot_?_eqns_201610.fits',count=nfiles)
;nfiles = n_elements(files)
undefine,newfitstr,newchipstr,newntstr
allfilters = ''
for i=0,nfiles-1 do begin
  print,'Loading ',files[i]
  fitstr1 = mrdfits(files[i],1)
  push,newfitstr,fitstr1
  ;chipstr1 = mrdfits(files[i],2)
  ;push,newchipstr,chipstr1
  ntstr1 = mrdfits(files[i],3)
  push,newntstr,ntstr1
  allfilters += strtrim(fitstr1[0].filter,2)
endfor

; Load the old "final" structures
outfile = '/data/smash/cp/red/photred/stdred/smashred_transphot_eqns.fits'
fitstr0 = mrdfits(outfile,1)
chipstr0 = mrdfits(outfile,2)
ntstr0 = mrdfits(outfile,3)
head0 = headfits(outfile,exten=0)

; add FITSTR and NTSTR, but leave CHIPSTR unchanged
fitstr = [fitstr0,newfitstr]
chipstr = chipstr0
ntstr = [ntstr0,newntstr]

FITS_WRITE,outfile,0,head0
MWRFITS,fitstr,outfile,/silent
MWRFITS,chipstr,outfile,/silent
MWRFITS,ntstr,outfile,/silent

stop

end
