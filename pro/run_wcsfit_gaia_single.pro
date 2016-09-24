pro run_wcsfit_gaia_single,file,redo=redo

dir = file_dirname(file)+'/'
cd,current=curdir
cd,dir

if file_test(file) eq 0 then begin
  print,file,' NOT FOUND'
  return
endif

base = file_basename(file,'.fits')
fitsfile = base+'.fits'
gaiafile = base+'_refcat_gaia.dat'
outfile = base+'.gaiawcs.head'

if file_test(gaiafile) eq 0 then begin
  print,gaiafile,' NOT FOUND'
  return
endif

if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  print,outfile,' EXISTS and /redo not set'
  return
endif

restore,gaiafile
undefine,error
wcsfit,fitsfile,refname='GAIA/GAIA',inprefcat=refcat,head=head,error=error,/noupdate

; There was an error
if n_elements(error) gt 0 then begin
  print,'An error occurred in wcsfit.pro'
  return
endif

; Save header
print,'Writing FITS header to ',outfile
writeline,outfile,head

cd,curdir

;stop

end
