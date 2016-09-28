pro run_wcsfit_gaia_single,file,redo=redo,rmslim=rmslim,loop=loop

dir = file_dirname(file)+'/'
cd,current=curdir
cd,dir
if keyword_set(loop) then nloop = 15 else nloop=1

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

; Get the GAIA reference file
restore,gaiafile

; Loop
str = replicate({rms:0.0,header:ptr_new()},nloop)
for i=0,nloop-1 do begin
  if keyword_set(loop) then print,'LOOP ITERATION = ',strtrim(i+1,2)
  undefine,error,head
  wcsfit,fitsfile,refname='GAIA/GAIA',inprefcat=refcat,head=head,error=error,/noupdate,rmslim=rmslim

  ; Get the RMS
  ;  HISTORY WCSFIT: RMS=0.027 arcsec on  Mon Sep 26 03:09:50 2016
  rmsind = first_el(where(stregex(head,'WCSFIT: RMS',/boolean) eq 1,nrmsind),/last)  ; last one
  wcsline = head[rmsind[0]]
  lo = strpos(wcsline,'RMS=')
  tmp = strmid(wcsline,lo+4)
  rms = float( first_el(strsplit(tmp,' ',/extract)) )

  ; Save the info
  str[i].rms = rms
  str[i].header = ptr_new(head)
endfor
; Get the best RMS
bestind = first_el(minloc(str.rms))
head = *str[bestind].header
print,'' & print,'BEST RMS = ',strtrim(str[bestind].rms,2),' arcsec'

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
