;+
;
; MAKE_CUSTOM_REFCAT
;
; INPUTS:
;  files      The array of FITS files to use for the custom reference
;               catalog.  FILE_cat.dat files need to exist.
;  outfile    The output IDL save file.
;  /clobber   If the output file already exists, overwrite it.
;
; OUTPUTS:
;  The reference catalog (in REFCATALL IDL structure) is written
;  to the filename specified in "outfile".
;
; USAGE:
;  IDL>make_custom_refcat,files,outfile,clobber=clobber
;
; By D.Nidever  March 2016
;-

pro make_custom_refcat,files,outfile,clobber=clobber

; Make custom reference catalog from exposures with good WCS

nfiles = n_elements(files)
noutfile = n_elements(outfile)
; Not enough inputs
if nfiles eq 0 or noutfile eq 0 then begin
  print,'Syntax - make_custom_refcat,files,outfile,clobber=clobber'
  return
endif

; Outfile already exists
if file_test(outfile) eq 1 and not keyword_set(clobber) then begin
  print,outfile,' ALREADY EXISTS and /clobber not set'
  return
endif

; Loop through the files
for i=0,nfiles-1 do begin
  base = file_basename(files[i],'.fits')
  dir = file_dirname(files[i])
  catfile = dir+'/'+base+'_cat.dat'
  if file_test(catfile) gt 0 then begin
    restore,catfile
    print,base,' ',strtrim(n_elements(cat),2)
    head = headfits(files[i])
    head_xyad,head,cat.x,cat.y,ra,dec,/deg

    ; Only keep "good" stars
    gd = where(cat.sharp ge 0.2 and cat.sharp le 1.0 and $
               cat.round ge -1.0 and cat.round le 1.0 and $
               cat.mag lt 50.0 and cat.err lt 1.0,ngd)
    print,' ',strtrim(ngd,2),' good sources'

    refcat1 = replicate({raj2000:0.0d0,dej2000:0.0d0,rmag:0.0,rerr:0.0},ngd)
    refcat1.raj2000 = ra[gd]
    refcat1.dej2000 = dec[gd]
    refcat1.rmag = cat[gd].mag
    refcat1.rerr = cat[gd].err

    ; Catalog already exists, look for matches
    if n_elements(refcatall) gt 0 then begin
      ; Check for matches
      srcmatch,refcatall.raj2000,refcatall.dej2000,refcat1.raj2000,refcat1.dej2000,0.5,ind1,ind2,/sph,count=nmatch

      ; Matches, remove duplicates
      if nmatch gt 0 then begin
        print,strtrim(nmatch,2),' matches'
        ; Remove duplicates
        if nmatch lt n_elements(refcat1) then remove,ind2,refcat1 else undefine,refcat1
        ; If some left, add them
        if n_elements(refcat1) gt 0 then push,refcatall,refcat1

      ; No matches, add all
      endif else push,refcatall,refcat1

    ; First catalog, use all
    endif else refcatall=refcat1

  endif else print,catfile,' NOT FOUND'

endfor

print,'Saving catalog to ',outfile
save,refcatall,file=outfile

;stop

end
