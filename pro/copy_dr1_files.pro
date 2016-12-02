pro copy_dr1_files

; Make the list of files to copy for the DR1 flat files

indir = '/data/smash/cp/red/photred/'
outdir = '/data/smash/dr1/'
;indir = '/datalab/users/dnidever/smash/cp/red/photred/'
;outdir = '/datalab/users/dnidever/smash/dr1/'

; Restore the DR1 fields
dr1 = mrdfits('smash_dr1_fields.fits',1)
ndr1 = n_elements(dr1)

; Get smash reduction info
smashred_getredinfo,allredinfo,/silent
bd = where(strtrim(allredinfo.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,allredinfo

; Loop through the fields and construct the list of files
for i=0,ndr1-1 do begin
  chstr = mrdfits(indir+'catalogs/final/v4/'+dr1[i].field+'_combined_chips.fits.gz',1)
  fstr = mrdfits(indir+'catalogs/final/v4/'+dr1[i].field+'_combined_exposures.fits.gz',1)

  allreddir = chstr.night+'/'+chstr.field
  ui = uniq(allreddir,sort(allreddir))
  reddir = allreddir[ui]
  nreddir = n_elements(reddir)

  match,allredinfo.field,dr1[i].field,ind1,ind2,/sort,count=nredinfo
  redinfo = allredinfo[ind1]

  ; Loop through each reduction directory
  for j=0,nreddir-1 do begin
    chind = where(allreddir eq reddir[j],nchind)
    ; Make the directory
    if file_test(outdir+reddir[j]) eq 0 then file_mkdir,outdir+reddir[j]

    ; Loop through the chips
    for k=0,nchind-1 do begin
      ; Regular files
      cfiles1 = indir+reddir[j]+'/'+chstr[chind[k]].base+$
                 ['.fits','.gaiawcs.head','s.fits.fz','.psf','.als','.opt','.als.opt','.log']
      file_copy,cfiles1,outdir+reddir[j]
      ; Allframe files
      if file_test(indir+reddir[j]+'/'+chstr[chind[k]].base+'.alf') eq 1 then begin
        cfiles2 = indir+reddir[j]+'/'+chstr[chind[k]].base+'.alf'
        file_copy,cfiles2,outdir+reddir[j]
      endif
      ; Reference exposures
      if chstr[chind[k]].expnum eq chstr[chind[k]].refexpnum then begin
        cfiles3 = indir+reddir[j]+'/'+chstr[chind[k]].base+$
                   ['.mch','.tfr','.raw']
      ; DO I NEED TO CHECK FOR TFR.ORIG FILES ???!!!
        file_copy,cfiles3,outdir+reddir[j]
      
        ; Allframe
        if file_test(indir+reddir[j]+'/'+chstr[chind[k]].base+'_comb.fits') eq 1 then begin
          cfiles4 = indir+reddir[j]+'/'+chstr[chind[k]].base+'_comb'+$
                    ['.fits','.bmp.fits','.psf','_allf.sex','_allf.als','']
          file_copy,cfiles4,outdir+reddir[j]
          cfiles5 = indir+reddir[j]+'/'+chstr[chind[k]].base+$
                    ['.zero','.weights','.scale','.shift']
          file_copy,cfiles5,outdir+reddir[j]
        endif
      endif

      stop
    endfor ; chips
    stop

    ; FIELD.fits and FIELD_summary.fits files
    ; "fields" file

  endfor  ; reduction directory

  stop

endfor  ; field

; There are fields per 

stop

end
