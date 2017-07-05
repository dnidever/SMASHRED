pro create_deep_multiplenight_symlinks

; Create symlinks in deep/ for Fields
; that were observed on multiple nights
; and were reprocessed.  This is to save
; space.

rootdir = smashred_rootdir()+'cp/red/photred/'

;check = mrdfits(rootdir+'catalogs/pro/check_calibrated_v5.fits',1)
check = mrdfits(rootdir+'catalogs/final/v5/check_calibrated_v5.fits',1)
check.field = strtrim(check.field,2)
ncheck = n_elements(check)
smashred_getredinfo,info,/silent
for i=0,ncheck-1 do begin
  gd = where(info.field eq check[i].field and info.sh eq 0,ngd)
  if ngd gt 1 then push,fields,check[i].field
endfor

;; Remove the short LMC fields
;fnum = long(strmid(check.field,5))
;;shlmc = where(check.nuchips eq 0 and check.ngchips lt 150 and check.nrchips lt 150 and check.nichips lt 150 and check.nzchips lt 150,nshlmc)
;gshlmc = where(fnum ge 184 and fnum lt 246 and check.calflag gt 0,ngshlmc)
;shlmcfields = check[gshlmc].field
;MATCH,shlmcfields,fields,ind1,ind2,/sort,count=nmatch
;if nmatch gt 0 then remove,ind2,fields

;fieldnum = [100,101,104,109,113,116,76,77,84,92,98]
;fields = 'Field'+strtrim(fieldnum,2)
nfields = n_elements(fields)

;smashred_getredinfo,info,/silent

for i=0,nfields-1 do begin
  ifield = fields[i]
  print,strtrim(i+1,2),' ',ifield
  ind = where(info.field eq ifield and info.sh eq 0,nind)
  if nind eq 1 then begin
    stop,'only one night'
  endif
  nnight = nind
  print,strtrim(i+1,2),'/',strtrim(nnight,2),' ',ifield,' ',strtrim(nnight,2),' '

  ;deep/Field110/F1/F1-00519030_01
  ;20160218/F2/F2-00519030_01
  ;The header of the FITS file is slightly modified because the AIRMASS line was updated
  ;(just the date in the comment changed), which changed the CHECKSUM line and the
  ;DATASUM comment line.
  ;All other files are identical except for the '.alf' and '.log' files.

  ; Loop over nights
  for j=0,nnight-1 do begin
    info1 = info[ind[j]]
    night1 = info1.night
    fstr1 = *info1.fstr
    chstr1 = mrdfits(info1.file,2,/silent)
    nexp = n_elements(fstr1)
    nch = n_elements(chstr1)


    ; Loop over all chip files (for all exposures)
    for k=0,nch-1 do begin
      ; Files to create symlinks for
      ;  33 files, everything except for .alf and .log
      names = ['.als','.als.inp','.als.opt','.ap','.cmn.ap','.cmn.coo','.cmn.log','.cmn.lst',$
               '.coo','.fits','.gaiawcs.head','.grp','.lst','.lst1','.lst1.chi','.lst2',$
               '.lst2.chi','.nei','.nst','.opt','.plst','.plst.chi','.psf','.psf.log',$
               '.psfini.ap','_cat.dat','_refcat.dat','_refcat_gaia.dat','a.als','a.als.inp',$
               'a.ap','a.log','s.fits.fz']

      ; Get the original and deep names
      orignames = rootdir+night1+'/'+chstr1[k].field+'/'+chstr1[k].base+names
      deepbase = repstr(chstr1[k].base,chstr1[k].field,'F1')
      deepnames = rootdir+'deep/'+ifield+'/F1/'+deepbase+names

      ; Check that the file sizes are the same
      originfo = file_info(orignames)
      deepinfo = file_info(deepnames)
      szmismatch = where(originfo.size ne deepinfo.size,nszmismatch)
      if nszmismatch gt 0 then stop,'size mismatch'

      ; 1) Temporarily rename the file (e.g., ".name")
      tempnames = file_dirname(deepnames)+'/.'+file_basename(deepnames)
      FILE_MOVE,deepnames,tempnames
      temptest = file_test(tempnames)
      if min(temptest) eq 0 then stop,'problem with temp files'
      ; 2) Create symlink
      FILE_LINK,orignames,deepnames
      ; 3) Check that the symlink worked and the file size, etc is correct
      linkinfo = file_infor(deepnames)
      szmismatch2 = where(originfo.size ne linkinfo.size or deepinfo.size ne linkinfo.size,nszmismatch2)
      if nszmismatch2 gt 0 then stop,'link size mismatch'
      ; 4) Delete original renamed temporary file
      FILE_DELETE,tempnames

      stop

    endfor ; chip files loop
  endfor  ; night loop

  stop

endfor  ; field loop

stop

end
