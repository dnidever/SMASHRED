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

; Load list of bad exposures
badexp = importascii('~/projects/SMASHRED/obslog/smash_badexposures.txt',/header)

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
  if (ifield eq 'Field169') or (ifield eq 'Field214') or (ifield eq 'Field246') then goto,fieldbomb
  ind = where(info.field eq ifield and info.sh eq 0,nind)
  print,strtrim(i+1,2),' ',ifield,' ',strtrim(nind,2)
  if nind eq 1 then begin
    stop,'only one night'
  endif
  nnight = nind

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
    print,strtrim(j+1,2),'/',strtrim(nnight,2),' ',ifield,' ',night1

    ; Loop over all chip files (for all exposures)
    for k=0,nch-1 do begin
      print,'  ',strtrim(k+1,2),'/',+strtrim(nch,2),' ',chstr1[k].base
      deepbase = repstr(chstr1[k].base,chstr1[k].field,'F1')
      ; Files to create symlinks for
      ;  33 files, everything except for .alf and .log
      names = ['.als','.als.inp','.als.opt','.ap','.cmn.ap','.cmn.coo','.cmn.log','.cmn.lst',$
               '.coo','.grp','.lst','.lst1','.lst1.chi','.lst2',$
               '.lst2.chi','.nei','.nst','.opt','.plst','.plst.chi','.psf','.psf.log',$
               '.psfini.ap','_cat.dat','_refcat.dat','a.als','a.als.inp',$
               'a.ap','a.log','s.fits.fz']
      ; not all files have gaia wcs info
      if file_test(rootdir+'deep/'+ifield+'/F1/'+deepbase+'.gaiawcs.head') eq 1 then names=[names,'.gaiawcs.head']
      if file_test(rootdir+'deep/'+ifield+'/F1/'+deepbase+'_refcat_gaia.dat') eq 1 then names=[names,'_refcat_gaia.dat']

      ; Get the original and deep names
      orignames = rootdir+night1+'/'+chstr1[k].field+'/'+chstr1[k].base+names
      deepnames = rootdir+'deep/'+ifield+'/F1/'+deepbase+names

      ; Bad exposures
      match,badexp.expnum,chstr1[k].expnum,badind1,badind2,count=nbad
      if total(file_test(deepnames)) eq 0 and nbad eq 1 then begin
        print,'  bad exposure. skipping'
        goto,BOMB
      endif

      ; Chip=02 files are sometimes missing
      if total(file_test(deepnames)) eq 0 and chstr1[k].chip eq 2 then begin
        print,'  chip 02 files are missing.  This is not unexpected.'
        goto,BOMB
      endif

      ; Check the files
      originfo = file_info(orignames)
      deepinfo = file_info(deepnames)
      ; Check if all the deep files exist
      ;   for now don't worry about missing ones
      deepmissing = where(deepinfo.exists eq 0,ndeepmissing)
      len = strlen(rootdir+'deep/'+ifield+'/F1/')
      if ndeepmissing gt 0 then print,'  ',strtrim(ndeepmissing,2),' deep file(s) missing: ',strjoin(strmid(deepnames[deepmissing],len),' ')
      ; Check if some files are already symlinks
      deeplinks = where(deepinfo.exists eq 1 and deepinfo.symlink eq 1,ndeeplinks)
      if ndeeplinks gt 0 then print,'  ',strtrim(ndeeplinks,2),' symlinks already exist'
      ; Get files that exist and are NOT symlinks
      gddeep = where(deepinfo.exists eq 1 and deepinfo.symlink eq 0,ngddeep,comp=bddeep,ncomp=nbddeep)
      if ngddeep eq 0 then begin
        print,'  No files to create symlinks for.  Skipping to FITS'
        goto,dofits
      endif
      ; Check that the file sizes are the same
      szmismatch = where(originfo[gddeep].size ne deepinfo[gddeep].size,nszmismatch)
      if nszmismatch gt 0 then stop,'size mismatch'

      ; Arrays for the deep files that exist
      names1 = names[gddeep]
      deepnames1 = deepnames[gddeep]
      orignames1 = orignames[gddeep]
      originfo1 = originfo[gddeep]
      deepinfo1 = deepinfo[gddeep]

      ; 1) Temporarily rename the file (e.g., ".name")
      tempnames1 = file_dirname(deepnames1)+'/.'+file_basename(deepnames1)
      FILE_MOVE,deepnames1,tempnames1
      temptest = file_test(tempnames1)
      if min(temptest) eq 0 then stop,'problem with temp files'
      ; 2) Create symlink
      origlink1 = '../../../'+night1+'/'+chstr1[k].field+'/'+chstr1[k].base+names1
      FILE_LINK,origlink1,deepnames1
      ; 3) Check that the symlink worked and the file size, etc is correct
      linkinfo1 = file_info(deepnames1)
      szmismatch2 = where(originfo1.size ne linkinfo1.size or deepinfo1.size ne linkinfo1.size,nszmismatch2)
      if nszmismatch2 gt 0 then stop,'link size mismatch'
      ; 4) Delete original renamed temporary file
      FILE_DELETE,tempnames1

      ; Deal with FITS file separately
      dofits:
      origfitsname = rootdir+night1+'/'+chstr1[k].field+'/'+chstr1[k].base+'.fits.fz'
      deepfitsname = rootdir+'deep/'+ifield+'/F1/'+deepbase+'.fits'
      if file_test(deepfitsname+'.fz') eq 1 then begin
        print,'  fits.fz file exists already'
        goto,bomb
      endif
      if file_test(origfitsname) eq 0 or file_test(deepfitsname) eq 0 then stop,'problem with FITS files'

      ; 1) Create symlink
      origfitslink = '../../../'+night1+'/'+chstr1[k].field+'/'+chstr1[k].base+'.fits.fz'
      deepfitslink = rootdir+'deep/'+ifield+'/F1/'+deepbase+'.fits.fz'
      FILE_LINK,origfitslink,deepfitslink
      ; 2) Check that it links to a real file
      linkinfo = FILE_INFO(deepfitslink)
      if linkinfo.exists eq 0 or linkinfo.dangling_symlink eq 1 then stop,'problem with symlink'
      ; 3) Delete original fits file
      FILE_DELETE,deepfitsname

      BOMB:

      ; Remove any comba.fits
      combafile = rootdir+'deep/'+ifield+'/F1/'+deepbase+'_comba.fits'
      if file_test(combafile) eq 1 then FILE_DELETE,combafile

      ;stop

    endfor ; chip files loop
  endfor  ; night loop

  ;stop
  FIELDBOMB:

endfor  ; field loop

stop

end
