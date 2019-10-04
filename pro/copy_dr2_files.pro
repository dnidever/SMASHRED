pro copy_dr2_files

; Make the list of files to copy for the DR1 flat files

catdir = '/dl1/users/dnidever/smash/cp/red/photred/'
indir = '/dl1/users/dnidever/smash/cp/red/photred/'
outdir = '/dl1/users/dnidever/smash/dr2/photred/'
version = 'v6'

; Restore the DR2 fields
dr2 = importascii('~/projects/SMASHRED/data/smash_fields_final.txt',/header)
add_tag,dr2,'field','',dr2
dr2.field = 'Field'+strtrim(dr2.num,2)

; Calibration information
calstr = mrdfits(catdir+'check_calibrated_v6.fits',1)
calstr.field = strtrim(calstr.field,2)
MATCH,dr2.field,calstr.field,ind1,ind2,/sort,count=nmatch
dr2 = dr2[ind1]  ; 197 fields
ndr2 = n_elements(dr2)

; Get smash reduction info
smashred_getredinfo,allredinfo,/silent
bd = where(strtrim(allredinfo.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,allredinfo

; Loop through the fields and construct the list of files
cinfiles = strarr(4000000L)
coutfiles = strarr(4000000L)
testinfiles = intarr(4000000L)
cnt = 0LL
undefine,reddirarr
for i=0,ndr2-1 do begin
;for i=143,ndr2-1 do begin
  print,strtrim(i+1,2),'/'+strtrim(ndr2,2),' ',dr2[i].field
  chstr = mrdfits(catdir+'catalogs/final/'+version+'/'+dr2[i].field+'_combined_chips.fits.gz',1)
  chstr.field = strtrim(chstr.field,2)
  chstr.base = strtrim(chstr.base,2)
  chstr.night = strtrim(chstr.night,2)
  chstr.expnum = strtrim(chstr.expnum,2)
  chstr.refexpnum = strtrim(chstr.refexpnum,2)
  if tag_exist(chstr,'PHOTDIR') then chstr.photdir = strtrim(chstr.photdir,2)
  fstr = mrdfits(catdir+'catalogs/final/'+version+'/'+dr2[i].field+'_combined_exposures.fits.gz',1)


;; JUST USE WHATEVER IS IN THE COMBINED CHIPS FILES!!!!

  if tag_exist(chstr,'PHOTDIR') then begin
    allreddir = strtrim(chstr.photdir,2)+'/'+strtrim(chstr.field,2)
    for j=0,n_elements(allreddir)-1 do begin  ; make relative to smash/cp/red/photred/
      lo = strpos(allreddir[j],'smash/cp/red/photred')
      allreddir[j] = strmid(allreddir[j],lo+21)
    endfor        
  endif else begin
    allreddir = chstr.night+'/'+chstr.field
  endelse

  ;; Unique REDDIR array
  ui = uniq(allreddir,sort(allreddir))
  reddir = allreddir[ui]
  nreddir = n_elements(reddir)
  push,reddirarr,reddir

;stop


;  ;; Use deep/ directory for long exposures
;  ;; ONE night:
;  ;;  no PHOTDIR in chstr
;  ;; Multiple nights:
;  ;;  PHOTDIR column and gives deep/ directory
;
;  ;; USE deep/ directory, same code works for both cases
;  gd = where(allredinfo.field eq dr2[i].field and allredinfo.sh eq 1,ngd)
;  dpsumfile = '/dl1/users/dnidever/smash/cp/red/photred/deep/'+dr2[i].field+'/'+dr2[i].field+'_summary.fits'
;  ;sumfiles = [dpsumfile,allredinfo[gd].file]
;  ;SMASHRED_GETREDINFO,redinfo,sumfiles=sumfiles,/silent
;
;
;  ;; ONE night, no PHOTRED in chstr 
;  singledeepnight = 0
;  if tag_exist(chstr,'PHOTDIR') eq 0 then begin
;    singledeepnight = 1
;  endif else begin
;    ndeepchstr = total(stregex(chstr.photdir,'deep',/boolean),/integer)
;    if ndeepchstr eq 0 then singledeepnight = 1
;  endelse
;  if dr2[i].field eq 'Field169' or dr2[i].field eq 'Field214' then goto,skipdeep
;
;  ;; ONE night, no PHOTRED in chstr 
;  if keyword_set(singledeepnight) then begin
;    allreddir = chstr.night+'/'+chstr.field
;    ;; Change to deep/ for long exposures
;    fstr1 = mrdfits(dpsumfile,1,/silent)
;    fstr1.expnum = strtrim(fstr1.expnum,2)
;    nfstr1 = n_elements(fstr1)
;    ;chstr1 = mrdfits(dpsumfile,2,/silent)  ;; use this to get the PHOTRED Field name
;    ;deepdir = 'deep/'+dr2[i].field+'/'+strtrim(chstr1[0].field,2)
;    for j=0,nfstr1-1 do begin
;      MATCH,chstr.expnum,fstr1[j].expnum,ind1,ind2,/sort,count=nmatch
;      allreddir[ind1] = 'deep/'+dr2[i].field+'/'+chstr[ind1].field
;    endfor
;
;  ;; Multiple nights: PHOTDIR column gives directory name
;  endif else begin
;    allreddir = strtrim(chstr.photdir,2)+'/'+strtrim(chstr.field,2)
;    for j=0,n_elements(allreddir)-1 do begin  ; make relative to smash/cp/red/photred/
;      lo = strpos(allreddir[j],'smash/cp/red/photred')
;      allreddir[j] = strmid(allreddir[j],lo+21)
;    endfor    
;  endelse
;
;  ;; Field169 does NOT have a deep/ directory
;  skipdeep:
;  if dr2[i].field eq 'Field169' or dr2[i].field eq 'Field214' then allreddir = chstr.night+'/'+chstr.field
;
;  ;; Field118 has deep/ multi-night directory, but it wasn't used to create final catalog
;
;  if singledeepnight eq 1 then begin
;    info = file_info(catdir+'deep/'+dr2[i].field+'F1')
;    if info.exists eq 1 and info.symlink eq 0 then print,'---',dr2[i].field,' DEEP PROBLEM --------------------------'
;  endif
;
;  ;; Unique REDDIR array
;  ui = uniq(allreddir,sort(allreddir))
;  reddir = allreddir[ui]
;  nreddir = n_elements(reddir)
;  push,reddirarr,reddir
;
;  ;; check if there are multiple directories in reddir that have deep
;  if total(stregex(reddir,'deep/',/boolean)) gt 1 then print,'-----',dr2[i].field,'  MULTIPLE DEEP DIRECTORIES --------------'


  ; Loop through each reduction directory, e.g. 20130317/F1/  or  deep/Field1/F1/
  for j=0,nreddir-1 do begin
    chind = where(allreddir eq reddir[j],nchind)
    ;nightdir = chstr[chind[0]].night+'/'
    nightdir = file_dirname(reddir[j])+'/'
    if chstr[chind[0]].alf_nsources gt 0 then alf=1 else alf=0
    print,' ',strtrim(j+1,2),'/',strtrim(nreddir,2),' ',reddir[j],' ',alf

    ; Loop through the chips in CHSTR
    for k=0,nchind-1 do begin
      undefine,cinfiles1,coutfiles1

      ; Regular files
      cfiles1 = reddir[j]+'/'+chstr[chind[k]].base+$
                 ['.fits.fz','s.fits.fz','.psf','.als','.opt','.als.opt','.log']
      ;; gaiawcs.head only exist for data taken before Sep 2016
      ;;   after that Gaia DR1 was the default astrometric catalog
      if long(chstr[chind[k]].night) lt 20160901L then push,cfiles1,reddir[j]+'/'+chstr[chind[k]].base+'.gaiawcs.head'
      push,cinfiles1,indir+cfiles1
      push,coutfiles1,outdir+cfiles1

      ;; The symlinks are missing from s.fits.fz for
      ;; deep/FieldX/Field1/F1-XXXXXXXX_XXs.fits.fz files
      ;; create them
      stest = file_test(catdir+reddir[j]+'/'+chstr[chind[k]].base+'s.fits.fz')
      if stest eq 0 then begin
        print,'no s.fits.fz file'
        linkname = file_readlink(catdir+reddir[j]+'/'+chstr[chind[k]].base+'.fits.fz')
        cd,current=curdir
        cd,catdir+reddir[j]
        file_link,repstr(linkname,'.fits.fz','s.fits.fz'),catdir+reddir[j]+'/'+chstr[chind[k]].base+'s.fits.fz'
        cd,curdir
      endif

      ; Allframe files
      if alf eq 1 then begin
        cfiles2 = reddir[j]+'/'+chstr[chind[k]].base+'.alf'
        ncfiles2 = n_elements(cfiles2)
        push,cinfiles1,indir+cfiles2
        push,coutfiles1,outdir+cfiles2
      endif
      ; Reference exposures
      if chstr[chind[k]].expnum eq chstr[chind[k]].refexpnum then begin
        cfiles3 = reddir[j]+'/'+chstr[chind[k]].base+$
                   ['.mch','.raw']
        push,cinfiles1,indir+cfiles3
        push,coutfiles1,outdir+cfiles3

        ; ALLSTAR TFR file
        if alf eq 1 then begin
          cinfiles4 = reddir[j]+'/'+chstr[chind[k]].base+'.tfr.orig'
          coutfiles4 = reddir[j]+'/'+chstr[chind[k]].base+'.tfr'
          push,cinfiles1,indir+cinfiles4
          push,coutfiles1,outdir+coutfiles4
        endif else begin
          cfiles4 = reddir[j]+'/'+chstr[chind[k]].base+'.tfr'
          push,cinfiles1,indir+cfiles4
          push,coutfiles1,outdir+cfiles4
        endelse

        ; Allframe
        if alf eq 1 then begin
          cfiles5 = reddir[j]+'/'+chstr[chind[k]].base+'_comb'+$
                    ['.fits','.bpm.fits.gz','.psf','_allf.sex','_allf.als']
          push,cinfiles1,indir+cfiles5
          push,coutfiles1,outdir+cfiles5
  
          ;; combination files
          if file_test(catdir+reddir[j]+'/'+chstr[chind[k]].base+'.zero') then begin
            cbase = reddir[j]+'/'+chstr[chind[k]].base
          endif else begin
            cbase = reddir[j]+'/'+chstr[chind[k]].base+'_comb'
          endelse
          cinfiles6 = cbase+['.zero','.weights','.scale']
          coutfiles6 = reddir[j]+'/'+chstr[chind[k]].base+'_comb'+$
                      ['.zero','.weights','.scale']
          if file_test(cbase+'.shift') eq 1 then begin
            ;; the recent alftiletype=WCS with allframe_combine.pro
            ;; does NOT create .shift files
             ;; this only affects: 20161029, 20161030, 20161031,
             ;; 20170308, 20170804
            push,cinfiles6,cbase+'.shift'
            push,coutfiles6,reddir[j]+'/'+chstr[chind[k]].base+'_comb.shift'
          endif
          push,cinfiles1,indir+cinfiles6
          push,coutfiles1,outdir+coutfiles6
        endif
      endif

      ;; Check that the files exist
      test = file_test(cinfiles1)
      b = where(test eq 0,nb)
      if nb gt 0 then begin
        print,'SOME FILES MISSING'
        printline,cinfiles1[b]
      endif

      ; Add to the final arrays
      nfiles1 = n_elements(cinfiles1)
      cinfiles[cnt:cnt+nfiles1-1] = cinfiles1
      coutfiles[cnt:cnt+nfiles1-1] = coutfiles1
      testinfiles[cnt:cnt+nfiles1-1] = test
      cnt += nfiles1

    endfor ; chips

    undefine,cinfiles1,coutfiles1  ; use for the "extra" 4 files

    ; The FIELD.fits.gz and FIELDS_summary.fits files
    ;  get info from the "fields" file
    if stregex(reddir[j],'deep',/boolean) eq 0 then begin 
      readcol,catdir+nightdir+'fields',shname,lngname,format='A,A',/silent
      gfield = where(strtrim(shname,2) eq chstr[chind[0]].field,ngfield)
      if ngfield eq 0 then stop,'CANNOT FIND FIELD '+chstr[chind[0]].field+' IN FIELDS FILE'
    endif else begin
      lngname = dr2[i].field
      gfield = 0
    endelse
    cfieldfiles = nightdir+lngname[gfield[0]]+['.fits.gz','_summary.fits']
    push,cinfiles1,indir+cfieldfiles
    push,coutfiles1,outdir+cfieldfiles

    ;; Create "fields" file for deep/
    if stregex(reddir[j],'deep',/boolean) eq 1 and file_test(catdir+nightdir+'fields') eq 0 then begin
      chstr1 = mrdfits(dpsumfile,2,/silent)  ;; use this to get the PHOTRED Field name  
      writeline,catdir+nightdir+'fields',strtrim(chstr1[0].field,2)+'   '+dr2[i].field
    endif

    ; The "fields" and photred.setup files
    push,cinfiles1,indir+nightdir+['fields','photred.setup']
    push,coutfiles1,outdir+nightdir+['fields','photred.setup']

    nfiles1 = n_elements(cinfiles1)
    cinfiles[cnt:cnt+nfiles1-1] = cinfiles1
    coutfiles[cnt:cnt+nfiles1-1] = coutfiles1
    testinfiles[cnt:cnt+nfiles1-1] = file_test(cinfiles1)
    cnt += nfiles1
  endfor  ; reduction directory

  ;stop

endfor  ; field

print,'Final CNT = ',cnt

; Trim the extra lines
cinfiles = cinfiles[0:cnt-1]
coutfiles = coutfiles[0:cnt-1]
testinfiles = testinfiles[0:cnt-1]

; Make sure they are unique, will get duplicate fields and
; photred.setup files
ui = uniq(cinfiles,sort(cinfiles))
cinfilesorig = cinfiles
coutfilesorig = coutfiles
cinfiles = cinfiles[ui]
coutfiles = coutfiles[ui]
testinfiles = testinfiles[ui]

; Checking that all the files exist
;print,'Checking that all ',strtrim(n_elements(cinfiles),2),' exist'
;t = catdir+strmid(cinfiles,strlen(indir))
;bd = where(file_test(t) eq 0,nbd)
;print,strtrim(nbd,2),' files do NOT exist'

; Construct the copy commands
;cmdcp = 'cp '+cinfiles+' '+coutfiles
cmdcp = 'ln -s '+cinfiles+' '+coutfiles   ; USE SYMLINKS!!

; Make commands to make the directories
cmddir = 'mkdir '+outdir
push,cmddir,'mkdir '+outdir+'deep'
nightdirs = outdir+file_dirname(reddirarr)
ui = uniq(nightdirs,sort(nightdirs))
nightdirs = nightdirs[ui]
push,cmddir,'mkdir '+nightdirs
push,cmddir,'mkdir '+outdir+reddirarr

cmd = [cmddir,cmdcp]
writeline,'copy_dr2_files.txt',cmd

stop

end
