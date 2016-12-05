pro copy_dr1_files

; Make the list of files to copy for the DR1 flat files

catdir = '/data/smash/cp/red/photred/'
;indir = '/data/smash/cp/red/photred/'
;outdir = '/data/smash/dr1/'
indir = '/datalab/users/dnidever/smash/cp/red/photred/'
outdir = '/datalab/users/dnidever/smash/dr1/'

; Restore the DR1 fields
dr1 = mrdfits('smash_dr1_fields.fits',1)
dr1.field = strtrim(dr1.field,2)
ndr1 = n_elements(dr1)

; Get smash reduction info
smashred_getredinfo,allredinfo,/silent
bd = where(strtrim(allredinfo.field,2) eq '',nbd)
if nbd gt 0 then remove,bd,allredinfo

; Loop through the fields and construct the list of files
cinfiles = strarr(2e6)
coutfiles = strarr(2e6)
cnt = 0LL
push,reddirarr
for i=0,ndr1-1 do begin
  print,strtrim(i+1,2),'/'+strtrim(ndr1,2),' ',dr1[i].field
  chstr = mrdfits(catdir+'catalogs/final/v4/'+dr1[i].field+'_combined_chips.fits.gz',1)
  chstr.field = strtrim(chstr.field,2)
  chstr.base = strtrim(chstr.base,2)
  chstr.night = strtrim(chstr.night,2)
  chstr.expnum = strtrim(chstr.expnum,2)
  chstr.refexpnum = strtrim(chstr.refexpnum,2)
  fstr = mrdfits(catdir+'catalogs/final/v4/'+dr1[i].field+'_combined_exposures.fits.gz',1)

  ; Get the reduction directories for this field
  allreddir = chstr.night+'/'+chstr.field
  ui = uniq(allreddir,sort(allreddir))
  reddir = allreddir[ui]
  nreddir = n_elements(reddir)
  push,reddirarr,reddir

  match,allredinfo.field,dr1[i].field,ind1,ind2,/sort,count=nredinfo
  redinfo = allredinfo[ind1]

  ; Loop through each reduction directory, e.g. 20130317/F1/
  for j=0,nreddir-1 do begin
    chind = where(allreddir eq reddir[j],nchind)
    nightdir = chstr[chind[0]].night+'/'
    if chstr[chind[0]].alf_nsources gt 0 then alf=1 else alf=0
    print,' ',strtrim(j+1,2),'/',strtrim(nreddir,2),' ',reddir[j],' ',alf

    ; Loop through the chips in CHSTR
    for k=0,nchind-1 do begin
      undefine,cinfiles1,coutfiles1

      ; Regular files
      cfiles1 = reddir[j]+'/'+chstr[chind[k]].base+$
                 ['.fits','.gaiawcs.head','s.fits.fz','.psf','.als','.opt','.als.opt','.log']
      push,cinfiles1,indir+cfiles1
      push,coutfiles1,outdir+cfiles1

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
                    ['.fits','.bpm.fits','.psf','_allf.sex','_allf.als']
          push,cinfiles1,indir+cfiles5
          push,coutfiles1,outdir+cfiles5
  
          ; combination files
          cinfiles6 = reddir[j]+'/'+chstr[chind[k]].base+$
                      ['.zero','.weights','.scale','.shift']
          coutfiles6 = reddir[j]+'/'+chstr[chind[k]].base+'_comb'+$
                      ['.zero','.weights','.scale','.shift']
          push,cinfiles1,indir+cinfiles6
          push,coutfiles1,outdir+coutfiles6
        endif
      endif

      ; Add to the final arrays
      nfiles1 = n_elements(cinfiles1)
      cinfiles[cnt:cnt+nfiles1-1] = cinfiles1
      coutfiles[cnt:cnt+nfiles1-1] = coutfiles1
      cnt += nfiles1

    endfor ; chips

    undefine,cinfiles1,coutfiles1  ; use for the "extra" 4 files

    ; The FIELD.fits.gz and FIELDS_summary.fits files
    ;  get info from the "fields" file
    readcol,catdir+nightdir+'fields',shname,lngname,format='A,A',/silent
    gfield = where(strtrim(shname,2) eq chstr[chind[0]].field,ngfield)
    if ngfield eq 0 then stop,'CANNOT FIND FIELD '+chstr[chind[0]].field+' IN FIELDS FILE'
    cfieldfiles = nightdir+lngname[gfield[0]]+['.fits.gz','_summary.fits']
    push,cinfiles1,indir+cfieldfiles
    push,coutfiles1,outdir+cfieldfiles

    ; The "fields" and photred.setup files
    push,cinfiles1,indir+nightdir+['fields','photred.setup']
    push,coutfiles1,outdir+nightdir+['fields','photred.setup']

    nfiles1 = n_elements(cinfiles1)
    cinfiles[cnt:cnt+nfiles1-1] = cinfiles1
    coutfiles[cnt:cnt+nfiles1-1] = coutfiles1
    cnt += nfiles1
  endfor  ; reduction directory

  ;stop

endfor  ; field

print,'Final CNT = ',cnt

; Trim the extra lines
cinfiles = cinfiles[0:cnt-1]
coutfiles = coutfiles[0:cnt-1]

; Make sure they are unique, will get duplicate fields and
; photred.setup files
ui = uniq(cinfiles,sort(cinfiles))
cinfilesorig = cinfiles
coutfilesorig = coutfiles
cinfiles = cinfiles[ui]
coutfiles = coutfiles[ui]

; Checking that all the files exist
t = catdir+strmid(cinfiles,strlen(indir))
bd = where(file_test(t) eq 0,nbd)
print,strtrim(nbd,2),' files do NOT exist'

; Construct the copy commands
cmdcp = 'cp '+cinfiles+' '+coutfiles

; Make commands to make the directories
cmddir = 'mkdir '+outdir
nightdirs = outdir+file_dirname(reddirarr)
ui = uniq(nightdirs,sort(nightdirs))
nightdirs = nightdirs[ui]
push,cmddir,'mkdir '+nightdirs
push,cmddir,'mkdir '+outdir+reddirarr

cmd = [cmddir,cmdcp]
writeline,'copy_dr1_files.txt',cmd

stop

end
