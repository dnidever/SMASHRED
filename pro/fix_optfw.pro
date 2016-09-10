pro fix_optfw,night,allchstr,allexpstr,pl=pl,verbose=verbose,redo=redo

; Check that all of the FWHM values in the daophot opt
; files are okay and consistent
; Use the summary files

;setdisp,/silent

if n_elements(pl) eq 0 then pl=0
if n_elements(redo) eq 0 then redo=0  ; don't redo by default

cd,current=curdir
dir = '/data/smash/cp/red/photred/'+night+'/'
if file_test(dir,/directory) eq 0 then begin
  print,dir,' NOT FOUND'
  return
endif
cd,dir

print,'Running FIX_OPTFW on ',night

readcol,'fields',shortname,fieldname,format='A,A',/silent
nfields = n_elements(fieldname)
print,strtrim(nfields,2),' fields'

; Loop through the fields and fix FWHM values
undefine,allchstr,allexpstr
for i=0,nfields-1 do begin
  cd,dir+shortname[i]

  ; Get the number of exposures
  exp1 = file_search(shortname[i]+'-*_01.fits',count=nexp)
  if nexp eq 0 then goto,BOMB
  dum = strsplitter(file_basename(exp1,'_01.fits'),'-',/extract)
  uexpnum = reform(dum[1,*])
  nexpnum = n_elements(uexpnum)

  print,' ',strtrim(i+1,2),' ',shortname[i],' ',fieldname[i],'   ',strtrim(nexpnum,2),' exposures'

  ; Exposure loop, load the opt files
  for j=0,nexpnum-1 do begin

    ; Check the FWHM values
    QACHECK_OPTFWEXP,shortname[i]+'-'+uexpnum[j],chstr,pl=pl,verbose=verbose,offstr=offstr
    if n_elements(chstr) eq 0 then goto,BOMB1
    chstr.dir = dir+shortname[i]
    chstr.field = shortname[i]

    ; Create exposure level structure
    expstr = {dir:'',night:'',field:'',expnum:'',filter:'',med_fwhm:0.0,sig_fwhm:0.0,nbad:0L}
    expstr.dir = dir+shortname[i]
    expstr.night = shortname[i]
    expstr.field = chstr[0].field
    expstr.expnum = chstr[0].expnum
    expstr.filter = chstr[0].filter
    expstr.med_fwhm = median(chstr.fwhm)
    expstr.sig_fwhm = mad(chstr.fwhm)
    expstr.nbad = total(chstr.bad)

    ; We have some bad ones to rerun
    bdind = where(chstr.bad eq 1,nbdind)
    if nbdind gt 0 and keyword_set(redo) then begin

      ; Rerun photred_mkopt on BAD chips and force
      ;   the fitted FWHM value.
      ;   NOTE: I reran some by hand and found that they gave
      ;   decent values but slightly lower (~15%).  Probably
      ;   because IMFWHM.PRO was modified and this is giving
      ;   systematically lower values.
      for k=0,nbdind-1 do begin
        file = file_basename(chstr[bdind[k]].file,'.opt')+'.fits'
        print,'   ',strtrim(k+1,2),' FIXING FWHM for ',file
        ; Check that VA and FITRADIUS_FWHM are okay
        va = chstr[bdind[k]].va
        if va lt 0 then va=2
        fitradius_fwhm = chstr[bdind[k]].fitradius_fwhm
        if fitradius_fwhm lt 0 then fitradius_fwhm = 1.0
        PHOTRED_MKOPT,file,va=va,fitradius_fwhm=fitradius_fwhm,inp_fwhm=chstr[bdind[k]].fitfwhm,/verbose
        chstr[bdind[k]].redo = 1   ; will need to redo this one
      endfor
    endif ; some bad
    PUSH,allchstr,chstr
    PUSH,allexpstr,expstr
    BOMB1:
  endfor  ; expnum loop
  BOMB:
endfor ; field loop

; If there are any bad ones to redo then run DAOPHOT on them.
bdind = where(allchstr.bad eq 1 and allchstr.redo eq 1,nbdind)
if nbdind gt 0 and keyword_set(redo) then begin
  print,'Rerunning DAOPHOT on ',strtrim(nbdind,2),' files'
  cd,dir
  fitsfiles = allchstr[bdind].dir+'/'+file_basename(allchstr[bdind].file,'.opt')+'.fits'
  alsfiles = allchstr[bdind].dir+'/'+file_basename(allchstr[bdind].file,'.opt')+'.als'

  ; Step 1.  Delete the existing ALS files, so PHOTRED_DAOPHOT
  ;           will naturally rerun them.
  print,'Step 1.  Deleting existing ALS files'
  FILE_DELETE,alsfiles,/allow,/quiet

  ; Step 2.  Add files to logs/DAOPHOT.inlist
  print,'Step 2. Adding files to logs/DAOPHOT.inlist'
  WRITELINE,'logs/DAOPHOT.inlist',fitsfiles,/append

  ; Step 3.  Run PHOTRED_DAOPHOT
  print,'Step 3.  Running PHOTRED_DAOPHOT'
  PHOTRED_DAOPHOT,redo=0,nmulti=30

  ; Should we check for failures?

  ; -- Run MATCH --

  ; Figure out all of the ALS files to input to MATCH
  undefine,allalsfiles
  for i=0,nbdind-1 do begin
    alsfiles1 = file_search(allchstr[bdind[i]].dir+'/'+allchstr[bdind[i]].field+'-*_'+$
                            string(allchstr[bdind[i]].chip,format='(i02)')+'.als',count=nalsfiles1)
    push,allalsfiles,alsfiles1
  endfor
  ui = uniq(allalsfiles,sort(allalsfiles))
  allalsfiles = allalsfiles[ui]

  ; Write to logs/MATCH.inlist
  WRITELINE,'logs/MATCH.inlist',allalsfiles,/append

  ; Run PHOTRED_MATCH
  PHOTRED_MATCH,/redo

  ; Rerun PHOTRED_ALLFRAME
  ;  it will know which fields to skip
  PHOTRED_ALLFRAME,/redo,nmulti=20

  ; Do I need to rerun APCOR?  Probably!

  ; Should we rerun ASTROM-SAVE as well or by hand??

  ;stop

endif

; Should double-check that the same reference frame was reused

cd,curdir

;stop

end
