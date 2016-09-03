pro fix_optfw,night,allchstr,pl=pl,verbose=verbose

; Check that all of the FWHM values in the daophot opt
; files are okay and consistent
; Use the summary files

setdisp,/silent

undefine,allchstr
if n_elements(pl) eq 0 then pl=0

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

; Loop through the fields
undefine,allchstr
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
    QACHECK_OPTFWEXP,shortname[i]+'-'+uexpnum[j],chstr1,pl=pl,verbose=verbose,offstr=offstr
    push,allchstr,chstr1

    ; We have some bad ones to rerun
    bdind1 = where(chstr1.bad eq 1,nbdind1)
redo = 0
    if nbdind1 gt 0 and keyword_set(redo) then begin

      ; Rerun photred_mkopt on BAD chips
      for k=0,nbdind1-1 do begin
        file = file_basename(chstr1[bdind1[k]].file,'.opt')+'.fits'
        print,'   ',strtrim(k+1,2),' FIXING FWHM for ',file
        ;file = chstr1[bdind1[k]].base+'_'+string(chstr1[bdind1[k]].chip,format='(i02)')+'.fits'
        ; Check that VA and FITRADIUS_FWHM are okay
        va = chstr1[bdind1[k]].va
        if va lt 0 then va=2
        fitradius_fwhm = chstr1[bdind1[k]].fitradius_fwhm
        if fitradius_fwhm lt 0 then fitradius_fwhm = 1.0
        PHOTRED_MKOPT,file,va=va,fitradius_fwhm=fitradius_fwhm,/verbose
        chstr1[bdind1[k]].redo = 1   ; will need to redo this one
      endfor

      ; Check FWHM again
      print,'Checking again'
      QACHECK_OPTFWEXP,shortname[i]+'-'+uexpnum[j],chstr2,pl=pl,verbose=verbose,offstr=offstr

      ; If any bad ones are left, force the fit
      bdind2 = where(chstr2.bad eq 1,nbdind2)
      for k=0,nbdind2-1 do begin
        ; Check that VA and FITRADIUS_FWHM are okay
        va = chstr2[bdind2[k]].va
        if va lt 0 then va=2
        fitradius_fwhm = chstr2[bdind2[k]].fitradius_fwhm
        if fitradius_fwhm lt 0 then fitradius_fwhm = 1.0
        PHOTRED_MKOPT,file,va=va,fitradius_fwhm=fitradius_fwhm,inp_fwhm=chstr2[bdind2[k]].fitfwhm,/verbose
        chstr2[bdind2[k]].redo = 1   ; will need to redo this one
      endfor

      stop

    endif ; some bad

  endfor  ; expnum loop

  ; If there are any to redo then run DAOPHOT on them.

  ; Figure out which files need to have MATCH rerun on them.

  ;stop

  BOMB:
endfor ; field loop

cd,curdir

;stop

end
