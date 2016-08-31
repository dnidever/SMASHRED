pro fix_optfw,night,verbose=verbose

; Check that all of the FWHM values in the daophot opt
; files are okay and consistent
; Use the summary files

setdisp

dir = '/data/smash/cp/red/photred/'+night+'/'
cd,dir

readcol,'fields',shortname,fieldname,format='A,A',/silent
nfields = n_elements(fieldname)

; Loop through the fields
for i=0,nfields-1 do begin
  print,' ',strtrim(i+1,2),' ',shortname[i],' ',fieldname[i]
  cd,dir+shortname[i]

  ; Get the number of exposures
  exp1 = file_search(shortname[i]+'-*_01.fits',count=nexp)
  if nexp eq 0 then goto,BOMB
  dum = strsplitter(file_basename(exp1,'_01.fits'),'-',/extract)
  uexpnum = reform(dum[1,*])
  nexpnum = n_elements(uexpnum)

  ; Exposure loop, load the opt files
  for j=0,nexpnum-1 do begin
    optfiles = file_search(shortname[i]+'-'+uexpnum[j]+'_??.opt',count=noptfiles)
    print,'  ',strtrim(j+1,2),' ',uexpnum[j],' ',strtrim(noptfiles,2),' opt files'
    ; Chip/optfile loop
    chstr = replicate({field:shortname[i],expnum:uexpnum[j],chip:0L,fwhm:0.0,bad:0B},noptfiles)
    dum = strsplitter(file_basename(optfiles,'.opt'),'_',/extract)
    chstr.chip = long(reform(dum[1,*]))
    for k=0,noptfiles-1 do begin
      READCOL,optfiles[k],key,equal,val,format='A,A,F',/silent
      fwind = where(key eq 'FW',nfwind)
      fwhm = val[fwind[0]]
      chstr[k].fwhm = fwhm
    endfor  ; chip loop
    med_fwhm = median(chstr.fwhm)
    sig_fwhm = mad(chstr.fwhm)
    bd_fwhm = where(abs(chstr.fwhm-med_fwhm)/sig_fwhm gt 5 or abs(chstr.fwhm-med_fwhm) gt 4,nbd_fwhm)
    ; Check if there are any problems
    plot,chstr.chip,chstr.fwhm,ps=1,xr=[0,63],xs=1,tit=uexpnum[j]
    oplot,[0,63],[0,0]+med_fwhm,linestyle=2,co=80
    if nbd_fwhm gt 0 then begin
      print,' bad chips. ',strjoin(chstr[bd_fwhm].chip,' ')
      oplot,[chstr[bd_fwhm].chip],[chstr[bd_fwhm].fwhm],ps=4,co=250,sym=2,thick=3
      bdchstr = chstr[bd_fwhm]
      writecol,-1,bdchstr.chip,bdchstr.fwhm,bdchstr.fwhm-med_fwhm,abs(bdchstr.fwhm-med_fwhm)/sig_fwhm,fmt='(I4,3F10.3)'
      chstr[bd_fwhm].bad = 1
    endif
    stop

  endfor  ; expnum loop

  stop


  fstr = mrdfits(sumfiles[j],1,/silent)
  chstr = mrdfits(sumfiles[j],2,/silent)
  nexp = n_elements(fstr)
  add_tag,fstr,'fwhm_sigma',0.0,fstr
  add_tag,fstr,'fwhm_nbad',0L,fstr
  add_tag,chstr,'delta_fwhm',0.0,chstr
  add_tag,chstr,'delta_fwhm_nsig',0.0,chstr
  add_tag,chstr,'fwhm_bad',0B,chstr
  ; Loop through exposures
  for k=0,nexp-1 do begin
    MATCH,chstr.expnum,fstr[k].expnum,ind1,ind2,/sort
    med_fwhm = median(chstr[ind1].fwhm)
    sig_fwhm = mad(chstr[ind1].fwhm)
    bd_fwhm = where(abs(chstr[ind1].fwhm-med_fwhm)/sig_fwhm gt 5 or abs(chstr[ind1].fwhm-med_fwhm) gt 4,nbd_fwhm)
    chstr[ind1].delta_fwhm = chstr[ind1].fwhm-med_fwhm
    chstr[ind1].delta_fwhm_nsig = abs(chstr[ind1].fwhm-med_fwhm)/sig_fwhm
    if nbd_fwhm gt 0 then print,'  ',strtrim(k+1,2),' ',fstr[k].expnum,'  ',strtrim(nbd_fwhm,2),' bad chips. ',strjoin(chstr[ind1[bd_fwhm]].base,' ') else $
         print,'  ',strtrim(k+1,2),' ',fstr[k].expnum,'  ',strtrim(nbd_fwhm,2),' bad chips.'
    fstr[k].fwhm_sigma = sig_fwhm
    fstr[k].fwhm_nbad = nbd_fwhm
      
    plot,chstr[ind1].chip,chstr[ind1].fwhm,ps=1,xr=[0,63],xs=1,tit=fstr[k].expnum
    oplot,[0,63],[0,0]+med_fwhm,linestyle=2,co=80
    if nbd_fwhm gt 0 then begin
      oplot,[chstr[ind1[bd_fwhm]].chip],[chstr[ind1[bd_fwhm]].fwhm],ps=4,co=250,sym=2,thick=3
      bdchstr = chstr[ind1[bd_fwhm]]
      writecol,-1,bdchstr.chip,bdchstr.fwhm,bdchstr.fwhm-med_fwhm,abs(bdchstr.fwhm-med_fwhm)/sig_fwhm,fmt='(I4,3F10.3)'

      chstr[ind1[bd_fwhm]].fwhm_bad = 1
;stop,'bad fwhm'
    endif
  endfor  ; exp loop
  push,allfstr,fstr
  push,allchstr,chstr
  BOMB:
endfor ; field loop

stop

end
