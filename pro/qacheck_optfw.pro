pro qacheck_optfw,night,verbose=verbose,allfields=allfields

; Check that all of the FWHM values in the daophot opt
; files are okay and consistent
; Use the summary files

setdisp

if n_elements(night) gt 0 then begin
  dirs = file_search('/data/smash/cp/red/photred/'+night,/test_directory,count=ndirs)
endif else begin
  dirs = file_search('/data/smash/cp/red/photred/20??????',/test_directory,count=ndirs)
endelse
print,strtrim(ndirs,2),' nights'

for i=0,ndirs-1 do begin
  print,strtrim(i+1,2),'  '+dirs[i]
  cd,dirs[i]

  ; Get the summary files
  sumfiles = file_search(dirs[i]+'/*_summary.fits',count=nsumfiles)

  ; Loop through the fields
  for j=0,nsumfiles-1 do begin
    print,' ',strtrim(j+1,2),' ',file_basename(sumfiles[j],'_summary.fits')
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
  endfor ; field loop
  print,''
  ;stop
endfor  ; directory loop

plotc,lindgen(n_elements(allfstr)),allfstr.fwhm_sigma,allfstr.fwhm_nbad,ps=1
plotc,lindgen(n_elements(allchstr)),allchstr.fwhm,allchstr.fwhm_bad,ps=1
bd = where(allchstr.fwhm_bad eq 1,nbd)
print,strtrim(nbd,2),' bad FWHM values'
; 902, ~215 short, ~690 long

rnd = randomu(seed,n_elements(allchstr))*0.4-0.2
plotc,allchstr.chip+rnd,allchstr.delta_fwhm,allchstr.fwhm_bad,ps=1
; there are some low and high ones that are getting missed here,
; probably from high sigma exposures.

;save,allfstr,allchstr,file='qacheck_optfw.dat'

stop

end
