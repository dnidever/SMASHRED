pro solve_transphot_night,arr,mfitstr,fitcolam=fitcolam,fixcolr=fixcolr,fixam=fixam,fixcolam=fixcolam,$
                        fixchipzp=fixchipzp,resid=resid,expstr=expstr,verbose=verbose,silent=silent,pl=pl,$
                        save=save,bootstrap=bootstrap,stp=stp,rejected=rejected,noreject=noreject

; Fit color and zpterm for each night+chip and airmass term
; per night. itertate until convergence
;
; INPUTS
;  arr         The structure of data to determine the transformation eqns. with.
;  =fixcolr    An input structure that contains color terms per chip
;                     to use instead of fitting them.
;  =fixchipzp  An input structure that contains zero-point terms per
;                     chip to use instead of fitting them.  A nightly
;                     zero-point term will still be determined.
;  =fixam      The airmass/extinction term to use instead of fitting one.
;  =fixcolam   The color*airmass term to use instead of fitting one.
;  /fitcolam   Fit the color*airmass term.  The default it so leave at
;                 0.0 and not fit this.
; /bootstrap   Use bootstrap analysis to determine coefficient uncertainties.
; /noreject    Don't do any rejection
;
; OUTPUTS:
;  mfitstr     The structure of fitted values, one element for each
;                chip for this night.
;  =resid      The residuals of the fit.
;  =expstr     Structure giving statistics for each exposure.
;  =rejected   The structure saying if a star was rejected or not.

undefine,mfitstr
undefine,expstr

narr = n_elements(arr)
imjd = arr[0].mjd

; Get unique chip tags
ui = uniq(arr.chip,sort(arr.chip))
uchips = arr[ui].chip
nchips = n_elements(uchips)

tags = tag_names(arr)
filter = strlowcase(tags[6])  ; this should be the filter name

; Initialize the fitting structure
mfitstr = replicate({mjd:-1L,chip:-1L,nightnum:0,nstars:0L,nrejected:0L,seeing:99.9,zpterm:999.0,zptermerr:999.0,$
                     colterm:999.0,coltermerr:999.0,amterm:0.0,amtermerr:999.0,colamterm:0.0,colamtermerr:999.0,$
                     rms:999.0,sig:999.0,chisq:999.0,medresid:999.0,nbrt:0L,brtrms:999.0,brtsig:999.0,brtchisq:999.0,status:-1,$
                     relzpterm:99.99,chiprelzpterm:99.99,chiprelcolterm:99.99,$
                     badsoln:-1,photometric:-1},nchips)


; Initialize the EXPOSURE structure
ui = uniq(arr.expnum,sort(arr.expnum))
uexpnum = arr[ui].expnum
nexpnum = n_elements(uexpnum)
expstr = replicate({expnum:'',mjd:-1L,nightnum:-1L,exptime:-1.0,airmass:-1.0,nstars:-1L,medresid:999999.0,sigresid:999999.0,rejected:0},nexpnum)
for i=0,nexpnum-1 do begin
  MATCH,arr.expnum,uexpnum[i],ind1,ind2,count=nmatch,/sort
  expstr[i].expnum = uexpnum[i]
  expstr[i].mjd = imjd
  expstr[i].exptime = arr[ind1[0]].exptime
  expstr[i].airmass = median([arr[ind1].airmass])
  expstr[i].nstars = nmatch
endfor


; Use input airmass term
if n_elements(fixam) gt 0 then begin
  mfitstr.amterm = fixam
  mfitstr.amtermerr = 0.0
endif
; Use input color*airmass term
if n_elements(fixcolam) gt 0 then begin
  mfitstr.colamterm = fixcolam
  mfitstr.colamterm = 0.0
endif

; Iterate until converge
endflag = 0
count = 0L
lastmfitstr = mfitstr
rejected = bytarr(narr)
ntzpterm = 0.0
WHILE (endflag eq 0) do begin

  ; Fit color term and zpterm for each chip separately
  ;---------------------------------------------------

  resid = fltarr(narr)+999.0       ; FULL residuals
  ntresid = fltarr(narr)+999.0     ; use to derive nightly zpterm, amterm, colterm, colamterm, chip-level zpterm removed
  amresid = fltarr(narr)+999.0     ; use to derive AIRMASS term, col and colam removed
  colamresid = fltarr(narr)+999.0  ; use to derive COLAM term, col and airmass removed
  FOR i=0,nchips-1 do begin
    ichip = uchips[i]
    MATCH,arr.chip,uchips[i],ind1,ind2,count=nmatch,/sort

    mfitstr[i].mjd = imjd
    mfitstr[i].chip = ichip
    mfitstr[i].nstars = nmatch
    if nmatch gt 0 then mfitstr[i].seeing=median([arr[ind1].seeing])

    ; Use input color term
    if n_elements(fixcolr) gt 0 then begin
      MATCH,fixcolr.chip,uchips[i],fixcolrind,dum,count=nfixcolr,/sort
      colterm = fixcolr[fixcolrind[0]].colterm
      coltermerr = fixcolr[fixcolrind[0]].coltermerr
      mfitstr[i].colterm = colterm
      mfitstr[i].coltermerr = coltermerr
    endif

    ; Use input chip-level zpterm
    if n_elements(fixchipzp) gt 0 then begin
      MATCH,fixchipzp.chip,uchips[i],fixchipzpind,dum,count=nfixchipzp,/sort
      zpterm = fixchipzp[fixchipzpind[0]].zpterm
      zptermerr = fixchipzp[fixchipzpind[0]].zptermerr
      ; chip-level relative zpterm has temporal dependence
      if tag_exist(fixchipzp,'zpmjdterm') then begin
        zpterm += imjd * fixchipzp[fixchipzpind[0]].zpmjdterm
        zptermerr = sqrt(zptermerr^2 + ( imjd * fixchipzp[fixchipzpind[0]].zpmjdtermerr )^2 )
      endif
      mfitstr[i].zpterm = zpterm          ; these two will get appended with the nightly zpterm/err below
      mfitstr[i].zptermerr = zptermerr
    endif

    ; Get non-rejected points
    MATCH,rejected[ind1],0,ind1b,ind2b,count=nmatch2,/sort

    ; Enough stars to fit the terms
    if nmatch gt 2 and nmatch2 gt 2 then begin

      ; non-rejected stars
      gdind1 = ind1[ind1b]

      ; Fit the color term
      ; CMAG - 6
      ; COL - 8
      x = arr[gdind1].(8)
      mag = arr[gdind1].mag
      y = arr[gdind1].mag-arr[gdind1].(6)
      yorig = y   ; with extinction still in 
      y -= arr[gdind1].airmass*mfitstr[i].amterm                        ; subtract AIRMASS term
      y -= arr[gdind1].(8) * arr[gdind1].airmass * mfitstr[i].colamterm  ; subtract COLOR*AIRMASS term
      err = arr[gdind1].err
      am = arr[gdind1].airmass

      ; Color term
      ;-------------
      ; Fit color term
      if n_elements(fixcolr) eq 0 then begin

        ; Enough color range to fit color term
        if range(x) gt 0.01 then begin
          ;coef0 = poly_fit(x,y,1,measure_errors=err,sigma=sigma0,yerror=yerror0,status=status0,yfit=yfit0)
          coef0 = robust_poly_fitq(x,y,1)
          yfit0 = poly(x,coef0)
          yerror0 = sqrt(mean((y-yfit0)^2))
          ; reject outliers and refit
          diff0 = y-yfit0
          gd = where(abs(diff0) lt 3*yerror0 and err lt 0.05,ngd)
          coef1 = robust_poly_fitq(x[gd],y[gd],1)  ; second robust fit
          ;yfit = poly(x,coef1)
          ;ydiff = y-yfit
          ; use poly_fit to get uncertainties
;invwt = err*am^2
          ;coef1 = poly_fit(x[gd],y[gd],1,measure_errors=invwt[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1)
          coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
          yfit = poly(x,coef1)
          ydiff = y-yfit

          ; The weighted fit actually gives better results

          zpterm = coef[0]       ; this might get redetermined below
          zptermerr = coeferr[0]
          colterm = coef[1]
          coltermerr = coeferr[1]
        ; Not enough range in color to fit color term
        endif else begin
          print,'Not enough range in color to fit color term.  Just measuring zero-point.'
          gd = where(y lt 100,ngd)
          coef = median(y[gd])
          ; use poly_fit to get uncertainties
          coef1 = dln_poly_fit(x[gd],y[gd],0,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
          ydiff = y-coef[0]

          zpterm = coef[0]       ; this might get redetermined below
          zptermerr = coeferr[0]
          colterm = 0.0
          coltermerr = 999.0
        endelse
      endif

      ; Zero-point term
      ;----------------
      ; Fit zpterm
      if n_elements(fixcolr) gt 0 and n_elements(fixchipzp) eq 0 then begin
        y -= arr[gdind1].(8) * colterm  ; subtract input color term

        coef0 = median(y)
        yfit0 = poly(x,coef0)
        yerror0 = sqrt(mean((y-yfit0)^2))
        ; reject outliers and refit
        diff0 = y-yfit0
        gd = where(abs(diff0) lt 3*yerror0 and err lt 0.05,ngd)
        coef = median(y[gd])
        yfit = poly(x,coef)
        ydiff = y-yfit
        ; use poly_fit to get uncertainties
        coef1 = dln_poly_fit(x[gd],y[gd],0,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
        zpterm = coef[0]
        zptermerr = coeferr[0]
      endif
      ; Use input chip-level zero-point term
      if n_elements(fixchipzp) gt 0 then begin
        y -= arr[gdind1].(8) * colterm  ; subtract input color term
        zpterm = mfitstr[i].zpterm      ; this is the input chip-level zpterm from above
        zptermerr = mfitstr[i].zptermerr
        ydiff = y-zpterm
        ;yerror = sqrt(mean((ydiff-median(ydiff))^2))  ; subtract out nightly offset
        yerror = sqrt(mean((ydiff-ntzpterm)^2))  ; subtract out nightly offset
        status = 1
      endif

      chisq = mean((ydiff/err)^2)
      ;if n_elements(fixchipzp) gt 0 then chisq = mean(((ydiff-median(ydiff))/err)^2)  ; subtract nightly zpterm
      if n_elements(fixchipzp) gt 0 then chisq = mean(((ydiff-ntzpterm)/err)^2)  ; subtract nightly zpterm

      ; Residuals
      yfit = arr[ind1].(8)*colterm + zpterm
      ; residuals for computing nightly zpterm
      ;   subtract amterm, colamterm, colterm and chip-level zpterm
      ntresid[ind1] = (arr[ind1].mag-arr[ind1].(6)) - arr[ind1].(8)*colterm - zpterm - arr[ind1].airmass*mfitstr[i].amterm $
                      - arr[ind1].(8)*arr[ind1].airmass*mfitstr[i].colamterm 
      ; residuals for computing airmass term
      ;   subtract colterm, zpterm
      ;     (nightly zpterm will be subtracted below if necessary)
      amresid[ind1] = (arr[ind1].mag-arr[ind1].(6)) - yfit
      ; residuals for computing col*airmass term
      ;   subtract colterm, zpterm (amterm done below)
      ;     (nightly zpterm will be subtracted below if necessary)
      colamresid[ind1] = (arr[ind1].mag-arr[ind1].(6)) - yfit

      ; Reject outliers
      sig = mad(ydiff)
      if count gt 0 and not keyword_set(noreject) then begin
        rejected1 = abs(ydiff-median(ydiff)) gt 5*sig
        rejected[gdind1] = rejected[gdind1] OR rejected1    ; combine the rejected bits
      endif

      ; Median residuals
      if n_elements(fixchipzp) gt 0 then medresid=median([ydiff-ntzpterm]) else medresid=median([ydiff])

      ; Stick information into chip structure
      mfitstr[i].zpterm = zpterm
      mfitstr[i].zptermerr = zptermerr
      mfitstr[i].colterm = colterm
      mfitstr[i].coltermerr = coltermerr
      mfitstr[i].rms = yerror
      mfitstr[i].sig = sig
      mfitstr[i].chisq = chisq
      mfitstr[i].medresid = medresid
      mfitstr[i].status = status
      mfitstr[i].nrejected = total(rejected[ind1])

      ; RMS of bright stars
      gbright = where(err lt 0.05 and mag lt min(mag)+3,nbright)
      mfitstr[i].nbrt = nbright
      if nbright gt 1 then begin
        ;brtrms = sqrt(mean((ydiff[gbright]-median(ydiff))^2))  ; subtract out nightly offset
        ;brtchisq = mean(((ydiff[gbright]-median(ydiff))/err[gbright])^2)
        if n_elements(fixchipzp) gt 0 then brtrms = sqrt(mean((ydiff[gbright]-ntzpterm)^2)) else $  ; subtract out nightly offset
          brtrms = sqrt(mean(ydiff[gbright]^2))
        if n_elements(fixchipzp) gt 0 then brtchisq = mean(((ydiff[gbright]-ntzpterm)/err[gbright])^2) else $
          brtchisq = mean(((ydiff[gbright]-ntzpterm)/err[gbright])^2)
        mfitstr[i].brtrms = brtrms
        mfitstr[i].brtsig = mad(ydiff[gbright])
        mfitstr[i].brtchisq = brtchisq
      endif
      if keyword_set(verbose) then print,strtrim(i+1,2),' ',ichip,' ',nmatch,' ',zpterm,' ',colterm,' ',yerror,' ',brtrms
    endif

  ENDFOR  ; chip loop


  ; Fit nightly zero-point term
  ;-----------------------------
  if n_elements(fixchipzp) gt 0 then begin
    ; NOT REJECTING ANY POINTS FOR NOW

    x = arr.mag*0
    y = ntresid  ; amterm, colamterm, colterm and chip-level zpterm removed
    err = arr.err

    coef0 = median(y)
    yfit0 = poly(x,coef0)
    yerror0 = sqrt(mean((y-yfit0)^2))
    ; reject outliers and refit
    diff0 = y-yfit0
    gd = where(abs(diff0) lt 3*yerror0 and err lt 0.05,ngd)
    coef = median([y[gd]])
    yfit = poly(x,coef)
    ydiff = y-yfit
    ; use poly_fit to get uncertainties
    coef1 = dln_poly_fit(x[gd],y[gd],0,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
    ntzpterm = coef[0]
    ntzptermerr = coeferr[0]

    ; Add this to the chip-level zpterm values
    mfitstr.zpterm += ntzpterm
    mfitstr.zptermerr = sqrt(mfitstr.zptermerr^2 + ntzptermerr^2)  ; add in quadrature

    ; Subtracted nightly zero-point term from residuals
    amresid -= ntzpterm
    colamresid -= ntzpterm
  endif

  ; Fit airmass term for all chips
  ;-------------------------------
  if n_elements(fixam) eq 0 then begin
    ; get non-rejected points
    MATCH,rejected,0,ind1,/sort
    x = arr[ind1].airmass
    y = amresid[ind1]
    mag = arr[ind1].mag
    err = arr[ind1].err
    ;mfitstr[i].minam = min(x)
    ;mfitstr[i].maxam = max(x)

    ; Enough airmass range to compute airmass term
    if range(x) gt 0.01 then begin

      ;coef0 = poly_fit(x,y,1,measure_errors=err,sigma=sigma0,yerror=yerror0,status=status0,yfit=yfit0)
      gd0 = where(y lt 100,ngd0)
      coef0 = robust_poly_fitq(x[gd0],y[gd0],1)
      yfit0 = poly(x,coef0)
      yerror0 = sqrt(mean((y-yfit0)^2))
      ; reject outliers and refit
      ydiff0 = y-yfit0
      gd = where(abs(ydiff0) lt 3*yerror0 and err lt 0.05,ngd)
      if yerror0 lt 0.0001 then gd = where(abs(ydiff0) lt 3*0.01 and err lt 0.05,ngd)
      coef1 = robust_poly_fitq(x[gd],y[gd],1)
      ;yfit = poly(x,coef)
      ;ydiff = y-yfit
      ; use weighted polynomial fit
      coef = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)
      yfit = poly(x,coef)
      ydiff = y-yfit
      ; the weighted fit actually gives better results

      mfitstr.amterm = coef[1]
      mfitstr.amtermerr = coeferr[1]

      ; Remove airmass term from residuals for colamterm
      colamresid -= arr.airmass*mfitstr.amterm   ; remove AIRMASS term

   ; Not enough airmass range
   endif else begin
     print,'Not enough of an airmass range to compute airmass term'
   endelse

  ; Use input airmass term
  endif else begin
    colamresid -= arr.airmass * mfitstr.amterm
  endelse

  ; Fit color*airmass term for all data
  ;------------------------------------
  if keyword_set(fitcolam) and n_elements(fixcolam) eq 0 then begin
    x = arr.(8) * arr.airmass
    y = colamresid
    mag = arr.mag
    err = arr.err

    coef0 = robust_poly_fitq(x,y,1)
    yfit0 = poly(x,coef0)
    yerror0 = sqrt(mean((y-yfit0)^2))
    ; reject outliers and refit
    ydiff0 = y-yfit0
    gd = where(abs(ydiff0) lt 3*yerror0 and err lt 0.05,ngd)
    if yerror0 lt 0.0001 then gd = where(abs(ydiff0) lt 3*0.01 and err lt 0.05,ngd)
    coef = robust_poly_fitq(x[gd],y[gd],1)
    yfit = poly(x,coef)
    ydiff = y-yfit
    ; use poly_fit to get uncertainties
    coef1 = dln_poly_fit(x[gd],y[gd],1,measure_errors=err[gd],sigma=coeferr,yerror=yerror,status=status,yfit=yfit1,bootstrap=bootstrap)

    mfitstr.colamterm = coef[1]
    mfitstr.colamtermerr = coeferr[1]
  endif

  ;resid2[ind1] = diff  ; save the residuals

  nrms = sqrt(mean((ydiff/err)^2))


  ; Final/full residuals
  for i=0,nchips-1 do begin
    MATCH,arr.chip,uchips[i],ind1,ind2,count=nmatch,/sort
    resid[ind1] = arr[ind1].mag - arr[ind1].(6) - mfitstr[i].zpterm - $
                  arr[ind1].(8)*mfitstr[i].colterm - arr[ind1].airmass*mfitstr[i].amterm
  endfor

  ; Get residual statistics for each EXPOSURE
  for i=0,nexpnum-1 do begin
    MATCH,arr.expnum,uexpnum[i],ind1,ind2,count=nmatch,/sort
    expstr[i].nstars = nmatch
    expstr[i].medresid = median([resid[ind1]])
    expstr[i].sigresid = mad([resid[ind1]-expstr[i].medresid],/zero)
  endfor

  ; Remove any outlier exposures
  ;  can't do it on the first iteration because zpterm
  ;  is messed up due to the contribution of amterm/colamterm
  if count gt 0 and not keyword_set(noreject) then begin
    sigexp = mad([expstr.medresid])
    bdexp = where(abs(expstr.medresid) gt (3*sigexp>0.1),nbdexp)
    if nbdexp gt 0 then begin
      print,'Rejected ',strtrim(nbdexp,2),' exposure(s). ',strjoin(expstr[bdexp].expnum,' ')
      for i=0,nbdexp-1 do begin
        MATCH,arr.expnum,expstr[bdexp[i]].expnum,ind1,ind2,count=nmatch,/sort
        rejected[ind1] = 1
        expstr[bdexp[i]].rejected = 1
      endfor
    endif
  endif  ; 2nd+ iteration

  ; Have we converged
  ; difference in zpterm, colterm, airmass terms
  maxdiff_thresh = 0.0001
  maxiter = 10
  diff = [ abs(mfitstr.zpterm-lastmfitstr.zpterm), $
           abs(mfitstr.colterm-lastmfitstr.colterm), $
           abs(mfitstr[0].amterm-lastmfitstr[0].amterm), abs(mfitstr[0].colamterm-lastmfitstr[0].colamterm) ]
  ;  need at least one iteration to get the zpterm right because
  ;  of the contribution of the airmass and colamterm terms
  if (count gt 1) and (max(diff) lt maxdiff_thresh or count ge maxiter) then endflag=1
  if not keyword_set(silent) then print,count,median(mfitstr.zpterm),median(mfitstr.colterm),mfitstr[0].amterm,$
          mfitstr[0].colamterm,max(diff[0:nchips-1]),max(diff[nchips:2*nchips-1]),max(diff[2*nchips]),max(diff[2*nchips+1]),format='(I5,8F11.5)'

  ; Save last values
  lastmfitstr = mfitstr

  count++

  ;stop

ENDWHILE

;if n_elements(fixchipzp) gt 0 then stop,'fixchipzp'

; Make a plot of residuals vs. color and residuals vs. airmass
; color-coded by the errors

if keyword_set(pl) or keyword_set(save) then begin
  if keyword_set(save) then begin
    !p.font = 0
    loadct,39,/silent
    psfile = 'solve_transphot_'+filter+'_'+strtrim(imjd,2)+'_resid'
    ps_open,psfile,thick=4,/color,/encap
  endif
  yr = [-4,4]*mad(resid)
  xr1 = [-0.5,3]  ;minmax(arr.(8))
  sicol = sort(arr.(8))
  xr1 = [arr[sicol[0.05*narr]].(8)-0.1,arr[sicol[0.95*narr]].(8)+0.1]  ; 5th and 95th percentile 
  xr2 = [1 < min(arr.airmass)-0.05, max(arr.airmass)+0.05]
  plotc,arr.(8),resid,arr.err,ps=3,xr=xr1,yr=yr,xs=1,ys=1,xtit='Color',ytit='Residuals',$
        pos=[0.08,0.08,0.53,0.90],colpos=[0.08,0.97,0.98,0.99]
  oplot,[-10,10],[0,0],linestyle=2
  plotc,arr.airmass,resid,arr.err,ps=3,xr=xr2,yr=yr,xs=1,ys=1,xtit='Airmass',$
        ytit=' ',ytickformat='(A1)',/noerase,/nocolorbar,pos=[0.53,0.08,0.98,0.90]
  oplot,[-10,10],[0,0],linestyle=2
  xyouts,0.5,0.91,'Residuals '+filter+'-band '+strtrim(imjd,2)+' RMS='+stringize(median(mfitstr.brtrms),ndec=3)+'  (color-coded by ERR)',align=0.5,charsize=1.2,/norm
  if keyword_set(save) then begin
    ps_close,/silent
  endif
endif

if keyword_set(stp) then stop

end

;-----------------

pro solve_transphot_allnights,arr,mstr,fitstr,fitcolam=fitcolam,fixcolr=fixcolr,fixam=fixam,fixcolam=fixcolam,$
                        fixchipzp=fixchipzp,resid=resid,expstr=expstr,verbose=verbose,silent=silent,pl=pl,save=save,$
                        bootstrap=bootstrap,rejected=rejected,noreject=noreject


tags = tag_names(arr)
filter = strlowcase(tags[6])  ; this should be the filter name
narr = n_elements(arr)

; Get unique nights
ui = uniq(arr.mjd,sort(arr.mjd))
umjds = arr[ui].mjd
nmjds = n_elements(umjds)

; Get unique chips
ui = uniq(arr.chip,sort(arr.chip))
uchips = arr[ui].chip
nchips = n_elements(uchips)

; Loop through nights and fit the color and airmass terms
mstr = replicate({mjd:0L,nightnum:0,date:'',nchips:0L,nstars:0L,nrejected:0L,seeing:99.9,rms:99.0,brtrms:99.0,$
                  airmass0:0.0,airmass1:0.0,amrange:0.0,zpterm:99.0,zptermsig:99.0,colterm:99.0,coltermsig:99.0,$
                  amterm:99.0,amtermsig:99.0,colamterm:99.0,colamtermsig:99.0,badsoln:-1,photometric:-1},nmjds)
observatory,'ctio',obs
undefine,fitstr
if not keyword_set(silent) then $
  print,'  NUM    MJD     DATE     NSTARS   NREJECT  ZPTERM    COLTERM     AMTERM    AMRANGE   COLAMTERM     RMS       BRTRMS     MEDSIG'
resid = fltarr(narr)+999.9
rejected = bytarr(narr)
for i=0,nmjds-1 do begin

  MATCH,arr.mjd,umjds[i],ind1,ind2,count=nmatch,/sort
  arr1 = arr[ind1]
  mstr[i].seeing = median(arr1.seeing)

  ; Fixing airmass term
  if n_elements(fixam) gt 0 then begin
    MATCH,fixam.mjd,umjds[i],ind3,ind4,count=nmatch2,/sort
    if nmatch2 gt 0 then fixam1=fixam[ind3].amterm else undefine,fixam1
  endif

  SOLVE_TRANSPHOT_NIGHT,arr1,mfitstr,/silent,fitcolam=fitcolam,fixcolr=fixcolr,fixam=fixam1,$
                      fixchipzp=fixchipzp,fixcolam=fixcolam,resid=resid1,expstr=expstr1,$
                      verbose=verbose,pl=pl,save=save,bootstrap=bootstrap,rejected=rejected1,$
                      noreject=noreject
  mfitstr.nightnum = i+1

  rejected[ind1] = rejected1
  resid[ind1] = resid1
  expstr1.mjd = umjds[i]
  expstr1.nightnum = i+1
  push,expstr,expstr1

  ;; Calculate residuals
  ;for j=0,nchips-1 do begin
  ;  MATCH,arr1.chip,uchips[j],chind,ind4,count=nmatch2,/sort
  ;  MATCH,mfitstr.chip,uchips[j],mind,/sort
  ;  if nmatch gt 0 then begin
  ;    resid[ind1[chind]] = arr[ind1[chind]].mag - arr[ind1[chind]].(6) - mfitstr[mind].zpterm - $
  ;          arr[ind1[chind]].(8)*mfitstr[mind].colterm - arr[ind1[chind]].airmass*mfitstr[mind].amterm
  ;  endif
  ;endfor

  ; Measure total RMS and NRMS per night
  totresid = (mfitstr.nstars-mfitstr.nrejected)*mfitstr.rms^2
  rms = sqrt(mean(total(totresid)/total(mfitstr.nstars-mfitstr.nrejected)))  ; rms for all stars of this night
  totbrtresid = mfitstr.nbrt*mfitstr.brtrms^2
  brtrms = sqrt(mean(total(totbrtresid)/total(mfitstr.nbrt)))
  medsig = median([mfitstr.sig])

  ; Save information to MJD structure
  mstr[i].mjd = umjds[i]
  ui = uniq(mfitstr.chip,sort(mfitstr.chip))
  mstr[i].nchips = n_elements(ui)
  mstr[i].nightnum = i+1
  mstr[i].nstars = nmatch
  mstr[i].nrejected = total(mfitstr.nrejected)
  mstr[i].rms = rms
  mstr[i].brtrms = brtrms
  mstr[i].zpterm = median([mfitstr.zpterm])
  mstr[i].zptermsig = mad(mfitstr.zpterm)
  mstr[i].colterm = median([mfitstr.colterm])
  mstr[i].coltermsig = mad(mfitstr.colterm)
  mstr[i].amterm = mfitstr[0].amterm
  mstr[i].amtermsig = mfitstr[0].amtermerr
  mstr[i].colamterm = mfitstr[0].colamterm
  mstr[i].colamtermsig = mfitstr[0].colamtermerr
  mstr[i].airmass0 = min(arr1.airmass)
  mstr[i].airmass1 = max(arr1.airmass)
  mstr[i].amrange = range(arr1.airmass)

  ; Convert MJD to YYYYMMDD date
  jd = umjds[i]+2400000.5d0-0.5+obs.tz/24.  
  caldat,jd,month,day,year,hour,min,sec
  date = strtrim(year,2)+string(month,format='(I02)')+string(day,format='(I02)')
  mstr[i].date = date

  ; Print summary info to screen
  print,i+1,umjds[i],date,nmatch,mstr[i].nrejected,median([mfitstr.zpterm]),median([mfitstr.colterm]),$
        mfitstr[0].amterm,mstr[i].amrange,mfitstr[0].colamterm,rms,brtrms,medsig,format='(I5,I8,A10,I8,I8,8F11.4)'
  ; Save structure
  push,fitstr,mfitstr

  ;stop

endfor

;stop

end

;----------------------

pro solve_transphot,file,fitcolam=fitcolam,errlim=errlim,arr=arr,mstr=mstr_fixcolzpam,resid=resid_fixcolzpam,$
                         rejected=rejected_fixcolzpam,noreject=noreject,fitchipzpmjdterm=fitchipzpmjdterm,$
                         fixchipzpterm=fixchipzpterm,nosmith=nosmith,stp=stp

; Check the STDRED data to see if we need separate color terms
; for each chip

;ui = uniq(arr.mjd,sort(arr.mjd))
;umjds = arr[ui].mjd
;nmjd = n_elements(umjds)
;tags = tag_names(arr)
;filter = strlowcase(tags[6])  ; this should be the filter name
;goto,averageamterm

if n_elements(file) eq 0 then begin
  print,'Syntax - solve_transphot,file'
  return
endif

; Outline of the program
; 1.) Removing outlier exposures
; 2.) Run on all the nights fitting everything individually
;  -set photometric tag
;  -flag bad solutions with not enough stars
; 3.) Determine median zpterm and colterm per chip
; 4.) Fixing colterm. Redetermining zpterm and amterm
;  -set photometric tag
;  -flag bad solutions with not enough stars
; 5.) Redetermine median zpterm per chip
; 6.) Fixing chip-level relative zpterm and colterm.  Retermining
;       nightly zpterm and amterm
; 7.) Averaging airmass terms
; 8.) Redetermine nightly zpterm holding everything else fixed
; 9.) Making final chip and night-level structures and write output files

; chip-level color terms
; chip-level zpterm
; nightly zpterm
; nightly airmass term

print,'RUNNING SOLVE_TRANSPHOT ON  >> ',file,' <<'

reduxdir = '/data/smash/cp/red/photred/'
plotsdir = reduxdir+'stdred/plots/'

; Defaults
if n_elements(noreject) eq 0 then noreject=1
if n_elements(nosmith) eq 0 then nosmith=1
if n_elements(fixchipzpterm) eq 0 then fixchipzpterm=0
if n_elements(fitchipzpmjdterm) eq 0 then fitchipzpmjdterm=0
; Setting
print,'NOREJECT = ',keyword_set(noreject)
print,'NOSMITH = ',keyword_set(nosmith)
print,'FIXCHIPZPTERM = ',keyword_set(fixchipzpterm)
print,'FITCHIPZPMJDTERM = ',keyword_set(fitchipzpmjdterm)


; Load the data
print,'Loading the data'
arr = MRDFITS(file,1)
narr = n_elements(arr)

tags = tag_names(arr)
filter = strlowcase(tags[6])  ; this should be the filter name
colorname = repstr(strlowcase(tags[8]),'_','-')

; Suffix tag for plots
pstag = ''
;if keyword_set(noreject) then pstag='_norej'

; Trim stars
if n_elements(errlim) eq 0 then errlim = 0.05
;if filter eq 'u' then errlim=0.08
print,'Imposing error cut, ERR<=',strtrim(errlim,2),' mag'
gd = where(arr.err lt errlim,ngd)
arr = arr[gd]
narr = n_elements(arr)

; Remove the Smith data
if keyword_set(nosmith) then begin
  print,'#############################'
  print,'REMOVING SMITH DATA!!!!!'
  print,'#############################'
  bd = where(strmid(strtrim(arr.id,2),0,5) eq 'ugriz',nbd)
  if nbd gt 0 then remove,bd,arr
  narr = n_elements(arr)
endif

; u-band, removing very blue and very red stars
;  and removing the quadratic color shape
if filter eq 'u' then begin
  print,'Removing very blue/red stars'
  gd = where(arr.(8) ge 1.00 and arr.u_g lt 2.5,ngd)  ; u-g
  ;gd = where(arr.(8) ge 0.73 and arr.u_g lt 2.5,ngd)  ; u-g
  ;gd = where(arr.(8) ge 1.3 and arr.u_g lt 2.5,ngd)  ; u-g
  arr = arr[gd]
  narr = n_elements(arr)

  ; Removing quadratic color shape
  print,'Removing quadratic color shape'
  ucoef = [ -0.19919922d0, 0.27770109d0, -0.088840717d0]
  arr.mag -= poly(arr.(8),ucoef)
endif

; r-band, removing very red stars
if filter eq 'r' then begin
  print,'Removing very red stars'
  gd = where(arr.(8) lt 1.20,ngd) ; g-r
  arr = arr[gd]
  narr = n_elements(arr)
endif

; Adding expnum to ARR
if tag_exist(arr,'expnum') eq 0 then begin
  print,'Adding EXPNUM column to ARR'
  add_tag,arr,'EXPNUM','',arr
  dum = strsplitter(arr.frame,'_',/extract)
  expnum = reform(dum[0,*])
  arr.expnum = expnum
endif

; Loading the weather conditions table
conditions = importascii('smash_observing_conditions.txt',/header)
print,'Weather conditions'
print,'  NUM    MJD       DATE  PHOTOMETRIC     COMMENTS'
WRITECOL,-1,conditions.num,conditions.mjd,conditions.date,conditions.photometric,'   '+conditions.comments,fmt='(I5,I8,I12,I5,A-50)'


; 1.) Removing outlier exposures
;---------------------------------
; there is also an exposure outlier rejection in solve_transphot_night.pro
print,'' & print,'' & print,'1.) Removing outlier exposures'
ui = uniq(arr.expnum,sort(arr.expnum))
uexp = arr[ui].expnum
nexp = n_elements(uexp)
colcoef = robust_poly_fitq(arr.(8),arr.mag-arr.(6),1,/silent)
amcoef = robust_poly_fitq(arr.airmass,arr.mag-arr.(6),1,/silent)
expstr = replicate({expnum:'',night:0L,mjd:0L,exptime:0.0,airmass:0.0,nsources:0L,zeropt:0.0,residzeropt:0.0,residzeroptnsig:0.0},nexp)
; Get indices for each group
expsi = sort(arr.expnum)
expbrk = where(shift(arr[expsi].expnum,1) ne arr[expsi].expnum,nexpbrk)
explo = expbrk
exphi = [expbrk[1:nexpbrk-1]-1,n_elements(arr)-1]
for i=0,nexp-1 do begin
  ;MATCH,arr.expnum,uexp[i],ind1,ind2,/sort,count=nmatch
  ind = expsi[explo[i]:exphi[i]]
  nind = exphi[i]-explo[i]+1
  medoff = median([arr[ind].mag-arr[ind].(6)-arr[ind].(8)*colcoef[1]-arr[ind].airmass*amcoef[1]])
  expstr[i].expnum = uexp[i]
  ;expstr[i].night = arr[ind[0]].night
  expstr[i].mjd = arr[ind[0]].mjd
  expstr[i].exptime = arr[ind[0]].exptime
  expstr[i].airmass = median([arr[ind].airmass])
  expstr[i].nsources = nind
  expstr[i].zeropt = medoff
  ;print,i+1,uexp[i],medoff,format='(I5,A12,F10.4)'
endfor
; Compute medoff per MJD
uimjd = uniq(expstr.mjd,sort(expstr.mjd))
umjd = expstr[uimjd].mjd
nmjd = n_elements(umjd)
mjdzeropt = fltarr(nmjd)
mjdzpsig = fltarr(nmjd)
for i=0,nmjd-1 do begin
  ind = where(expstr.mjd eq umjd[i],nind)
  expstr[ind].night = i+1
  mjdzeropt[i] = median([expstr[ind].zeropt])
  mjdzpsig[i] = mad([expstr[ind].zeropt])
  expstr[ind].residzeropt = expstr[ind].zeropt - mjdzeropt[i]
endfor
; Measure sigma for each MJD
bdsig = where(mjdzpsig lt 1e-8,nbdsig,comp=gdsig)
zpsig = median(mjdzpsig[gdsig])
if nbdsig gt 0 then mjdzpsig[bdsig] = zpsig
; Measure deviation in Nsigma
;  constrain the sigma to [0.7,1.5]*median(sigma)
for i=0,nmjd-1 do begin
  ind = where(expstr.mjd eq umjd[i],nind)
  expstr[ind].residzeroptnsig = abs(expstr[ind].residzeropt) / (0.7*zpsig > mjdzpsig[i] < zpsig*1.5)
endfor
; Find bad exposures
bdexp = where(expstr.residzeroptnsig gt 5,nbdexp)
if nbdexp gt 0 then begin
  for i=0,nmjd-1 do begin
    match,umjd[i],expstr[bdexp].mjd,ind1,ind2,/sort,count=nmatch1
    if nmatch1 gt 0 then print,strtrim(umjd[i],2),'  ',strtrim(nmatch1,2),' outlier exposures: ',expstr[bdexp[ind2]].expnum
  endfor
  ;print,'Removing ',strtrim(nbdexp,2),' bad exposures: ',expstr[bdexp].expnum
endif
undefine,bdind
for i=0,nbdexp-1 do begin
  ind = expsi[explo[bdexp[i]]:exphi[bdexp[i]]]
  push,bdind,ind
endfor
if not keyword_set(noreject) then REMOVE,bdind,arr else print,'/NOREJECT.  Not rejecting any exposures.'
narr = n_elements(arr)

; 57097  huge airmass dependence, maybe only use low-airmass exps to
;          get zpterm
; 57099  ~4 expnums are way off the others

; Other bad exposures
badexp = ['00187877','00187878',$  ; there are ~10 chips that have problems, night 5+6, i-band
          '00422929','00423077','00423026','00422153','00422939','00422158']   ; 57097, 57099 bad data
nbadexp = n_elements(badexp)
for i=0,nbadexp-1 do begin
  MATCH,arr.expnum,badexp[i],ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 and not keyword_set(noreject) then begin
    remove,ind1,arr
    print,'Removing bad exposure: ',badexp[i]
  endif
endfor


; 2.) Run on all the nights fitting everything individually
;------------------------------------------------------------
print,'' & print,'2.) Run on all the nights fitting everything individually'
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr0,fitstr0,resid=resid0,expstr=expstr0,rejected=rejected0,/bootstrap,noreject=noreject
nmjds = n_elements(mstr0)
uimjd = uniq(mstr0.mjd,sort(mstr0.mjd))
umjd = mstr0[uimjd].mjd
nfitstr = n_elements(fitstr0)


; Set PHOTOMETRIC tag
;--------------------
;   2   56370  20130318   24080     0.1266,  zpterm high by 0.5 mag
;        log says there were intermittent clouds
;  40   57099  20150317   64949     0.6022,  zpterm WAY HIGH by 1.0 mag
;        log says non-photometric
;  50   57435  20160216    5579     0.0409,  zpterm high by 0.3 mag
;        log says thin cirrus, non-photometric
MATCH,mstr0.mjd,conditions.mjd,ind1,ind2
mstr0[ind1].photometric = conditions[ind2].photometric
for i=0,n_elements(conditions)-1 do begin
  MATCH,fitstr0.mjd,conditions[i].mjd,ind1,ind2,count=nmatch,/sort
  if nmatch gt 0 then fitstr0[ind1].photometric=conditions[i].photometric
endfor

; SHOULD I REMOVE THESE FROM "ARR" so they don't POLLUTE the rest???

; Flag "bad" solution with not Not enough stars
;--------------------------------
; 8   56664  20140106, only 16 stars,  SO FEW STARS, DON'T TRUST ANY OF THE RESULTS
;      maybe determine single global zpterm for the night
; 9   56665  20140107, only 51 stars,  not much better
bdmjd = where(mstr0.nstars lt 100,nbdmjd,comp=gdmjd,ncomp=ngdmjd)
if nbdmjd gt 0 then mstr0[bdmjd].badsoln = 1    ; not enough stars
if ngdmjd gt 0 then mstr0[gdmjd].badsoln = 0    ; enough stars
fitstr0.badsoln = 0
for i=0,nbdmjd-1 do begin
  MATCH,fitstr0.mjd,mstr0[bdmjd[i]].mjd,ind1,ind2,count=nmatch,/sort
  fitstr0[ind1].badsoln = 1
endfor


; 3.) Determine median zpterm and color terms per chip
;------------------------------------------------------
print,'' & print,'3.) Determining median zpterm and colterm per chip'
ui = uniq(fitstr0.chip,sort(fitstr0.chip))  ; unique chips
uchips = fitstr0[ui].chip
nuchips = n_elements(uchips)
; Measure median and RMS per chip
if keyword_set(fitchipzpmjdterm) then begin
  chipstr = replicate({chip:0L,zpterm:99.0,zptermerr:99.0,zpmjdterm:99.0,zpmjdtermerr:99.0,zptermsig:99.0,colterm:99.0,coltermerr:99.0,$
                       coltermsig:99.0,rms:99.0,medrms:99.0},nuchips)
endif else begin
  chipstr = replicate({chip:0L,zpterm:99.0,zptermerr:99.0,zptermsig:99.0,colterm:99.0,coltermerr:99.0,coltermsig:99.0,rms:99.0,medrms:99.0},nuchips)
endelse
fitstr0.relzpterm = fitstr0.zpterm        ; initialize with zpterm
print,' CHIP   ZPTERM   ZPTERMERR  COLTERM  COLTERMER    RMS'
for i=0,nuchips-1 do begin
  chipstr[i].chip = uchips[i]
  MATCH,fitstr0.chip,uchips[i],ind,ind0,count=nind,/sort
  if nind gt 0 then begin
    ; Remove median zpterm for each night
    MATCH,mstr0.mjd,fitstr0[ind].mjd,ind1,ind2
    fitstr0[ind[ind2]].relzpterm -= mstr0[ind1].zpterm  ; remove the nightly median zeropoint
    ; Only use "good" value to estimate median and scatter
    gdphot = where(fitstr0[ind].photometric eq 1 and fitstr0[ind].badsoln eq 0,ngdphot)
    ; Fit linear temporal variations in chip relative zeropoint term
    if keyword_set(fitchipzpmjdterm) then begin
      zpmjdcoef1 = robust_poly_fitq(fitstr0[ind[gdphot]].mjd,fitstr0[ind[gdphot]].relzpterm,1)
      ydiff1 = fitstr0[ind[gdphot]].relzpterm-poly(fitstr0[ind[gdphot]].mjd,zpmjdcoef1)
      relzptermsig1 = mad(ydiff1)
      gd2 = where(abs(ydiff1) lt 5*relzptermsig1,ngd2)
      zpmjdcoef = dln_poly_fit(fitstr0[ind[gdphot[gd2]]].mjd,fitstr0[ind[gdphot[gd2]]].relzpterm,1,$
                               measure_errors=fitstr0[ind[gdphot[gd2]]].zptermerr,sigma=zpmjdcoeferr,yerror=yerror,status=status,yfit=yfit1,/bootstrap)
      ydiff = fitstr0[ind[gdphot]].relzpterm-poly(fitstr0[ind[gdphot]].mjd,zpmjdcoef)
      relzptermsig = mad(ydiff)
      chipstr[i].zpterm = zpmjdcoef[0]
      chipstr[i].zptermerr = zpmjdcoeferr[0]
      chipstr[i].zpmjdterm = zpmjdcoef[1]
      chipstr[i].zpmjdtermerr = zpmjdcoeferr[1]
      chipstr[i].zptermsig = relzptermsig
      fitstr0[ind].chiprelzpterm = fitstr0[ind].relzpterm - (chipstr[i].zpterm+fitstr0[ind].mjd*chipstr[i].zpmjdterm) ; relative to night and chip median
    endif else begin
      chipstr[i].zpterm = median(fitstr0[ind[gdphot]].relzpterm)
      chipstr[i].zptermsig = mad(fitstr0[ind[gdphot]].relzpterm)
      chipstr[i].zptermerr = chipstr[i].zptermsig / sqrt(ngdphot)
      fitstr0[ind].chiprelzpterm = fitstr0[ind].relzpterm - chipstr[i].zpterm   ; relative to night and chip medians
    endelse
    chipstr[i].colterm = median([fitstr0[ind[gdphot]].colterm])
    chipstr[i].coltermsig = mad([fitstr0[ind[gdphot]].colterm])
    chipstr[i].coltermerr = chipstr[i].coltermsig / sqrt(ngdphot)
    totresid = fitstr0[ind[gdphot]].nbrt*fitstr0[ind[gdphot]].brtrms^2
    rms = sqrt(mean(total(totresid)/total(fitstr0[ind[gdphot]].nbrt)))  ; brt rms
    chipstr[i].rms = rms
    chipstr[i].medrms = median([fitstr0[ind[gdphot]].rms])
    fitstr0[ind].chiprelcolterm = fitstr0[ind].colterm - chipstr[i].colterm
    print,chipstr[i].chip,chipstr[i].zpterm,chipstr[i].zptermerr,chipstr[i].colterm,chipstr[i].coltermerr,chipstr[i].rms,format='(I5,5F10.5)'
  endif
endfor

;stop

; Making some figures
print,'Making some figures'
setdisp,/silent
!p.font = 0
; Relative ZPTERM vs. chip
psfile = plotsdir+'transphot_'+filter+'_zpterm_chip'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
rnd = randomu(seed,nfitstr)*0.6-0.3
yr = [-0.1,0.1]
gd = where(fitstr0.photometric eq 1 and fitstr0.badsoln eq 0,ngd)
plotc,fitstr0[gd].chip+rnd[gd],fitstr0[gd].relzpterm,fitstr0[gd].mjd,ps=8,sym=0.5,xr=[0,63],yr=yr,xs=1,ys=1,xtit='Chip',ytit='Relative ZPTERM',$
      tit=filter+'-band zero-points   (color-coded by MJD)',pos=[0.08,0.08,0.98,0.40],colpos=[0.08,0.47,0.98,0.48]
oplot,chipstr.chip,chipstr.zpterm,ps=-1
xyouts,2,yr[1]-0.08*range(yr),'RMS = '+stringize(mad(fitstr0[gd].chiprelzpterm),ndec=4),align=0,charsize=1.1,charthick=4
plotc,fitstr0[gd].chip+rnd[gd],fitstr0[gd].relzpterm,fitstr0[gd].zptermerr,ps=8,sym=0.5,xr=[0,63],yr=yr,xs=1,ys=1,max=0.01,xtit='Chip',ytit='Relative ZPTERM',$
      tit=filter+'-band zero-points   (color-coded by ZPTERMERR)',pos=[0.08,0.58,0.98,0.90],colpos=[0.08,0.97,0.98,0.98],/noerase
oplot,chipstr.chip,chipstr.zpterm,ps=-1
xyouts,2,yr[1]-0.08*range(yr),'RMS = '+stringize(mad(fitstr0[gd].chiprelzpterm),ndec=4),align=0,charsize=1.1,charthick=4
ps_close
ps2png,psfile+'.eps',/eps
; COLTERM vs. chip
psfile = plotsdir+'transphot_'+filter+'_colterm_chip'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
rnd = randomu(seed,nfitstr)*0.6-0.3
;yr = [-0.3,0.3]
yr = [min(-3*chipstr.coltermsig+chipstr.colterm),max(3*chipstr.coltermsig+chipstr.colterm)]
plotc,fitstr0[gd].chip+rnd[gd],fitstr0[gd].colterm,fitstr0[gd].mjd,ps=8,sym=0.5,xr=[0,63],yr=yr,xs=1,ys=1,xtit='Chip',ytit='COLTERM',$
      tit=filter+'-band color terms   (color-coded by MJD)',pos=[0.08,0.08,0.98,0.40],colpos=[0.08,0.47,0.98,0.48]
oplot,chipstr.chip,chipstr.colterm,ps=-1
xyouts,2,yr[1]-0.08*range(yr),'RMS = '+stringize(mad(fitstr0[gd].chiprelcolterm),ndec=4),align=0,charsize=1.1,charthick=4
plotc,fitstr0[gd].chip+rnd[gd],fitstr0[gd].colterm,fitstr0[gd].coltermerr,ps=8,sym=0.5,xr=[0,63],yr=yr,xs=1,ys=1,max=0.01,xtit='Chip',ytit='COLTERM',$
      tit=filter+'-band color terms   (color-coded by COLTERMERR)',pos=[0.08,0.58,0.98,0.90],colpos=[0.08,0.97,0.98,0.98],/noerase
oplot,chipstr.chip,chipstr.colterm,ps=-1
xyouts,2,yr[1]-0.08*range(yr),'RMS = '+stringize(mad(fitstr0[gd].chiprelcolterm),ndec=4),align=0,charsize=1.1,charthick=4
ps_close
ps2png,psfile+'.eps',/eps
; CHIPRELZPTERM vs. NIGHTNUM
psfile = plotsdir+'transphot_'+filter+'_chiprelzpterm_nightnum'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
rnd = randomu(seed,nfitstr)*0.6-0.3
gg = where(abs(fitstr0.chiprelzpterm) lt 2)
yr = [ min(fitstr0[gg].chiprelzpterm) > (-3*stddev(fitstr0[gg].chiprelzpterm)),$
       max(fitstr0[gg].chiprelzpterm) < 3*stddev(fitstr0[gg].chiprelzpterm)]
plotc,fitstr0.nightnum+rnd,fitstr0.chiprelzpterm,fitstr0.chip,ps=8,sym=0.5,xr=[0,53],yr=yr,xs=1,ys=1,xtit='Night Number',ytit='Relative ZPTERM',$
      tit=filter+'-band relative zero-point term   (color-coded by CHIP)',pos=[0.08,0.08,0.98,0.40],colpos=[0.08,0.47,0.98,0.48]
plotc,fitstr0.nightnum+rnd,fitstr0.chiprelzpterm,fitstr0.zptermerr,ps=8,sym=0.5,xr=[0,53],yr=yr,xs=1,ys=1,max=0.01,xtit='Night Number',ytit='Relative ZPTERM',$
      tit=filter+'-band relative zero-point term   (color-coded by ZPTERMERR)',pos=[0.08,0.58,0.98,0.90],colpos=[0.08,0.97,0.98,0.98],/noerase
ps_close
ps2png,psfile+'.eps',/eps
; CHIPRELCOLTERM vs. NIGHTNUM
psfile = plotsdir+'transphot_'+filter+'_chiprelcolterm_nightnum'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
rnd = randomu(seed,nfitstr)*0.6-0.3
gg = where(abs(fitstr0.chiprelcolterm) lt 2)
yr = [ min(fitstr0[gg].chiprelcolterm) > (-5*stddev(fitstr0[gg].chiprelcolterm)),$
       max(fitstr0[gg].chiprelcolterm) < 5*stddev(fitstr0[gg].chiprelcolterm)]
plotc,fitstr0.nightnum+rnd,fitstr0.chiprelcolterm,fitstr0.chip,ps=8,sym=0.5,xr=[0,53],yr=yr,xs=1,ys=1,xtit='Night Number',ytit='Relative COLTERM',$
      tit=filter+'-band relative color term   (color-coded by CHIP)',pos=[0.08,0.08,0.98,0.40],colpos=[0.08,0.47,0.98,0.48]
plotc,fitstr0.nightnum+rnd,fitstr0.chiprelcolterm,fitstr0.coltermerr,ps=8,sym=0.5,xr=[0,53],yr=yr,xs=1,ys=1,max=0.01,xtit='Night Number',ytit='Relative COLTERM',$
      tit=filter+'-band relative color term   (color-coded by COLTERMERR)',pos=[0.08,0.58,0.98,0.90],colpos=[0.08,0.97,0.98,0.98],/noerase
ps_close
ps2png,psfile+'.eps',/eps



; 4.) Fixing colterm, Redetermine zpterm and airmass term
;---------------------------------------------------------
print,'' & print,'4.) Fixing colterm.  Redetermining zpterm and amterm'
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolr,fitstr_fixcolr,fixcolr=chipstr,resid=resid_fixcolr,rejected=rejected_fixcolr,$
                          /bootstrap,noreject=noreject

; Set PHOTOMETRIC tag
MATCH,mstr_fixcolr.mjd,conditions.mjd,ind1,ind2
mstr_fixcolr[ind1].photometric = conditions[ind2].photometric
for i=0,n_elements(conditions)-1 do begin
  MATCH,fitstr_fixcolr.mjd,conditions[i].mjd,ind1,ind2,count=nmatch,/sort
  if nmatch gt 0 then fitstr_fixcolr[ind1].photometric=conditions[i].photometric
endfor

; Flag "bad" solution with not Not enough stars
bdmjd = where(mstr_fixcolr.nstars lt 100,nbdmjd,comp=gdmjd,ncomp=ngdmjd)
if nbdmjd gt 0 then mstr_fixcolr[bdmjd].badsoln = 1    ; not enough stars
if ngdmjd gt 0 then mstr_fixcolr[gdmjd].badsoln = 0    ; enough stars
fitstr_fixcolr.badsoln = 0
for i=0,nbdmjd-1 do begin
  MATCH,fitstr_fixcolr.mjd,mstr_fixcolr[bdmjd[i]].mjd,ind1,ind2,count=nmatch,/sort
  fitstr_fixcolr[ind1].badsoln = 1
endfor


; 5.) Re-Determine median zpterm per chip
;------------------------------------------
print,'' & print,'5.) Re-determine median zpterm per chip'
; Measure median chip zpterm again
if keyword_set(fitchipzpmjdterm) then begin
  chipstr_fixcolr = replicate({chip:0L,zpterm:99.0,zptermerr:99.0,zpmjdterm:99.0,zpmjdtermerr:99.0,zptermsig:99.0,colterm:99.0,coltermerr:99.0,$
                               coltermsig:99.0,rms:99.0,medrms:99.0},nuchips)
endif else begin
  chipstr_fixcolr = replicate({chip:0L,zpterm:99.0,zptermerr:99.0,zptermsig:99.0,colterm:99.0,coltermerr:99.0,coltermsig:99.0,rms:99.0,medrms:99.0},nuchips)
endelse
fitstr_fixcolr.relzpterm = fitstr_fixcolr.zpterm        ; initialize with zpterm
print,' CHIP   ZPTERM   ZPTERMERR  COLTERM  COLTERMER    RMS'
for i=0,nuchips-1 do begin
  chipstr_fixcolr[i].chip = uchips[i]
  MATCH,fitstr_fixcolr.chip,uchips[i],ind,ind0,count=nind,/sort
  if nind gt 0 then begin
    ; Remove median zpterm for each night
    MATCH,mstr_fixcolr.mjd,fitstr_fixcolr[ind].mjd,ind1,ind2
    fitstr_fixcolr[ind[ind2]].relzpterm -= mstr_fixcolr[ind1].zpterm  ; remove the nightly median zeropoint
    ; Only use "good" value to estimate median and scatter
    gdphot = where(fitstr_fixcolr[ind].photometric eq 1 and fitstr_fixcolr[ind].badsoln eq 0,ngdphot)
    ; Fit linear temporal variations in chip relative zeropoint term
    if keyword_set(fitchipzpmjdterm) then begin
      zpmjdcoef1 = robust_poly_fitq(fitstr_fixcolr[ind[gdphot]].mjd,fitstr_fixcolr[ind[gdphot]].relzpterm,1)
      ydiff1 = fitstr_fixcolr[ind[gdphot]].relzpterm-poly(fitstr_fixcolr[ind[gdphot]].mjd,zpmjdcoef1)
      relzptermsig1 = mad(ydiff1)
      gd2 = where(abs(ydiff1) lt 5*relzptermsig1,ngd2)
      zpmjdcoef = dln_poly_fit(fitstr_fixcolr[ind[gdphot[gd2]]].mjd,fitstr_fixcolr[ind[gdphot[gd2]]].relzpterm,1,$
                               measure_errors=fitstr_fixcolr[ind[gdphot[gd2]]].zptermerr,sigma=zpmjdcoeferr,yerror=yerror,status=status,yfit=yfit1,/bootstrap)
      ydiff = fitstr_fixcolr[ind[gdphot]].relzpterm-poly(fitstr_fixcolr[ind[gdphot]].mjd,zpmjdcoef)
      relzptermsig = mad(ydiff)
      chipstr_fixcolr[i].zpterm = zpmjdcoef[0]
      chipstr_fixcolr[i].zptermerr = zpmjdcoeferr[0]
      chipstr_fixcolr[i].zpmjdterm = zpmjdcoef[1]
      chipstr_fixcolr[i].zpmjdtermerr = zpmjdcoeferr[1]
      chipstr_fixcolr[i].zptermsig = relzptermsig
      fitstr_fixcolr[ind].chiprelzpterm = fitstr_fixcolr[ind].relzpterm - (chipstr_fixcolr[i].zpterm+fitstr_fixcolr[ind].mjd*chipstr_fixcolr[i].zpmjdterm) ; relative to night and chip median
    endif else begin
      chipstr_fixcolr[i].zpterm = median(fitstr_fixcolr[ind[gdphot]].relzpterm)
      chipstr_fixcolr[i].zptermsig = mad(fitstr_fixcolr[ind[gdphot]].relzpterm)
      chipstr_fixcolr[i].zptermerr = chipstr_fixcolr[i].zptermsig / sqrt(ngdphot)
      fitstr_fixcolr[ind].chiprelzpterm = fitstr_fixcolr[ind].relzpterm - chipstr_fixcolr[i].zpterm   ; relative to night and chip medians
    endelse
    chipstr_fixcolr[i].colterm = median([fitstr_fixcolr[ind[gdphot]].colterm])
    chipstr_fixcolr[i].coltermsig = mad([fitstr_fixcolr[ind[gdphot]].colterm])
    chipstr_fixcolr[i].coltermerr = chipstr_fixcolr[i].coltermsig / sqrt(ngdphot)
    totresid = fitstr_fixcolr[ind[gdphot]].nbrt*fitstr_fixcolr[ind[gdphot]].brtrms^2
    rms = sqrt(mean(total(totresid)/total(fitstr_fixcolr[ind[gdphot]].nbrt)))  ; brt rms
    chipstr_fixcolr[i].rms = rms
    chipstr_fixcolr[i].medrms = median([fitstr_fixcolr[ind[gdphot]].rms])
    fitstr_fixcolr[ind].chiprelcolterm = fitstr_fixcolr[ind].colterm - chipstr_fixcolr[i].colterm
    print,chipstr_fixcolr[i].chip,chipstr_fixcolr[i].zpterm,chipstr_fixcolr[i].zptermerr,chipstr_fixcolr[i].colterm,$
          chipstr_fixcolr[i].coltermerr,chipstr_fixcolr[i].rms,format='(I5,5F10.5)'
    ;; Checking for temporal trends in relative zpterm
    ;zpcoef = robust_poly_fitq(fitstr_fixcolr[ind[gdphot]].mjd,fitstr_fixcolr[ind[gdphot]].chiprelzpterm,1)
    ;plotc,fitstr_fixcolr[ind[gdphot]].mjd,fitstr_fixcolr[ind[gdphot]].chiprelzpterm,ps=8,yerr=fitstr_fixcolr[ind[gdphot]].zptermerr,xs=1,tit='chip = '+strtrim(uchips[i],2)
    ;xx = scale_vector(findgen(100),min(umjd),max(umjd))
    ;oplot,xx,poly(xx,zpcoef),co=250
    ;print,i+1,zpcoef[1],zpcoef[1]*range(umjd)
    ;wait,1
  endif
endfor
; COLTERMERR = 0.0 because it was fixed

;stop


; 6.) Redetermine NIGHTLY zpterm and amterm with chip-level relative zpterm 
;      and colterm fixed
;---------------------------------------------------------------------------
print,'' & print,'6.) Fixing chip-level relative zpterm and colterm.  Redetermining NIGHTLY zpterm and amterm'
if keyword_set(fixchipzpterm) then begin
  SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolzp,fitstr_fixcolzp,fixcolr=chipstr,fixchipzp=chipstr_fixcolr,resid=resid_fixcolzp,$
                            rejected=rejected_fixcolzp,/bootstrap,noreject=noreject
endif else begin
  print,'NOT FIXING CHIP-LEVEL ZERO-POINT TERMS'
  mstr_fixcolzp = mstr_fixcolr
  fitstr_fixcolzp = fitstr_fixcolr
  resid_fixcolzp = resid_fixcolr
  rejected_fixcolzp = rejected_fixcolr
endelse

; Reset PHOTOMETRIC tag
MATCH,mstr_fixcolzp.mjd,conditions.mjd,ind1,ind2
mstr_fixcolzp[ind1].photometric = conditions[ind2].photometric
for i=0,n_elements(conditions)-1 do begin
  MATCH,fitstr_fixcolzp.mjd,conditions[i].mjd,ind1,ind2,count=nmatch,/sort
  if nmatch gt 0 then fitstr_fixcolzp[ind1].photometric=conditions[i].photometric
endfor

; Flag "bad" solution with not Not enough stars
bdmjd = where(mstr_fixcolzp.nstars lt 100,nbdmjd,comp=gdmjd,ncomp=ngdmjd)
if nbdmjd gt 0 then mstr_fixcolzp[bdmjd].badsoln = 1    ; not enough stars
if ngdmjd gt 0 then mstr_fixcolzp[gdmjd].badsoln = 0    ; enough stars
fitstr_fixcolzp.badsoln = 0
for i=0,nbdmjd-1 do begin
  MATCH,fitstr_fixcolzp.mjd,mstr_fixcolzp[bdmjd[i]].mjd,ind1,ind2,count=nmatch,/sort
  fitstr_fixcolzp[ind1].badsoln = 1
endfor

;stop

; 7.) Averaging airmass terms
;----------------------------
; Bad for non-photometric nights 2, 40, 50
; Bad for bad solution nights 8, 9
; Bad for 38
; 6, 15, 22, 24 are a bit high

; i-band: 6, 14, 23, 37, 39, 49  these are off
;      56510 (clear/phot), 56681 (cloudy), 56702 (clear), 57097 (a few
;      cloud, crazy standard data), 57099 (nonphot), 57435 (some clouds)
; most are non-phot or bad soln, but only one that's "good" is off 56510 (clear/phot)

averageamterm:

; the airmass terms DO get more consistent/better after the various
; parameter "fixings"

; Fix airmass terms for nights with low airmass range
print,'' & print,'7.) Averaging airmass terms'
ntstr_fixcolzp = mstr_fixcolzp
add_tag,ntstr_fixcolzp,'amavgflag',0,ntstr_fixcolzp  ; amterm averaged
add_tag,ntstr_fixcolzp,'namavg',0,ntstr_fixcolzp  ; amterm averaged
print,'  NUM   MJD    AMTERM   AMTERMSIG FLAG  NAVG   COMMENT'
for i=0,nmjds-1 do begin
  ; Only look at photometric nights
  if mstr_fixcolzp[i].photometric eq 1 then begin
    ; Low airmass range or bad solution
    if mstr_fixcolzp[i].amrange lt 0.35 or mstr_fixcolzp[i].badsoln eq 1 then begin
      ; Get nights close in time with good airmass range
      si = sort(abs(mstr_fixcolzp.mjd-mstr_fixcolzp[i].mjd))
      gdnei = where(mstr_fixcolzp[si].amrange gt 0.35 and mstr_fixcolzp[si].badsoln eq 0,ngdnei)
      nei = si[gdnei[0:3]]
      ; Weights, time difference and amterm uncertainty
      wt = 1.0 / abs(mstr_fixcolzp[nei].mjd-mstr_fixcolzp[i].mjd)
      wt *= 1.0/mstr_fixcolzp[nei].amtermsig^2
      wt /= total(wt)   ; normalize
      ; Computed weighted mean
      xmn = total(mstr_fixcolzp[nei].amterm * wt) / total(wt)
      ; error of weighted mean
      ;  https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
      ;  "Statistical Properties" section
      xsig = sqrt( total( wt^2 * mstr_fixcolzp[nei].amtermsig^2 ) )
      ; Stuff it in the structure
      ntstr_fixcolzp[i].amterm = xmn
      ntstr_fixcolzp[i].amtermsig = xsig
      ntstr_fixcolzp[i].amavgflag = 1      ; we averaged it the AMTERM values
      ntstr_fixcolzp[i].namavg = n_elements(nei)    ; number averaged
      comment = '  low amrange'

    ; Refit AMTERM with neighboring nights's data
    endif else begin ; low amrange

      ; Get nights close in time with good airmass range
      diff_mjd = abs(mstr_fixcolzp.mjd-mstr_fixcolzp[i].mjd)
      gdnei1 = where(mstr_fixcolzp.amrange gt 0.35 and diff_mjd le 30 and mstr_fixcolzp.badsoln eq 0 and mstr_fixcolzp.photometric eq 1,ngdnei1)

      ; Enough to average
      if ngdnei1 gt 1 then begin

        ; Take the nearest 4 neighbor nights (and the current one)
        si = sort(diff_mjd[gdnei1])
        gdnei = gdnei1[si[0:(4<(ngdnei1-1))]]
        ngdnei = n_elements(gdnei)

        undefine,useind,resid
        for j=0,ngdnei-1 do begin
          MATCH,arr.mjd,mstr_fixcolzp[gdnei[j]].mjd,ind1,ind2,count=nmatch,/sort
          if nmatch gt 0 then begin
            push,useind,ind1
            ; Add AMTERM back into residuals, but subtract 1
            ;   so that airmass=1 points will be around zero
            push,resid,resid_fixcolzp[ind1]+(arr[ind1].airmass-1)*mstr_fixcolzp[gdnei[j]].amterm
          endif
        endfor
        coef1 = robust_poly_fitq(arr[useind].airmass,resid,1)
        sig1 = mad(resid-poly(arr[useind].airmass,coef1))
        ; cut outliers and stars previously rejected
        gdind = where(abs(resid-poly(arr[useind].airmass,coef1)) lt 4.0*sig1 and rejected_fixcolzp[useind] eq 0,ngdind)
        coef = dln_poly_fit(arr[useind[gdind]].airmass,resid[gdind],1,measure_errors=arr[useind[gdind]].err,sigma=coeferr,yerror=yerror,status=status,/bootstrap)
        ntstr_fixcolzp[i].amterm = coef[1]
        ntstr_fixcolzp[i].amtermsig = coeferr[1]
        ntstr_fixcolzp[i].amavgflag = 2      ; we determine AMTERM from actual points
        ntstr_fixcolzp[i].namavg = ngdnei    ; number averaged
        comment = ''
      endif else comment='  not enough nearby nights to average'
    endelse  ; average

  ; Non-photometric
  ;   Do NOT determine/set airmass or zpterm for NON-PHOTOMETRIC nights
  endif else begin
    ntstr_fixcolzp[i].amterm = 0.0
    ntstr_fixcolzp[i].amtermsig = 0.0
    comment = '  non-photometric'
  endelse
  print,i+1,ntstr_fixcolzp[i].mjd,ntstr_fixcolzp[i].amterm,ntstr_fixcolzp[i].amtermsig,ntstr_fixcolzp[i].amavgflag,$
        ntstr_fixcolzp[i].namavg,comment,format='(I5,I7,2F10.5,I5,I5,A-30)'
endfor

;plotc,mstr_fixcolzp.nightnum,mstr_fixcolzp.amterm,mstr_fixcolzp.amrange,ps=1
;oplot,mstr_fixcolr.nightnum,mstr_fixcolr.amterm,ps=4,co=250
;oplot,mstr0.nightnum,mstr0.amterm,ps=6,co=150
;
;plotc,mstr_fixcolzp.mjd,mstr_fixcolzp.amterm,mstr_fixcolzp.amrange,ps=1,xs=1
;oplot,mstr_fixcolr.mjd,mstr_fixcolr.amterm,ps=4,co=250
;oplot,mstr0.mjd,mstr0.amterm,ps=6,co=150


; Make a figure of AMTERM
psfile = plotsdir+'transphot_'+filter+'_amterm'+pstag
ps_open,psfile,/color,thick=4,/encap
device,/inches,xsize=9.5,ysize=9.5
gphot = where(ntstr_fixcolzp.photometric eq 1,ngphot)
xr = minmax(ntstr_fixcolzp[gphot].mjd)
xr = [xr[0]-10,xr[1]+10]
plotc,ntstr_fixcolzp[gphot].mjd,ntstr_fixcolzp[gphot].amterm,ntstr_fixcolzp[gphot].amrange,ps=7,sym=1.5,pos=[0.08,0.52,0.97,0.91],colpos=[0.08,0.97,0.97,0.99],$
      xr=xr,xs=1,yerr=ntstr_fixcolzp[gphot].amtermsig,xtit='MJD',ytit='AMTERM',tit='AMTERM (photometric nights, color-coded by airmass range)'
oplot,mstr_fixcolzp[gphot].mjd,mstr_fixcolzp[gphot].amterm,ps=1
al_legend,['Original value','Averaged value'],psym=[1,7],color=[0,250],/top,/left,charsize=1.0
xr2 = minmax(ntstr_fixcolzp[gphot].nightnum)
xr2 = [xr2[0]-1,xr2[1]+1]
plotc,ntstr_fixcolzp[gphot].nightnum,ntstr_fixcolzp[gphot].amterm,ntstr_fixcolzp[gphot].amrange,ps=7,sym=1.5,pos=[0.08,0.07,0.97,0.46],/nocolorbar,/noerase,$
      xr=xr2,xs=1,yerr=ntstr_fixcolzp[gphot].amtermsig,xtit='Night Number',ytit='AMTERM'
oplot,mstr_fixcolzp[gphot].nightnum,mstr_fixcolzp[gphot].amterm,ps=1
ps_close
ps2png,psfile+'.eps',/eps

; Combine filter-level summary plots
; transphot_g_zpterm_chip.png
; transphot_g_colterm_chip.png
; transphot_g_chiprelzpterm_nightnum.png
; transphot_g_chiprelcolterm_nightnum.png
; transphot_g_amterm.png
cd,current=curdir
cd,plotsdir
figfiles = 'transphot_'+filter+'_'+['zpterm_chip','colterm_chip','chiprelzpterm_nightnum','chiprelcolterm_nightnum','amterm']
for i=0,n_elements(figfiles)-1 do spawn,['epstopdf',figfiles[i]+'.eps'],/noshell
; Combine
print,'Writing combined filter-level summary plots to transphot_'+filter+'_comb'+pstag+'.pdf'
cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=transphot_'+filter+'_comb'+pstag+'.pdf '
cmd += strjoin(figfiles+'.pdf',' ')
spawn,cmd,out,errout
cd,curdir


; 8.) Redetermine NIGHTLY zpterm holding everything else fixed
;--------------------------------------------------------------
print,'' & print,'8.) Redetermine NIGHTLY zpterm holding everything else fixed'
if keyword_set(fixchipzpterm) then begin
  SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolzpam,fitstr_fixcolzpam,fixcolr=chipstr,fixchipzp=chipstr_fixcolr,$
                          fixam=ntstr_fixcolzp,resid=resid_fixcolzpam,rejected=rejected_fixcolzpam,/bootstrap,noreject=noreject
endif else begin
  print,'NOT FIXING CHIP-LEVEL ZERO-POINT TERMS'
  SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolzpam,fitstr_fixcolzpam,fixcolr=chipstr,$
                          fixam=ntstr_fixcolzp,resid=resid_fixcolzpam,rejected=rejected_fixcolzpam,/bootstrap,noreject=noreject
endelse

; Reset PHOTOMETRIC tag
MATCH,mstr_fixcolzpam.mjd,conditions.mjd,ind1,ind2
mstr_fixcolzpam[ind1].photometric = conditions[ind2].photometric
for i=0,n_elements(conditions)-1 do begin
  MATCH,fitstr_fixcolzpam.mjd,conditions[i].mjd,ind1,ind2,count=nmatch,/sort
  if nmatch gt 0 then fitstr_fixcolzpam[ind1].photometric=conditions[i].photometric
endfor

; Flag "bad" solution with Not enough stars
bdmjd = where(mstr_fixcolzpam.nstars lt 100,nbdmjd,comp=gdmjd,ncomp=ngdmjd)
if nbdmjd gt 0 then mstr_fixcolzpam[bdmjd].badsoln = 1    ; not enough stars
if ngdmjd gt 0 then mstr_fixcolzpam[gdmjd].badsoln = 0    ; enough stars
fitstr_fixcolzpam.badsoln = 0
for i=0,nbdmjd-1 do begin
  MATCH,fitstr_fixcolzpam.mjd,mstr_fixcolzpam[bdmjd[i]].mjd,ind1,ind2,count=nmatch,/sort
  fitstr_fixcolzpam[ind1].badsoln = 1
endfor
bdfitstr = where(fitstr_fixcolzpam.nstars lt 3 and fitstr_fixcolzpam.zpterm gt 100,nbdfitstr)
if nbdfitstr gt 0 then fitstr_fixcolzpam[bdfitstr].badsoln = 1

; Make final residual figures
print,'Making final residual plots'
undefine,residfigs
for i=0,nmjds-1 do begin

  MATCH,arr.mjd,mstr_fixcolzpam[i].mjd,ind1,ind2,count=nmatch,/sort
  gd = where(rejected_fixcolzp[ind1] eq 0,ngd)
  if ngd gt 0 then ind = ind1[gd] else ind=-1
  nind = ngd

  ; Make resid plots, one per night
  ;  resid vs. color, resid vs. airmass, resid vs. chip number

  setdisp,/silent
  psfile = plotsdir+'transphot_'+filter+'_'+strtrim(mstr_fixcolzpam[i].mjd,2)+'_resid'+pstag
  ps_open,psfile,/color,thick=4,/encap
  device,/inches,xsize=9.5,ysize=9.5
  x0 = 0.08 & x1=0.97
  y0 = 0.05 & dy = 0.30
  sig = mad(resid_fixcolzpam[ind])
  ;yr = minmax(resid_fixcolzpam[ind])
  yr = [-1,1]*(3.5*sig > 0.2)
  yr = [yr[0] > (-1.0), yr[1] < 1.0]
  zr = [0.0, max(arr.err)] 
  psym = 3
  symsize = 1.0
  if nind lt 10000 then begin
    psym = 8
    symsize = 0.1
  endif
  if nind lt 1000 then begin
    psym = 8
    symsize = 0.4
  endif
  ; Resid vs. color
  xr = minmax(arr[ind].(8))
  xr = [-0.5,2.0 < (xr[1]+0.05) > 3.5]
  plotc,[arr[ind].(8)],[resid_fixcolzpam[ind]],[arr[ind].err],ps=psym,symsize=symsize,xr=xr,yr=yr,xs=1,ys=1,xtit=colorname,ytit='Residuals',$
        pos=[x0,y0,x1,dy],colpos=[x0,0.97,x1,0.99],xticklen=0.04,min=zr[0],max=zr[1]
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'<COLTERM> = '+stringize(mstr_fixcolzpam[i].colterm,ndec=4)+'+/-'+$
         stringize(mstr_fixcolzpam[i].coltermsig,ndec=4),align=0,charsize=1.3
  oplot,xr,[0,0],linestyle=2
  ; Resid vs. airmass
  xr = minmax(arr[ind].airmass)
  xr = [0.9,(xr[1]+0.05) > 2.0]
  plotc,[arr[ind].airmass],[resid_fixcolzpam[ind]],[arr[ind].err],ps=psym,symsize=symsize,xr=xr,yr=yr,xs=1,ys=1,xtit='Airmass',ytit='Residuals',$
        pos=[x0,y0+dy,x1,2*dy],/nocolorbar,xticklen=0.04,/noerase,min=zr[0],max=zr[1]
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'AMTERM = '+stringize(mstr_fixcolzpam[i].amterm,ndec=4)+'+/'+$
         stringize(ntstr_fixcolzp[i].amtermsig,ndec=4),align=0,charsize=1.3
  if mstr_fixcolzpam[i].photometric eq 1 then photcom='PHOTOMETRIC' else photcom='NON-PHOTOMETRIC'
  if mstr_fixcolzpam[i].photometric eq 1 then col=0 else col=250
  xyouts,xr[1]-0.03*range(xr),yr[1]-0.12*range(yr),photcom,align=1,charsize=1.3,color=col
  oplot,xr,[0,0],linestyle=2
  ; Resid vs. chip number
  if nind gt 0 then rnd = randomu(seed,nind)*0.6-0.3 else rnd=0.0
  xr = [0,63]
  plotc,[arr[ind].chip+rnd],[resid_fixcolzpam[ind]],[arr[ind].err],ps=psym,symsize=symsize,xr=xr,yr=yr,xs=1,ys=1,xtit='Chip',ytit='Residuals',/noerase,$
        xticklen=0.04,tit=strtrim(mstr_fixcolzpam[i].mjd,2)+' - Residuals vs. color/airmass/chip (color-coded by ERR)',pos=[x0,y0+2*dy,x1,3*dy],$
        min=zr[0],max=zr[1],/nocolorbar
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'ZPTERM = '+stringize(mstr_fixcolzpam[i].zpterm,ndec=4)+'+/-'+$
         stringize(mstr_fixcolzpam[i].zptermsig,ndec=4),align=0,charsize=1.3
  xyouts,xr[1]-0.03*range(xr),yr[1]-0.12*range(yr),'RMS = '+stringize(mstr_fixcolzpam[i].rms,ndec=4),align=1,charsize=1.3
  oplot,xr,[0,0],linestyle=2

  ps_close
  ps2png,psfile+'.eps',/eps
  push,residfigs,psfile

  ;stop

endfor

; Convert png to pdf and combine
print,'Converting residual plots to PDF'
cd,current=curdir
cd,plotsdir
for i=0,n_elements(residfigs)-1 do begin
  residfigs1 = file_basename(residfigs[i])
  spawn,['convert',residfigs1+'.png',residfigs1+'.pdf'],/noshell
endfor
; Combine
print,'Writing combined residual plots to transphot_'+filter+'_resid_all'+pstag+'.pdf'
cmd = 'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=transphot_'+filter+'_resid_all'+pstag+'.pdf '
cmd += strjoin(file_basename(residfigs)+'.pdf',' ')
spawn,cmd,out,errout
cd,curdir

; --- Creating summary statistics for each standard star field ---
ui = uniq(arr.stdfield,sort(arr.stdfield))
ustdfield = arr[ui].stdfield
nstdfield = n_elements(ustdfield)
fieldstr = replicate({stdfield:'',stdnum:0L,ra:0.0d0,dec:0.0d0,sdss:0,nexp:0L,nmjd:0L,$
                      nsources:0L,nrejected:0L,medresid:0.0,sigresid:0.0},nstdfield)
; Get indices for each group
stdsi = sort(arr.stdfield)
stdbrk = where(shift(arr[stdsi].stdfield,1) ne arr[stdsi].stdfield,nstdbrk)
stdlo = stdbrk
stdhi = [stdbrk[1:nstdbrk-1]-1,n_elements(arr)-1]
for i=0,nstdfield-1 do begin
  ind = stdsi[stdlo[i]:stdhi[i]]
  nind = stdhi[i]-stdlo[i]+1
  if strmid(strupcase(ustdfield[i]),0,4) eq 'SDSS' or strmid(strupcase(ustdfield[i]),0,2) eq 'SA' then fieldstr[i].sdss=1
  uiexp = uniq(arr[ind].expnum,sort(arr[ind].expnum))
  uimjd = uniq(arr[ind].mjd,sort(arr[ind].mjd))
  fieldstr[i].stdfield = ustdfield[i]
  fieldstr[i].stdnum = i+1
  fieldstr[i].ra = median([arr[ind].ra])
  fieldstr[i].dec = median([arr[ind].dec])
  fieldstr[i].nexp = n_elements(uiexp)
  fieldstr[i].nmjd = n_elements(uimjd)
  fieldstr[i].nsources = nind
  fieldstr[i].nrejected = total(rejected_fixcolzpam[ind])
  fieldstr[i].medresid = median(resid_fixcolzpam[ind])
  fieldstr[i].sigresid = mad(resid_fixcolzpam[ind])
endfor


; --- Creating summary statistics for each exposure ---
ui = uniq(arr.expnum,sort(arr.expnum))
uexp = arr[ui].expnum
nexp = n_elements(uexp)
expstr = replicate({expnum:'',mjd:0L,nightnum:0L,exptime:0.0,airmass:0.0,seeing:0.0,ra:0.0d0,dec:0.0d0,$
                    stdfield:'',stdnum:0L,sdss:0,nsources:0L,nrejected:0L,medresid:0.0,sigresid:0.0},nexp)
; Get indices for each group
expsi = sort(arr.expnum)
expbrk = where(shift(arr[expsi].expnum,1) ne arr[expsi].expnum,nexpbrk)
explo = expbrk
exphi = [expbrk[1:nexpbrk-1]-1,n_elements(arr)-1]
for i=0,nexp-1 do begin
  ind = expsi[explo[i]:exphi[i]]
  nind = exphi[i]-explo[i]+1
  expstr[i].expnum = uexp[i]
  expstr[i].mjd = arr[ind[0]].mjd
  MATCH,mstr0.mjd,expstr[i].mjd,ind1,ind2,/sort
  expstr[i].nightnum = mstr0[ind1[0]].nightnum
  expstr[i].exptime = arr[ind[0]].exptime
  expstr[i].airmass = median([arr[ind].airmass])
  expstr[i].seeing = median([arr[ind].seeing])
  expstr[i].ra = median([arr[ind].ra])
  expstr[i].dec = median([arr[ind].dec])
  expstr[i].stdfield = arr[ind[0]].stdfield
  if strmid(strupcase(expstr[i].stdfield),0,4) eq 'SDSS' or strmid(strupcase(expstr[i].stdfield),0,2) eq 'SA'  then expstr[i].sdss=1
  MATCH,fieldstr.stdfield,arr[ind[0]].stdfield,ind1,ind2,/sort
  expstr[i].stdnum = fieldstr[ind1[0]].stdnum
  expstr[i].nsources = nind
  expstr[i].nrejected = total(rejected_fixcolzpam[ind])
  expstr[i].medresid = median(resid_fixcolzpam[ind])
  expstr[i].sigresid = mad(resid_fixcolzpam[ind])
endfor

;stop


; 9.) Making final output transformation equation structures
;------------------------------------------------------------
print,'' & print,'9.) Making final transformation equation structures'

; Information for each chip+night combination
; most of the info is in FITSTR_FIXCOLZPAM
fitstr = replicate({mjd:-1L,chip:-1L,nightnum:0,filter:'',color:'',colband:'',colsign:0,nstars:0L,nrejected:0L,seeing:99.9,zpterm:999.0,zptermerr:999.0,$
                     colterm:999.0,coltermerr:999.0,amterm:0.0,amtermerr:999.0,colamterm:0.0,colamtermerr:999.0,$
                     rms:999.0,sig:999.0,chisq:999.0,medresid:999.0,nbrt:0L,brtrms:999.0,brtsig:999.0,brtchisq:999.0,$
                     badsoln:-1,photometric:-1,amavgflag:0,namavg:0},n_elements(fitstr_fixcolzpam))
struct_assign,fitstr_fixcolzpam,fitstr
fitstr.filter = filter
fitstr.colamtermerr = 0.0
; Copy amtermerr over
for i=0,nmjds-1 do begin
  MATCH,fitstr.mjd,ntstr_fixcolzp[i].mjd,ind1,ind2,/sort,count=nmatch
  fitstr[ind1].amtermerr = ntstr_fixcolzp[ind2[0]].amtermsig
  fitstr[ind1].amavgflag = ntstr_fixcolzp[ind2[0]].amavgflag
  fitstr[ind1].namavg = ntstr_fixcolzp[ind2[0]].namavg
endfor
fitstr.color = colorname
dum = strsplit(colorname,'-',/extract)
if dum[0] eq filter then begin  ; filter - colband
  colband = dum[1]
  colsign = 1
endif else begin                ; colband - filter
  colband = dum[0]
  colsign = -1
endelse
fitstr.colband = colband
fitstr.colsign = colsign

; Chip-level color and zero-point terms
;   CHIPSTR has the fixed color terms
;   CHIPSTR_FIXCOLR has the relative zeropoint terms
fchipstr = chipstr_fixcolr
fchipstr.colterm = chipstr.colterm
fchipstr.coltermsig = chipstr.coltermsig
fchipstr.coltermerr = chipstr.coltermerr
add_tag,fchipstr,'filter',filter,fchipstr
add_tag,fchipstr,'color',colorname,fchipstr
add_tag,fchipstr,'colband',colband,fchipstr
add_tag,fchipstr,'colsign',colsign,fchipstr

; Nightly zeropoint and airmass terms
;  NTSTR_FIXCOLZP has the fixed airmass terms
;  MSTR_FIXCOLZPAM has the nightly zero-point terms
fntstr = replicate({mjd:0L,nightnum:0,date:'',nchips:0L,filter:'',nstars:0L,seeing:99.9,rms:99.0,brtrms:99.0,airmass0:0.0,airmass1:0.0,amrange:0.0,$
                    zpterm:99.0,zptermsig:99.0,amterm:99.0,amtermsig:99.0,colamterm:99.0,colamtermsig:99.0,$
                    badsoln:-1,photometric:-1,amavgflag:0,namavg:0},nmjds)
struct_assign,mstr_fixcolzpam,fntstr
fntstr.filter = filter
; Copy the AMTERM info from NTSTR_FIXCOLZP
fntstr.amterm = ntstr_fixcolzp.amterm
fntstr.amtermsig = ntstr_fixcolzp.amtermsig
fntstr.amavgflag = ntstr_fixcolzp.amavgflag
fntstr.namavg = ntstr_fixcolzp.namavg
add_tag,fntstr,'color',colorname,fntstr
add_tag,fntstr,'colband',colband,fntstr
add_tag,fntstr,'colsign',colsign,fntstr

; NEED TO ADD IN INFO FOR NIGHTS WITH NO STANDARD STAR DATA!!!
;  or do we?  how should we handle these??
; 20140105, could take nightly zero-point and amterm from neighboring
; nights, average them.


; DO I NEED A LINE FOR EVERY NIGHT+CHIP???

; Deal with nights that don't have enough stars for their own solutions


; -- FINAL OUTPUT --
; 1st extension: chip-level structure with color term (+err) and
;     zpterm (+err)
; 2nd extension: night-level structure with zpter (+err) and
;     amterm (+err) and photometric flag
outfile = 'transphot_'+filter+'_eqns.fits'
print,'Writing final photometric transformation equations to ',outfile
print,'HDU1: chip+night information'
print,'HDU2: unique chip-level information'
print,'HDU3: unique night-level information'
MKHDR,head,0
sxaddhist,'SMASH photometric transformation equuations for filter='+filter,head
sxaddhist,'HDU1: chip+night information (colterm, zeropoint)',head
sxaddhist,'HDU2: unique chip-level information (colterm, relative zeropoint)',head
sxaddhist,'HDU3: unique night-level information (amterm, nightly zeropoint)',head
sxaddpar,head,'filter',filter
sxaddhist,'NOREJECT = '+strtrim(keyword_set(noreject),2),head
sxaddhist,'NOSMITH = '+strtrim(keyword_set(nosmith),2),head
sxaddhist,'FITCHIPZPMJDTERM = '+strtrim(keyword_set(fitchipzpmjdterm),2),head
sxaddhist,'FIXCHIPZPTERM = '+strtrim(keyword_set(fixchipzpterm),2),head
FITS_WRITE,outfile,0,head
MWRFITS,fitstr,outfile,/silent
MWRFITS,fchipstr,outfile,/silent
MWRFITS,fntstr,outfile,/silent

save,arr,fitstr0,fitstr_fixcolr,fitstr_fixcolzp,fitstr_fixcolzpam,$
     mstr_fixcolzpam,resid_fixcolzpam,rejected_fixcolzpam,$
     expstr,fieldstr,file='solve_transphot_'+filter+'_resid.dat'

if keyword_set(stp) then stop

end
