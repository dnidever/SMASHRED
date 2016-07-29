pro solve_transphot_night,arr,mfitstr,fitcolam=fitcolam,fixcolr=fixcolr,fixam=fixam,fixcolam=fixcolam,$
                        fixchipzp=fixchipzp,resid=resid,expstr=expstr,verbose=verbose,silent=silent,pl=pl,$
                        save=save,bootstrap=bootstrap,stp=stp,rejected=rejected

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
mfitstr = replicate({mjd:-1L,chip:-1L,nightnum:0,nstars:-1L,seeing:99.9,zpterm:999.0,zptermerr:999.0,$
                     colterm:999.0,coltermerr:999.0,amterm:0.0,amtermerr:999.0,colamterm:0.0,colamtermerr:999.0,$
                     rms:999.0,sig:999.0,chisq:999.0,nbrt:0L,brtrms:999.0,brtsig:999.0,brtchisq:999.0,status:-1,$
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
        yerror = sqrt(mean((ydiff-median(ydiff))^2))  ; subtract out nightly offset
        status = 1
      endif

      chisq = mean((ydiff/err)^2)
      if n_elements(fixchipzp) gt 0 then chisq = mean(((ydiff-median(ydiff))/err)^2)  ; subtract nightly zpterm

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

      ; Stick information into chip structure
      mfitstr[i].zpterm = zpterm
      mfitstr[i].zptermerr = zptermerr
      mfitstr[i].colterm = colterm
      mfitstr[i].coltermerr = coltermerr
      mfitstr[i].rms = yerror
      mfitstr[i].sig = mad(ydiff)
      mfitstr[i].chisq = chisq
      mfitstr[i].status = status

      ; RMS of bright stars
      gbright = where(err lt 0.05 and mag lt min(mag)+3,nbright)
      mfitstr[i].nbrt = nbright
      if nbright gt 1 then begin
        brtrms = sqrt(mean((ydiff[gbright]-median(ydiff))^2))  ; subtract out nightly offset
        brtchisq = mean(((ydiff[gbright]-median(ydiff))/err[gbright])^2)
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
    coef = median(y[gd])
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
  if count gt 0 then begin
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
                        bootstrap=bootstrap,rejected=rejected


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
mstr = replicate({mjd:0L,nightnum:0,date:'',nchips:0L,nstars:0L,seeing:99.9,rms:99.0,brtrms:99.0,airmass0:0.0,airmass1:0.0,amrange:0.0,$
                  zpterm:99.0,zptermsig:99.0,colterm:99.0,coltermsig:99.0,amterm:99.0,amtermsig:99.0,colamterm:99.0,colamtermsig:99.0,$
                  badsoln:-1,photometric:-1},nmjds)
observatory,'ctio',obs
undefine,fitstr
if not keyword_set(silent) then $
  print,'  NUM    MJD     DATE     NSTARS    ZPTERM    COLTERM     AMTERM    AMRANGE   COLAMTERM     RMS       BRTRMS     MEDSIG'
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
                      verbose=verbose,pl=pl,save=save,bootstrap=bootstrap,rejected=rejected1
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
  totresid = mfitstr.nstars*mfitstr.rms^2
  rms = sqrt(mean(total(totresid)/nmatch))  ; rms for all stars of this night
  totbrtresid = mfitstr.nbrt*mfitstr.brtrms^2
  brtrms = sqrt(mean(total(totbrtresid)/total(mfitstr.nbrt)))
  medsig = median(mfitstr.sig)

  ; Save information to MJD structure
  mstr[i].mjd = umjds[i]
  ui = uniq(mfitstr.chip,sort(mfitstr.chip))
  mstr[i].nchips = n_elements(ui)
  mstr[i].nightnum = i+1
  mstr[i].nstars = nmatch
  mstr[i].rms = rms
  mstr[i].brtrms = brtrms
  mstr[i].zpterm = median(mfitstr.zpterm)
  mstr[i].zptermsig = mad(mfitstr.zpterm)
  mstr[i].colterm = median(mfitstr.colterm)
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
  print,i+1,umjds[i],date,nmatch,median(mfitstr.zpterm),median(mfitstr.colterm),$
        mfitstr[0].amterm,mstr[i].amrange,mfitstr[0].colamterm,rms,brtrms,medsig,format='(I5,I8,A10,I8,8F11.4)'
  ; Save structure
  push,fitstr,mfitstr

  ;stop

endfor

;stop

end

;----------------------

pro solve_transphot,file,fitcolam=fitcolam,errlim=errlim,arr=arr,mstr_fixcolzp=mstr_fixcolzp,resid_fixcolzp=resid_fixcolzp

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
; 4.) Fixing colterm. Redetermining zpterm and amter
;  -set photometric tag
;  -flag bad solutions with not enough stars
; 5.) Redetermining median zpterm per chip
; 6.) Fixing chip-level relative zpterm and colterm.  Retermining
;       nightly zpterm

; chip-level color terms
; chip-level zpterm
; nightly zpterm
; nightly airmass term

; -fit everything (colterm, zpterm, amterm)
; -determine chip-level color terms across everything
; -fix chip-level color terms, redetermine everything else
; -determine chip-level zpterm  + AMTERM???
; -fix chip-level color and zpterm, redetermine nightly zpterm and
;    airmass term
; -determine nightly airmass terms
; -fix color-level color and zpterm and nightly airmass term, redetermine
;    nightly zpterm 

reduxdir = '/data/smash/cp/red/photred/'

; Load the data
arr = MRDFITS(file,1)
narr = n_elements(arr)

tags = tag_names(arr)
filter = strlowcase(tags[6])  ; this should be the filter name

; Trim stars
if n_elements(errlim) eq 0 then errlim = 0.05
;if filter eq 'u' then errlim=0.08
print,'Imposing error cut, ERR<=',strtrim(errlim,2),' mag'
gd = where(arr.err lt errlim,ngd)
arr = arr[gd]
narr = n_elements(arr)

; u-band, removing very blue and very red stars
if filter eq 'u' then begin
  print,'Removing very blue/red stars'
  gd = where(arr.(8) ge 0.73 and arr.u_g lt 2.5,ngd)  ; u-g
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
for i=0,nmjd-1 do begin
  ind = where(expstr.mjd eq umjd[i],nind)
  expstr[ind].residzeroptnsig = abs(expstr[ind].residzeropt) / mjdzpsig[i]
endfor
; Find bad exposures
bdexp = where(expstr.residzeroptnsig gt 5,nbdexp)
if nbdexp gt 0 then print,'Removing ',strtrim(nbdexp,2),' bad exposures: ',expstr[bdexp].expnum
undefine,bdind
for i=0,nbdexp-1 do begin
  ind = expsi[explo[bdexp[i]]:exphi[bdexp[i]]]
  push,bdind,ind
endfor
REMOVE,bdind,arr
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
  if nmatch gt 0 then begin
    remove,ind1,arr
    print,'Removing bad exposure: ',badexp[i]
  endif
endfor


; 2.) Run on all the nights fitting everything individually
;------------------------------------------------------------
print,'' & print,'2.) Run on all the nights fitting everything individually'
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr0,fitstr0,resid=resid0,expstr=expstr0,/bootstrap  ;,/save
nmjds = n_elements(mstr0)
uimjd = uniq(mstr0.mjd,sort(mstr0.mjd))
umjd = mstr0[uimjd]
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
chipstr = replicate({chip:0L,zpterm:99.0,zptermsig:99.0,zptermerr:99.0,colterm:99.0,coltermsig:99.0,coltermerr:99.0,rms:99.0,medrms:99.0},nuchips)
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
    chipstr[i].zpterm = median(fitstr0[ind[gdphot]].relzpterm)
    chipstr[i].zptermsig = mad(fitstr0[ind[gdphot]].relzpterm)
    chipstr[i].zptermerr = chipstr[i].zptermsig / sqrt(ngdphot)
    chipstr[i].colterm = median([fitstr0[ind[gdphot]].colterm])
    chipstr[i].coltermsig = mad([fitstr0[ind[gdphot]].colterm])
    chipstr[i].coltermerr = chipstr[i].coltermsig / sqrt(ngdphot)
    totresid = fitstr0[ind[gdphot]].nbrt*fitstr0[ind[gdphot]].brtrms^2
    rms = sqrt(mean(total(totresid)/total(fitstr0[ind[gdphot]].nbrt)))  ; brt rms
    chipstr[i].rms = rms
    chipstr[i].medrms = median([fitstr0[ind[gdphot]].rms])
    fitstr0[ind].chiprelzpterm = fitstr0[ind].relzpterm - chipstr[i].zpterm   ; relative MJD and chip medians
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
file = 'transphot_'+filter+'_zpterm_chip'
ps_open,file,/color,thick=4,/encap
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
ps2png,file+'.eps',/eps
; COLTERM vs. chip
file = 'transphot_'+filter+'_colterm_chip'
ps_open,file,/color,thick=4,/encap
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
ps2png,file+'.eps',/eps
; CHIPRELZPTERM vs. NIGHTNUM
file = 'transphot_'+filter+'_chiprelzpterm_nightnum'
ps_open,file,/color,thick=4,/encap
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
ps2png,file+'.eps',/eps
; CHIPRELCOLTERM vs. NIGHTNUM
file = 'transphot_'+filter+'_chiprelcolterm_nightnum'
ps_open,file,/color,thick=4,/encap
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
ps2png,file+'.eps',/eps



; 4.) Fixing colterm, Redetermine zpterm and airmass term
;---------------------------------------------------------
print,'' & print,'4.) Fixing colterm.  Redetermining zpterm and amterm'
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolr,fitstr_fixcolr,fixcolr=chipstr,/bootstrap

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
chipstr_fixcolr = replicate({chip:0L,zpterm:99.0,zptermsig:99.0,zptermerr:99.0,colterm:99.0,coltermsig:99.0,coltermerr:99.0,rms:99.0,medrms:99.0},nuchips)
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
    chipstr_fixcolr[i].zpterm = median(fitstr_fixcolr[ind[gdphot]].relzpterm)
    chipstr_fixcolr[i].zptermsig = mad(fitstr_fixcolr[ind[gdphot]].relzpterm)
    chipstr_fixcolr[i].zptermerr = chipstr_fixcolr[i].zptermsig / sqrt(ngdphot)
    chipstr_fixcolr[i].colterm = median([fitstr_fixcolr[ind[gdphot]].colterm])
    chipstr_fixcolr[i].coltermsig = mad([fitstr_fixcolr[ind[gdphot]].colterm])
    chipstr_fixcolr[i].coltermerr = chipstr_fixcolr[i].coltermsig / sqrt(ngdphot)
    totresid = fitstr_fixcolr[ind[gdphot]].nbrt*fitstr_fixcolr[ind[gdphot]].brtrms^2
    rms = sqrt(mean(total(totresid)/total(fitstr_fixcolr[ind[gdphot]].nbrt)))  ; brt rms
    chipstr_fixcolr[i].rms = rms
    chipstr_fixcolr[i].medrms = median([fitstr_fixcolr[ind[gdphot]].rms])
    fitstr_fixcolr[ind].chiprelzpterm = fitstr_fixcolr[ind].relzpterm - chipstr_fixcolr[i].zpterm   ; relative MJD and chip medians
    fitstr_fixcolr[ind].chiprelcolterm = fitstr_fixcolr[ind].colterm - chipstr_fixcolr[i].colterm
    print,chipstr_fixcolr[i].chip,chipstr_fixcolr[i].zpterm,chipstr_fixcolr[i].zptermerr,chipstr_fixcolr[i].colterm,$
          chipstr_fixcolr[i].coltermerr,chipstr_fixcolr[i].rms,format='(I5,5F10.5)'
  endif
endfor
; COLTERMERR = 0.0 because it was fixed




; 6.) Redetermine NIGHTLY zpterm and amterm with chip-level relative zpterm 
;      and colterm fixed
;---------------------------------------------------------------------------
print,'' & print,'6.) Fixing chip-level relative zpterm and colterm.  Redetermining NIGHTLY zpterm and amterm'
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolzp,fitstr_fixcolzp,fixcolr=chipstr,fixchipzp=chipstr_fixcolr,resid=resid_fixcolzp,/bootstrap

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
print,'' & print,'7.) Averagine airmass terms'
ntstr_fixcolzp = mstr_fixcolzp
add_tag,ntstr_fixcolzp,'amavgflag',0,ntstr_fixcolzp  ; amterm averaged
add_tag,ntstr_fixcolzp,'namavg',0,ntstr_fixcolzp  ; amterm averaged
print,'  NUM   MJD    AMTERM   AMTERMSIG FLAG  NAVG   COMMENT'
for i=0,nmjd-1 do begin
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
        gdind = where(abs(resid-poly(arr[useind].airmass,coef1)) lt 4.0*sig1,ngdind)
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
file = 'transphot_'+filter+'_amterm'
ps_open,file,/color,thick=4,/encap
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
ps2png,file+'.eps',/eps


; 7.) Redetermine NIGHTLY zpterm holding everything else fixed
;--------------------------------------------------------------
print,'' & print,'7.) Redetermine NIGHTLY zpterm holding everything else fixed'
SOLVE_TRANSPHOT_ALLNIGHTS,arr,mstr_fixcolzpam,fitstr_fixcolzpam,fixcolr=chipstr,fixchipzp=chipstr_fixcolr,$
                        fixam=ntstr_fixcolzp,resid=resid_fixcolzpam,/bootstrap

; Reset PHOTOMETRIC tag
MATCH,mstr_fixcolzpam.mjd,conditions.mjd,ind1,ind2
mstr_fixcolzpam[ind1].photometric = conditions[ind2].photometric
for i=0,n_elements(conditions)-1 do begin
  MATCH,fitstr_fixcolzpam.mjd,conditions[i].mjd,ind1,ind2,count=nmatch,/sort
  if nmatch gt 0 then fitstr_fixcolzpam[ind1].photometric=conditions[i].photometric
endfor

; Flag "bad" solution with not Not enough stars
bdmjd = where(mstr_fixcolzpam.nstars lt 100,nbdmjd,comp=gdmjd,ncomp=ngdmjd)
if nbdmjd gt 0 then mstr_fixcolzpam[bdmjd].badsoln = 1    ; not enough stars
if ngdmjd gt 0 then mstr_fixcolzpam[gdmjd].badsoln = 0    ; enough stars
fitstr_fixcolzpam.badsoln = 0
for i=0,nbdmjd-1 do begin
  MATCH,fitstr_fixcolzpam.mjd,mstr_fixcolzpam[bdmjd[i]].mjd,ind1,ind2,count=nmatch,/sort
  fitstr_fixcolzpam[ind1].badsoln = 1
endfor

; Make final residual figures
print,'Making final residual plots'
for i=0,nmjds-1 do begin

  MATCH,arr.mjd,mstr_fixcolzpam[i].mjd,ind1,ind2,count=nmatch,/sort

  ; Make resid plots, one per night
  ;  resid vs. color, resid vs. airmass, resid vs. chip number

  file = 'transphot_'+filter+'_'+strtrim(mstr_fixcolzpam[i].mjd,2)+'_resid'
  ps_open,file,/color,thick=4,/encap
  device,/inches,xsize=9.5,ysize=9.5
  x0 = 0.08 & x1=0.97
  y0 = 0.05 & dy = 0.30
  sig = mad(resid_fixcolzpam[ind1])
  ;yr = minmax(resid_fixcolzpam[ind1])
  yr = [-1,1]*(3.5*sig > 0.2)
  yr = [yr[0] > (-1.0), yr[1] < 1.0]
  ; Resid vs. color
  xr = minmax(arr[ind1].(8))
  xr = [-0.5,2.0 < (xr[1]+0.05) > 3.5]
  plotc,arr[ind1].(8),resid_fixcolzpam[ind1],arr[ind1].err,ps=3,xr=xr,yr=yr,xs=1,ys=1,xtit='Color',ytit='Residuals',$
        pos=[x0,y0,x1,dy],colpos=[x0,0.97,x1,0.99],xticklen=0.04
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'<COLTERM> = '+stringize(mstr_fixcolzpam[i].colterm,ndec=4)+'+/-'+$
         stringize(mstr_fixcolzpam[i].coltermsig,ndec=4),align=0,charsize=1.3
  oplot,xr,[0,0],linestyle=2
  ; Resid vs. airmass
  xr = minmax(arr[ind1].airmass)
  xr = [0.9,(xr[1]+0.05) > 2.0]
  plotc,arr[ind1].airmass,resid_fixcolzpam[ind1],arr[ind1].err,ps=3,xr=xr,yr=yr,xs=1,ys=1,xtit='Airmass',ytit='Residuals',/noerase,$
        pos=[x0,y0+dy,x1,2*dy],/nocolorbar,xticklen=0.04
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'AMTERM = '+stringize(mstr_fixcolzpam[i].amterm,ndec=4)+'+/'+$
         stringize(ntstr_fixcolzp[i].amtermsig,ndec=4),align=0,charsize=1.3
  if mstr_fixcolzpam[i].photometric eq 1 then photcom='PHOTOMETRIC' else photcom='NON-PHOTOMETRIC'
  xyouts,xr[1]-0.03*range(xr),yr[1]-0.12*range(yr),photcom,align=1,charsize=1.3
  oplot,xr,[0,0],linestyle=2
  ; Resid vs. chip number
  rnd = randomu(seed,nmatch)*0.6-0.3
  xr = [0,63]
  plotc,arr[ind1].chip+rnd,resid_fixcolzpam[ind1],arr[ind1].err,ps=3,xr=xr,yr=yr,xs=1,ys=1,xtit='Chip',ytit='Residuals',/noerase,xticklen=0.04,$
        tit=strtrim(mstr_fixcolzpam[i].mjd,2)+' - Residuals vs. color/airmass/chip (color-coded by ERR)',pos=[x0,y0+2*dy,x1,3*dy],/nocolorbar
  xyouts,xr[0]+0.03*range(xr),yr[1]-0.12*range(yr),'ZPTERM = '+stringize(mstr_fixcolzpam[i].zpterm,ndec=4)+'+/-'+$
         stringize(mstr_fixcolzpam[i].zptermsig,ndec=4),align=0,charsize=1.3
  xyouts,xr[1]-0.03*range(xr),yr[1]-0.12*range(yr),'RMS = '+stringize(mstr_fixcolzpam[i].rms,ndec=4),align=1,charsize=1.3
  oplot,xr,[0,0],linestyle=2

  ps_close
  ps2png,file+'.eps',/eps

endfor

stop


; 8.) Making final chip and night-level structures
;-------------------------------------------------
print,'' & print,'8.) Making final chip and night-level structures'


; Chip-level color and zero-point terms

; Nightly zeropoint and airmass terms

; DO I NEED A LINE FOR EVERY NIGHT+CHIP???





; Deal with nights that don't have enough stars for their own solutions
;----------------------------------------------------------------------



; Deal with non-photometric nights
;---------------------------------
; use median chip-level color terms
; set zpterm, amterm and colamter to ZERO


stop



for i=0,nbdmjd-1 do begin

  ; Fit global zpterm, use colterm and amterm from neighbors

  stop

endfor

; FINAL OUTPUT
; 1st extension: chip-level structure with color term (+err) and
;     zpterm (+err)
; 2nd extension: night-level structure with zpter (+err) and
;     amterm (+err) and photometric flag


stop

rnd = randomu(seed,nfitstr)*0.6-0.3
plotc,fitstr.nightnum+rnd,fitstr.zpterm,fitstr.sig,ps=3,xr=[0,nmjds+1],yr=[-1,1],xs=1,ys=1,max=0.1,$
      xtit='Night Number',ytit='ZP term',tit='Zero-point term (color-coded by zptermerr)'


; Photometric night, based on zpterm scatter
bdmjd = where(mstr.brtrms gt median(mstr.brtrms)+2.5*mad(mstr.brtrms),nbdmjd)
print,strtrim(nbdmjd,2),' nights have high RMS'
if nbdmjd gt 0 then print,umjds[bdmjd]


stop

; Unique CHIPS
ui = uniq(fitstr.chip,sort(fitstr.chip))
uchips = fitstr[ui].chip
nuchips = n_elements(uchips)

; Measure median and RMS per chip
chipmedzpterm = fltarr(nuchips)+99.99
chipzptermsig = fltarr(nuchips)+99.99
chipzptermerr = fltarr(nuchips)+99.99
chiprelzpterm = fltarr(nfitstr)+99.99
chipmedcolterm = fltarr(nuchips)+99.99
chipcoltermsig = fltarr(nuchips)+99.99
chipcoltermerr = fltarr(nuchips)+99.99
chiprelcolterm = fltarr(nfitstr)+99.99
chipmedrms = fltarr(nuchips)+99.99

for i=0,nuchips-1 do begin
  ind = where(fitstr.chip eq uchips[i],nind)
  if nind gt 0 then begin
    chipmedzpterm[i] = median([fitstr[ind].zpterm])
    chipzptermsig[i] = mad([fitstr[ind].zpterm])
    chipzptermerr[i] = chipzptermsig[i] / sqrt(nind)
    chiprelzpterm[ind] = fitstr[ind].zpterm - chipmedzpterm[i]
    chipmedcolterm[i] = median([fitstr[ind].colterm])
    chipcoltermsig[i] = mad([fitstr[ind].colterm])
    chipcoltermerr[i] = chipcoltermsig[i] / sqrt(nind)
    chiprelcolterm[ind] = fitstr[ind].colterm - chipmedcolterm[i]
    chipmedrms[i] = median([fitstr[ind].rms])
  endif
endfor

plotc,fitstr.chip+rnd,fitstr.zpterm,fitstr.mjd,ps=3,yr=[-0.5,0.5],xs=1
plotc,fitstr.chip+rnd,fitstr.zpterm,fitstr.zptermerr,ps=3,yr=[-0.5,0.5],xs=1,max=0.01
oplot,uchips,chipmedzpterm,ps=-1
plotc,fitstr.chip+rnd,fitstr.colterm,fitstr.mjd,ps=3,yr=[-0.3,0.1],xs=1
plotc,fitstr.chip+rnd,fitstr.colterm,fitstr.coltermerr,ps=3,yr=[-0.3,0.1],xs=1,max=0.01
oplot,uchips,chipmedcolterm,ps=-1


; Measure median relative chip zpterm/colterm per night
mjdchipmedzpterm = fltarr(nmjds)
mjdchipzptermsig = fltarr(nmjds)
mjdchipzptermerr = fltarr(nmjds)
mjdchipmedcolterm = fltarr(nmjds)
mjdchipcoltermsig = fltarr(nmjds)
mjdchipcoltermerr = fltarr(nmjds)
nightnum = lonarr(nfitstr)
for i=0,nmjds-1 do begin
  ind = where(fitstr.mjd eq umjds[i],nind)
  ; ZP term
  mjdchipmedzpterm[i] = median( chiprelzpterm[ind] )
  mjdchipzptermsig[i] = mad(chiprelzpterm[ind])
  mjdchipzptermerr[i] = mjdchipzptermsig[i] / sqrt(nind)
  ; Color term
  mjdchipmedcolterm[i] = median( chiprelcolterm[ind] )
  mjdchipcoltermsig[i] = mad(chiprelcolterm[ind])
  mjdchipcoltermerr[i] = mjdchipcoltermsig[i] / sqrt(nind)
  nightnum[ind] = i+1   ; night number
endfor

; Zero-point
plotc,fitstr.mjd+rnd,chiprelzpterm,fitstr.coltermerr,ps=3,yr=[-0.1,0.1],max=0.01,xtit='MJD',ytit='Relative ZP Term',tit='ZP relative to median per chip'
oplot,umjds,mjdchipmedzpterm,ps=1,sym=2

; Color term vs. night
;plotc,fitstr.mjd+rnd,chiprelcolterm,fitstr.coltermerr,ps=3,yr=[-0.02,0.02],max=0.01,xtit='MJD',ytit='Relative Color Term',tit='Color term relative to median per chip'
;oplot,umjds,mjdchipmedcolterm,ps=1,sym=2
z = fitstr.zpterm ; fitstr.coltermerr
;zr = minmax(z)
zr = [-0.5,0.5]
ztit = 'zpterm'
plotc,nightnum+rnd,chiprelcolterm,z,ps=8,sym=0.6,xs=1,ys=1,xr=[0,max(nightnum)+1],yr=[-0.06,0.06],min=zr[0],max=zr[1],$
      xtit='Night Number',ytit='Relative Color Term',tit='Color term relative to median per chip (color-coded by '+ztit+')'
;col = scale(fitstr.coltermerr,[min(fitstr.coltermerr),0.01],[50,250])<250
col = 50 > scale(z,zr,[50,250]) < 250
for i=0,nfitstr-1 do oplot,[0,0]+nightnum[i]+rnd[i],fitstr[i].coltermerr*[-1,1]+chiprelcolterm[i],co=col[i]
oplot,indgen(nmjds)+1,mjdchipmedcolterm,ps=1,sym=2,thick=2
oplot,[0,100],[-1,-1]*mad(mjdchipmedcolterm),linestyle=2
oplot,[0,100],[1,1]*mad(mjdchipmedcolterm),linestyle=2

; These two nights have lower color-terms on average by 0.01
;  43   57362  20151205   43397 
;  44   57363  20151206   37811 
; This one is a "normal" color-term night
;   42   57336  20151109   36156 
; if you plot resid vs. color for this night it is very tight,
; the same plot for 57362 shows much more variation due to airmass
; could it be a color*extinction term that we need????

; fit color terms for 57362, there is a small difference with airmass
; at the 0.01 level
; coef for airmass<1.2
;    0.0556403     0.113632  (intercept, slope)
; coef for airmass>1.6
;   -0.0722108     0.121178

; same for 57336, smaller effect
;    0.0646907     0.112354
; -0.000309065     0.118200

;    38   57097  20150315   14386    -2.4869    -0.1065, very low
;         zpterm (??), large scatter in color terms

; There definitively are mean color term variations for a few
; nights, it's about a ~3% variation
; COULD THESE BE NON-PHOTOMETRIC NIGHTS??


; Color term vs. chip
z = nightnum
zr = minmax(z)
ztit = 'nightnum'
plotc,fitstr.chip+rnd,chiprelcolterm,z,ps=8,sym=0.6,xs=1,ys=1,xr=[0,max(uchips)+1],yr=[-0.06,0.06],min=zr[0],max=zr[1],$
      xtit='Chip Number',ytit='Relative Color Term',tit='Color term relative to median per chip (color-coded by '+ztit+')'
col = 50 > scale(z,zr,[50,250]) < 250
for i=0,nfitstr-1 do oplot,[0,0]+fitstr[i].chip+rnd[i],fitstr[i].coltermerr*[-1,1]+chiprelcolterm[i],co=col[i]
oplot,indgen(nmjds)+1,chipmedcolterm,ps=1,sym=2,thick=2
oplot,[0,100],[0,0],linestyle=1
oplot,[0,100],[-1,-1]*mad(chiprelcolterm),linestyle=2
oplot,[0,100],[1,1]*mad(chiprelcolterm),linestyle=2
; I don't see any real systematics here


; Discrepant nights
bd = where(abs(mjdchipmedcolterm) gt 0.04,nbd)
print,strtrim(nbd,2),' nights with discrepant median color terms'
if nbd gt 0 then print,umjds[bd]

stop

; some chip-to-chip variations
; some night-to-night variations, especially between two years
;   could some of those be from using SDSS and the Smith standards?
;  maybe it was a non-photometric night???

; These groups of days have lower color terms then the ones after that
; 56369-56372 and 56509+56510
; March 2013 and Aug 2013
; you see this in g-band and r-band but not i-band or u-band
; z-band doesn't show a constant offset, but the chip-to-chip patterns
; seem to change

; 56677-56686 color-terms are a bit high, but these are the Jan 21-28, 2014
; days when only small-field Smith/DES fields were observed and so
; only a handful of chips have observations
; the one day in that range that has low RMS seems to have a similar
; color-term in the r-band.

plotc,mchfitstr.chip+rnd,mchfitstr.rms,mchfitstr.mjd,ps=3,yr=[-2,2],xs=1  


; airmass terms
ploterror,mfitstr.mjd-min(mfitstr.mjd),mfitstr.coef[0],mfitstr.mjd*0,mfitstr.coeferr[1],ps=1,yr=[-1,1]

; 56678 and 56680 have less than 50 stars
; 56677, 56678, 56680, 56683-56686, less than 600 stars, not good fits
; maybe use the "standard" value for these


rndmjd = randomu(seed,numjdchip)*0.8-0.4
plotc,mchfitstr.mjd+rndmjd,mchfitstr.coef[0],mchfitstr.chip,ps=3,yr=[-1,1],xs=1
plotc,mchfitstr.mjd+rndmjd,mchfitstr.coef[1],mchfitstr.chip,ps=3,yr=[-0.3,0.1],xs=1
; plot using unique night numbers



; Measure median and RMS per chip
chipmed = fltarr(2,nuchips)+99.99
chiprms = fltarr(nuchips)+99.99
for i=0,nuchips-1 do begin
  ind = where(mchfitstr2.chip eq uchips[i],nind)
  if nind gt 0 then begin
    chipmed[0,i] = median([mchfitstr2[ind].coef[0]])
    chipmed[1,i] = median([mchfitstr2[ind].coef[1]])
    chiprms[i] = median([mchfitstr2[ind].rms])
  endif
endfor


plothist,mfitstr.medrms,bin=0.04,xr=[0,2]

bdmjd = where(mfitstr.rms gt median(mfitstr.rms)+2.5*mad(mfitstr.rms),nbdmjd)
print,strtrim(nbdmjd,2),' nights have high RMS'
if nbdmjd gt 0 then print,umjds[bdmjd]

; There are definitely some outlier stars, should remove them
plot,arr.mag,resid3,ps=3,yr=[-1,1]

; Compare the two RMS values, should improve slightly with airmass
; term removed
plot,mchfitstr.rms,ps=1,yr=[0,1]
oplot,mchfitstr2.rms,ps=4,co=250

stop

end
