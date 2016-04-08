;+
;
; STDRED_TRANSPHOT_SMASH
;
; This determined the photometric transformation equations from
; the standard star data via STDRED.  Modifications for SMASH.
;
; This program finds the photometric transformation equations
; with iterative outlier rejection given the proper input file
; (same as SKAWDPHOT.PRO).
;
; INPUTS:
;  input    Filename of input data (i.e. g.cat).
;  =fitzp   The type of fitting to use for the zerpoint term(s).
;             fitam=1  separate value for each night (default, unless /sepchip)
;             fitam=2  global nightly variations and chip-specific variations
;                           fixed across nights (default if /sepchip)
;  =fitam   The type of fitting to use for the airmass/extinction terms(s).
;             fitam=1  one global value for all nights (default)
;             fitam=2  separate extinction terms for each night
;  =fitcolr The type of fitting to use for the color terms(s).
;             fitam=1  one global value for all nights and chips (default)
;             fitam=2  separate color terms for each night
;  /fitac   Fit the airmass*color term.  The default is keep it fixed at 0.0.
;  /fitcolsq Fit the color*color term.  The default is keep it fixed at 0.0.
;  =fixam   Fix the airmass term to this value.
;  =fixcolr Fix the color term to this value.
;  =fixac   Fix the airmass*color term to this value.  The default is 0.0.
;  =fixcolsq Fix the color*color term to this value.  The default is 0.0.
;  =errlim  Only use stars with observed errors lower than this value.
;  =inparr  Input the structure instead of the filename
;  /plotresid  Plot the residuals.
;  =yrange     Yrange for residual plotting.
;  /nooutput  Don't print out any of the transformation files.
;  /sepchip  Each chip separately.
;  /silent   Don't print anything to the screen.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  The transformation equations are written to a file called
;  input+'.trans'.  Separate transformation equations for each
;  night are also put in the file.  Also, each night's transformation
;  equations are written to their own files called 'n#MAG.trans'
;  (i.e. n1M.trans').
;
; =arr      The data structure.
; =tlines   An array of the transformation equations for each night.
; =trans    The transformation structure, one element for each night.  
; =rms      The final rms of the fit.
;
; USAGE:
;  IDL>stdred_transphot_smash,'g.cat'
;
; By D. Nidever   Jan. 2008
;       April 2016, modifications
;-

pro std_trans_dummy
; This makes it so that you don't have to compile before running
FORWARD_FUNCTION std_devtrans, std_transfunc
end

;--------------------

function std_devtrans,p,y=y,mag=mag,col=col,am=am,night=night,chip=chip,$
                      mapnight=mapnight,weight=weight

model = std_transfunc(p,mag=mag,col=col,am=am,night=night,chip=chip,mapnight=mapnight)

return,(y-model)*weight

end

;---------------------------

function std_transfunc,p,mag=mag,col=col,am=am,night=night,chip=chip,mapnight=mapnight


;nnights = max(night)-min(night)+1
;zeropar = P[0:nnights-1]
;ampar = P[nnights]
;colpar = P[nnights+1]
;nightzero = zeropar[night-min(night)]

npar = n_elements(p)

;zeropar = P[0:npar-3]
;ampar = P[npar-2]
;colpar = P[npar-1]
zeropar = P[0:npar-5]
ampar = P[npar-4]
colpar = P[npar-3]
amcolpar = P[npar-2]
colsqpar = P[npar-1]

nightzero = zeropar[mapnight[night]]

;   V = mV - v1 - v2 * XV - v3 * (B-V) - v4 * XV*(B-V) - v5 * (B-V) * (B-V)
;         +(aperture correction) + (time correction)
;model = mag + nightzero + col*colpar + am*ampar
model = mag - nightzero - col*colpar - am*ampar - col*am*amcolpar - col*col*colsqpar

return,model

end 

;---------------------

pro stdred_transphot_smash,input,stp=stp,arr=arr,plotresid=plotresid,$
                  yrange=yrange,fitzp=fitzp,fitam=fitam,fitcolr=fitcolr,fitac=fitac,fitcolsq=fitcolsq,$
                  fixam=fixam,fixcolr=fixcolr,fixac=fixac,fixcolsq=fixcolsq,trans=trans,$
                  tlines=tlines,rms=rms,errlim=errlim,nooutput=nooutput,$
                  inparr=inparr,sepchip=sepchip,silent=silent

; Not enough inputs
if n_elements(input) eq 0 and n_elements(inparr) eq 0 then begin
  print,'Syntax - stdred_transphot,input,inparr=inparr,stp=stp,plotresid=plotresid,arr=arr'
  print,'                          fitzp=fitzp,fitam=fitam,fitcolr=fitcolr,fitac=fitac,fitcolsq=fitcolsq'
  print,'                          fixam=fixam,fixcolr=fixcolr,fixac=fixac,fixcolsq=fixcolsq'
  print,'                          errlim=errlim,sepchip=sepchip,nooutput=nooutout,silent=silent'
  return
endif

; Loading the data
if n_elements(input) gt 0 then if input ne '' then begin

  print,'Loading the data'
  ext = first_el(strsplit(input[0],'.',/extract),/last)
  if ext eq 'gz' then begin
    dum = strsplit(input[0],'.',/extract)
    ndum = n_elements(dum)
    ext = strjoin(dum[ndum-2:ndum-1],'.')  
  endif 
  ; Go through the possible cases
  case ext of
  'fits': arr=MRDFITS(input[0],1)
  'fits.gz': arr=MRDFITS(input[0],1)
  'cat': begin
    ; Get the fieldnames and fieldtypes
    ; We need ID to be STRING not LONG
    tempfile = MKTEMP('temp')
    SPAWN,'head '+input[0]+' >> '+tempfile,out,errout
    temp = IMPORTASCII(tempfile,/header,/noprint)
    FILE_DELETE,tempfile,/allow,/quiet
    fieldnames = TAG_NAMES(temp)
    nfieldnames = n_elements(fieldnames)
    fieldtypes = lonarr(nfieldnames)
    for k=0,nfieldnames-1 do fieldtypes[k] = SIZE(temp[0].(k),/type)
    idind = where(fieldnames eq 'ID',nidind)
    fieldtypes[idind[0]] = 7
    arr = IMPORTASCII(input[0],fieldnames=fieldnames,fieldtypes=fieldtypes,skip=1,/noprint)
    ;arr = importascii(input[0],/header,/noprint)
  end
  else: begin
    print,'Extension ',ext,' NOT SUPPORTED'
    return
  end
  endcase

; Using structure input
endif else begin
  arr = inparr
endelse
narr = n_elements(arr)
;orig = arr

; Defaults
if n_elements(errlim) eq 0 then errlim=0.03
if n_elements(fitzp) eq 0 and keyword_set(sepchip) then fitzp=2
if n_elements(fitzp) eq 0 then fitzp=1
if n_elements(fitam) eq 0 then fitam=1
if n_elements(fitcolr) eq 0 then fitcolr=1
if n_elements(fitac) eq 0 then fitac=0
if n_elements(fitcolsq) eq 0 then fitcolsq=0

if not keyword_set(silent) then begin
  case fitzp of
     1: print,'FITZP=1  Nightly variations'
     2: print,'FITZP=2  Global nightly variations and chip-specific variations fixed across nights'
     else: print,'FITZP=',strtrim(fitzp,2),'  ???'
  endcase
  case fitam of
     1: print,'FITAM=1  One airmass term for all nights'
     2: print,'FITAM=2  Separate airmass term for each night'
     else: print,'FITAM=',strtrim(fitam,2),'  ???'
  endcase
  case fitcolr of
     1: print,'FITCOLR=1  One global color term for all chips and nights'
     2: print,'FITCOLR=2  Separate color term for each chip'
     else: print,'FITCOLR=',strtrim(fitcolr,2),'  ???'
  endcase
end


; Each chip separately
;----------------------
;  This is a quick KLUDGE
if (keyword_set(sepchip) or fitzp eq 2 or fitam eq 2 or fitcolr eq 2) and tag_exist(arr,'CHIP') then begin
  print,'' & print,'Solving all CHIPS separately'

  ; What we are fitting:  zeropoint, airmass, color
  ; can vary: night, chip
  ;
  ; zeropoint:
  ;  1. nightly variations
  ;  2. global nightly variations and chip-specific variations
  ;      fixed across nights
  ;
  ; airmass:
  ;  1. one airmass term for all nights
  ;  2. separate airmass term for each night
  ;
  ; color:
  ;  1. one color term for all chips and nights
  ;  2. separate color term for each chip
  ;  3. separate color term for each chip/night, Probably don't want this!

  ; 1.) Initial solution for all data combined (zp=1,am=1,col=1)
  ; 2.) Initial separate solutions for each chip
  ; 3.) Fit nightly airmass terms
  ; 4.) Remove airmass terms and fit chip-specfic color-terms and chip/night
  ;      specific zero-point terms
  ; 5.) Measure nightly zero-point term
  ; 6.) Remove airmass terms and nightly zpterm and fit for
  ;       chip-specific zero-point and color terms
  ; 7.) Construct final transformation equations
  
  ; Nights/MJD
  ui = uniq(arr.mjd,sort(arr.mjd))
  unights = arr[ui].mjd
  nnights = n_elements(unights)
  ;ui = uniq(arr.night,sort(arr.night))
  ;unights = arr[ui].night
  ;nnights = n_elements(unights)

  ; Chips
  ui = uniq(arr.chip,sort(arr.chip))
  chips = arr[ui].chip
  nchips = n_elements(chips)

  ;; Map from chip number to index in alltrans
  ;mapchip = lonarr(max(chips)+1)-1
  ;mapchip[chips] = indgen(nchips)

  ; ---Check how many data points there are for each chip/night combination---
  print,'' & print,'The matrix of datapoints per night and chip'
  npts_chipnight = lonarr(nnights,nchips)
  chipnight = strtrim(arr.chip,2)+'-'+strtrim(arr.mjd,2)
  for i=0,nnights-1 do begin
    for j=0,nchips-1 do begin
      ;gd = where(arr.night eq unights[i] and arr.chip eq chips[j] and arr.err le errlim,ngd)
      ;gd = where(arr.mjd eq unights[i] and arr.chip eq chips[j] and arr.err le errlim,ngd)
      ichipnight = strtrim(chips[j],2)+'-'+strtrim(unights[i],2)
      MATCH,chipnight,ichipnight,ind1,ind2,count=ngd,/sort
      npts_chipnight[i,j] = ngd
    endfor
  endfor
  ; print to the screen
  print,'NIGHT   ',long(unights)
  for i=0,nchips-1 do print,'CHIP = ',strtrim(chips[i],2),reform(npts_chipnight[*,i])

;stop

  ; ---Check airmass coverage for each night
  print,'' & print,'Airmass coverage per night'
  print,'Night  RANGE    MIN   MAX   Individual airmasses per frame'
  for i=0,nnights-1 do begin
    ;gd = where(arr.night eq unights[i] and arr.err le errlim,ngd)
    MATCH,arr.mjd,unights[i],ind1,ind2,/sort
    gd1 = where(arr[ind1].err le errlim,ngd)
    gd = ind1[gd1]
    ;gd = where(arr.mjd eq unights[i] and arr.err le errlim,ngd)
    arr0 = arr[gd]
    ; get unique frames, trim chip number
    frame1 = reform( (strsplitter(arr0.frame,'_',/extract))[0,*] )
    ui = uniq(frame1,sort(frame1))  
    uframe1 = frame1[ui]
    nuframe1 = n_elements(uframe1)
    medam_uframe1 = fltarr(nuframe1)
    for j=0,nuframe1-1 do begin
      gd2 = where(frame1 eq uframe1[j],ngd2)
      medam_uframe1[j] = median([arr0[gd2].airmass])
    endfor
    si = sort(medam_uframe1)
    print,'Night = ',strtrim(unights[i],2),'  ',strtrim(range(medam_uframe1),2),' ',strtrim(min(medam_uframe1),2),' ',strtrim(max(medam_uframe1),2),$
          '  ',strtrim(medam_uframe1[si],2)
  endfor

stop
  
  ; 1.) Initial solution for all data combined
  ;-------------------------------------------
  print,'' & print,'1.) Initial global solution for all data'
  undefine,trans0,tlines0
  inparr = arr
  stdred_transphot_smash,'',stp=stp,plotresid=plotresid,$
                 yrange=yrange,fixam=fixam,fixcolr=fixcolr,fixac=fixac,fixcolsq=fixcolsq,$
                 fitac=fitac,fitcolsq=fitcolsq,trans=trans0,$
                 tlines=tlines0,rms=rms0,errlim=errlim,$
                 inparr=inparr,arr=arr0,/nooutput,/silent
  colterm0 = trans0[0].colterm
  colsqterm0 = trans0[0].colsqterm
  print,'RMS      = '+string(rms0,format='(F9.6)')
  fmt2 = '(F7.4)'
  for i=0,n_elements(trans0)-1 do $
    print,'Night '+strtrim(trans0[i].night,2)+' zpoint = '+string(trans0[i].zpterm,format=fmt2)+$
             ' ('+string(trans0[i].zperr,format=fmt2)+')'
  print,'Airmass term   = '+string(trans0[0].amterm,format=fmt2)+' ('+string(trans0[0].amerr,format=fmt2)+')'
  print,'Color term     = '+string(trans0[0].colterm,format=fmt2)+' ('+string(trans0[0].colerr,format=fmt2)+')'
  print,'Am*Color term  = '+string(trans0[0].amcolterm,format=fmt2)+' ('+string(trans0[0].amcolerr,format=fmt2)+')'
  print,'ColorSq term   = '+string(trans0[0].colsqterm,format=fmt2)+' ('+string(trans0[0].colsqerr,format=fmt2)+')'


  
stop
  
  ; 2.) Initial separate solutions for each chip  
  ;---------------------------------------------
  print,'' & print,'2.) Initial solution of each chip separately'
  print,' Chip  Night ZPterm    AMterm   COLterm      RMS'
  undefine,alltrans1,allarr1
  for i=0,nchips-1 do begin
    ;print,'' & print,'Chip = ',strtrim(chips[i],2)
    ind = where(arr.chip eq chips[i],nind)
    inparr = arr[ind]
    undefine,arr1,trans,tlines
    stdred_transphot_smash,'',stp=stp,plotresid=plotresid,$
                 yrange=yrange,fixam=fixam,fixcolr=fixcolr,fixac=fixac,fixcolsq=colsqterm0,$
                 fitac=fitac,trans=trans,$
                 tlines=tlines,rms=rms,errlim=errlim,$
                 arr=arr1,inparr=inparr,/nooutput,/silent
    ; the err limit has been applied to ARR1

;print,'chip = ',strtrim(i+1,2)
;;yr1 = minmax(arr1.mag-arr1.cmag-trans[0].colterm*arr1.col)
;yr1 = [-0.64,-0.1]
;plotc,arr1.airmass,arr1.mag-arr1.cmag-trans[0].colterm*arr1.col,arr1.night,ps=1,$
;      xr=[1,2.4],yr=yr1,xs=1,ys=1,pos=[0.08,0.55,0.98,0.90],colpos=[0.08,0.97,0.98,0.99]
;oplot,[0,3],poly([0,3],[mean(trans.zpterm),trans[0].amterm])
;xyouts,mean([1.0,2.4]),yr1[0]+0.05*range(yr1),'chip '+strtrim(i+1,2)+' RMS='+stringize(rms,ndec=4)+$
;       ' amterm='+stringize(trans[0].amterm,ndec=4),align=0.5,charsize=1.3
;;yr2 = minmax(arr1.mag-arr1.cmag-trans[0].amterm*arr1.airmass)
;yr2 = [-0.9,-0.2]
;plotc,arr1.col,arr1.mag-arr1.cmag-trans[0].amterm*arr1.airmass,arr1.airmass,ps=1,$
;      xr=[-0.4,2.5],yr=yr2,xs=1,ys=1,pos=[0.08,0.05,0.98,0.40],colpos=[0.08,0.47,0.98,0.49],/noerase
;oplot,[-1,5],poly([-1,5],[mean(trans.zpterm),trans[0].colterm])
;xyouts,mean([-0.4,2.5]),yr2[0]+0.05*range(yr2),'colterm='+stringize(trans[0].colterm,ndec=4),align=0.5,charsize=1.3
;stop
    push,allarr1,arr1
    add_tag,trans,'chip',chips[i],trans
    PUSH,alltrans1,trans
    for j=0,n_elements(trans)-1 do $
       print,trans[j].chip,trans[j].night,trans[j].zpterm,trans[j].amterm,$
             trans[j].colterm,trans[j].rms,format='(2I5,4F10.5)'
  endfor

stop

  ; 3.) Fit nightly airmass terms
  ;------------------------------
  print,'' & print,'3.) Fitting airmass/extinction term(s)'
  CASE fitam of
  ; Same value for all nights
  1: begin
    print,'FITAM=1  One airmass value for all nights'
    ; Use previously found global value
    amtermarr = fltarr(nnights)+trans0[0].amterm
    amerrarr = fltarr(nnights)+trans0[0].amerr
  end
  ; Nightly variation of airmass
  2: begin
    print,'FITAM=2  Fitting nightly airmass/extinction terms'
    ;medcolterm = median(alltrans1.colterm)
    amtermarr = fltarr(nnights)
    amerrarr = fltarr(nnights)
    For i=0,nnights-1 do begin
      ;ind = where(allarr1.night eq unights[i],nind)
      ind = where(allarr1.mjd eq unights[i],nind)
      arr1 = allarr1[ind]
      resid = arr1.mag-arr1.cmag-colterm0*arr1.col
      gd = where(arr1.rejected eq 0,ngd)
      coef = robust_poly_fitq(arr1[gd].airmass,resid[gd],1)
      ; get parameter errors from POLY_FIT
      coef2 = poly_fit(arr1[gd].airmass,resid[gd],1,measure_errors=arr1[gd].err,sigma=sigma)
      amtermarr[i] = coef[1]
      amerrarr[i] = sigma[1]
      print,'Night ',strtrim(unights[i],2),' amterm=',stringize(coef[1],ndec=4)
      ;print,coef
      ;plot,arr1[gd].airmass,resid[gd],ps=3,xtit='airmass',tit='night '+strtrim(unights[i],2)
      ;oplot,[0,3],poly([0,3],coef),co=250    
      ;stop
    Endfor
  end
  ; FITAM case not supported yet
  else: begin
    print,'FITAM = ',strtrim(fitam,2),' NOT supported yet.'
    return
  endelse
  ENDCASE  ; fitam case

  ; Fit Color Squared Term
  ;  this also sets the GLOBAL COLTERM if FITCOL=1
  if keyword_set(fitcolsq) then begin
    print,'' & print,'Fitting Color Squared term'
    arr_nozpamterm = arr       ; amterm removed
    for i=0,nnights-1 do begin
      ind = where(arr_nozpamterm.night eq unights[i],nind)
      ind2 = where(trans0.night eq unights[i],nind2)
      arr_nozpamterm[ind].mag -= amtermarr[i]*arr_nozpamterm[ind].airmass
      arr_nozpamterm[ind].mag -= trans0[ind2[0]].zpterm
    endfor
    ; Keep ZP=0, AM=0 and fit only COLR and COLSQ
    arr_nozpamterm.night = 1  ; no chip or nightly variations
    arr_nozpamterm.chip = 1
    inparr = arr_nozpamterm
    stdred_transphot_smash,'',stp=stp,plotresid=plotresid,$
                     yrange=yrange,fixam=0.0,fixac=0.0,fitcolr=1,fitcolsq=1,$
                     trans=trans_colsq,tlines=tlines,rms=rms,errlim=errlim,$
                     arr=arr2,inparr=inparr,/nooutput,/silent
    colterm = trans_colsq.colterm
    colerr = trans_colsq.colerr
    colsqterm = trans_colsq.colsqterm
    colsqerr = trans_colsq.colsqerr
    print,'Color term   = ',strtrim(colterm,2),' +/- ',strtrim(colerr,2)
    print,'ColorSq term = ',strtrim(colsqterm,2),' +/- ',strtrim(colsqerr,2)
  endif else begin
    colterm = colterm0
    colerr = trans0[0].colerr
    colsqterm = 0.0
    colsqerr = 0.0
 endelse
    
stop
  
  ; 4.) Remove airmass terms and fit chip-specfic color-terms and chip/night
  ;       specific zero-point terms
  ;-------------------------------------------------------------------------
  print,'' & print,'4.) Solving for chip-specific color and zero-point terms'
  CASE fitcolr of
  ; One color term for all chips and nights
  1: begin
    print,'FITCOLR=1  Single color-term for all chips/nights'
    arr_noamterm = arr            ; amterm removed
    for i=0,nnights-1 do begin
      ind = where(arr_noamterm.night eq unights[i],nind)
      arr_noamterm[ind].mag -= amtermarr[i]*arr_noamterm[ind].airmass
    endfor
    print,' Chip  Night ZPterm    AMterm   COLterm      RMS'
    undefine,alltrans2,allarr2
    for i=0,nchips-1 do begin
      ind = where(arr_noamterm.chip eq chips[i],nind)
      inparr = arr_noamterm[ind]
      undefine,arr2
      stdred_transphot_smash,'',stp=stp,plotresid=plotresid,$
                   yrange=yrange,fixam=0.0,fixcolr=colrterm,fixac=fixac,fixcolsq=colsqterm,$
                   fitac=fitac,trans=trans,$
                   tlines=tlines,rms=rms,errlim=errlim,$
                   arr=arr2,inparr=inparr,/nooutput,/silent
      push,allarr2,arr2
      add_tag,trans,'chip',chips[i],trans
      PUSH,alltrans2,trans
      for j=0,n_elements(trans)-1 do $
         print,trans[j].chip,trans[j].night,trans[j].zpterm,trans[j].amterm,$
               trans[j].colterm,trans[j].rms,format='(2I5,4F10.5)'
    endfor
    ; ALLARR2 has the amterm REMOVED
  end
  ; Separate color term for each chip
  2: begin
    print,'FITCOLR=2  Chip-specific color terms'
    print,'Solving for chip-specific color and zero-point terms'
    arr_noamterm = arr            ; amterm removed
    for i=0,nnights-1 do begin
      ind = where(arr_noamterm.night eq unights[i],nind)
      arr_noamterm[ind].mag -= amtermarr[i]*arr_noamterm[ind].airmass
    endfor
    print,' Chip  Night ZPterm    AMterm   COLterm      RMS'
    undefine,alltrans2,allarr2
    for i=0,nchips-1 do begin
      ind = where(arr_noamterm.chip eq chips[i],nind)
      inparr = arr_noamterm[ind]
      undefine,arr2
      stdred_transphot_smash,'',stp=stp,plotresid=plotresid,$
                   yrange=yrange,fixam=0.0,fixcolr=fixcolr,fixac=fixac,fixcolsq=colsqterm,$
                   fitac=fitac,trans=trans,$
                   tlines=tlines,rms=rms,errlim=errlim,$
                   arr=arr2,inparr=inparr,/nooutput,/silent
      push,allarr2,arr2
      add_tag,trans,'chip',chips[i],trans
      PUSH,alltrans2,trans
      for j=0,n_elements(trans)-1 do $
         print,trans[j].chip,trans[j].night,trans[j].zpterm,trans[j].amterm,$
               trans[j].colterm,trans[j].rms,format='(2I5,4F10.5)'
    endfor
    ; ALLARR2 has the amterm REMOVED
  end
  ; FITCOLR case not supported yet
  else: begin
    print,'FITCOLR = ',strtrim(fitcolr,2),' NOT supported yet.'
    return
  endelse
  ENDCASE  ; fitcolr case

stop
  
  ; 5.) Measure nightly zero-point term
  ;------------------------------------
  print,'' & print,'5.) Fitting nightly zero-point term'
  zptermarr = fltarr(nnights)
  zperrarr = fltarr(nnights)
  for i=0,nnights-1 do begin
    ind = where(allarr2.night eq unights[i],nind)
    arr2 = allarr2[ind]
    resid = arr2.mag-arr2.cmag
    for j=0,nchips-1 do begin
      indch = where(arr2.chip eq chips[j],nindch)
      indtrans = where(alltrans2.night eq unights[i] and alltrans2.chip eq chips[i],nindtrans)
      resid[indch] -= alltrans2[indtrans[0]].colterm*arr2[indch].col
    endfor
    gd = where(arr2.rejected eq 0,ngd)
    zpterm = median(resid[gd])
    ; get parameter errors from POLY_FIT
    coef = poly_fit(arr2[gd].airmass,resid[gd],0,measure_errors=arr2[gd].err,sigma=sigma)
    zptermarr[i] = zpterm
    zperrarr[i] = sigma[0]
    print,'Night ',strtrim(unights[i],2),' zpterm=',stringize(zpterm,ndec=4)
  endfor


  
stop
  
  ; 6.)  Remove airmass terms and nightly zpterm and fit for
  ;       chip-specific zero-point and color terms
  ;---------------------------------------------------------
  print,'' & print,'6.) Solving for chip-specific zero-point and color terms'
  tarr = arr            ; amterm and nightly zpterm removed
  for i=0,nnights-1 do begin
    ind = where(tarr.night eq unights[i],nind)
    tarr[ind].mag -= amtermarr[i]*tarr[ind].airmass
    tarr[ind].mag -= zptermarr[i]
  endfor
  if fitcolr eq 1 then fixcolr=colterm0  ; global color term
  print,' Chip  Night ZPterm    AMterm   COLterm      RMS'
  undefine,alltrans3,allarr3
  for i=0,nchips-1 do begin
    ind = where(tarr.chip eq chips[i],nind)
    inparr = tarr[ind]
    inparr.night = min(inparr.night)  ; solve as ONE night, one solution
    undefine,arr3,fixcolr1
    if fitcolr eq 1 then fixcolr1=colterm   ; use global color term
    stdred_transphot_smash,'',stp=stp,plotresid=plotresid,$
                 yrange=yrange,fixam=0.0,fixcolr=fixcolr1,fixac=fixac,fixcolsq=colsqterm,$
                 fitac=fitac,trans=trans,$
                 tlines=tlines,rms=rms,errlim=errlim,$
                 arr=arr3,inparr=inparr,/nooutput,/silent
    push,allarr3,arr3
    add_tag,trans,'chip',chips[i],trans
    PUSH,alltrans3,trans
    for j=0,n_elements(trans)-1 do $
       print,trans[j].chip,trans[j].night,trans[j].zpterm,trans[j].amterm,$
             trans[j].colterm,trans[j].rms,format='(2I5,4F10.5)'
;if trans.zperr lt 0.00001 then stop
  endfor
  ; ALLARR3 has the amterm and nightly zpterm REMOVED
  ; Sometimes ZPERR is 0, MPFIT problem, set this to a minium value
  gdtrans = where(alltrans3.zperr gt 0,ngdtrans,comp=bdtrans,ncomp=nbdtrans)
  if nbdtrans gt 0 then alltrans3[bdtrans].zperr = min(alltrans3[gdtrans].zperr)

  
stop
  
  ; 7.) Construct final transformation equations
  ;---------------------------------------------
  temp = alltrans3[0]
  struct_assign,{hello:''},temp
  alltrans = replicate(temp,nnights*nchips)  ; not all chip/night combinations might be represented
  ;alltrans = alltrans1            ; initialize
  allarr = allarr1
  for i=0,nnights-1 do begin
    for j=0,nchips-1 do begin
      ;indtrans = where(alltrans.night eq unights[i] and alltrans.chip eq chips[j],nindtrans)
      indtrans = i*nchips+j

      ; start with final chip-specific color/zpterm and errors
      indtrans3 = where(alltrans3.chip eq chips[j],nindtrans3)
      trans = alltrans3[indtrans3]
      ; Add nightly zpterm
      trans.zpterm += zptermarr[i]
      trans.zperr = sqrt(trans.zperr^2 + zperrarr[i]^2)  ; add in nightly zpterm error
      ; Add nightly amterm
      trans.amterm += amtermarr[i]
      trans.amerr = amerrarr[i]
      ; Nightly color term
      if fitcolr eq 1 then begin
         trans.colterm = colterm
         trans.colerr = colerr
      endif
      ; Add am*col term
      trans.amcolterm = trans0[0].amcolterm
      trans.amcolerr = trans0[0].amcolerr
      ; Add col*col term
      trans.colsqterm = colsqterm
      trans.colsqerr = colsqerr
      ; Night number
      trans.night = unights[i]
      ; Chip number
      trans.chip = chips[j]

; AC and colsq terms!!!!
      
      ; Calculate final residuals and rms
      ind = where(allarr1.night eq unights[i] and allarr1.chip eq chips[j],nind)
      if nind gt 0 then begin
        arr1 = allarr1[ind]
        arr1.resid = arr1.mag-arr1.cmag-trans.zpterm-$
                     trans.amterm*arr1.airmass-trans.colterm*arr1.col-trans.colsqterm*arr1.col*arr1.col
        allarr[ind] = arr1
        gd = where(arr1.rejected eq 0,npts)
        totwt = total(arr1[gd].weight)
        rms = sqrt( total(arr1[gd].weight*arr1[gd].resid^2.)*npts/((npts-1.)*totwt) )
        trans.rms = rms
      endif else trans.rms=!values.f_nan
        
      ; Put in final array
      alltrans[indtrans] = trans
;stop
    endfor
  endfor

  ; Calculate median offset for each frame
  print,' ' & print,'Statistics for each frame'
  print,'NUM   Frame    MedResid    RMS'
  dum = strsplitter(allarr.frame,'_',/extract)
  frames = reform(dum[0,*])
  ui = uniq(frames,sort(frames))
  uframe = frames[ui]
  nframe = n_elements(uframe)
  print,strtrim(nframe,2),' frames'
  medresid_frame = fltarr(nframe)
  rmsresid_frame = fltarr(nframe)
  for i=0,nframe-1 do begin
    ind = where(frames eq uframe[i] and allarr.rejected eq 0,nind)
    medresid_frame[i] = median([allarr[ind].resid])
    totwt = total(allarr[ind].weight)
    rms = sqrt( total(allarr[ind].weight*allarr[ind].resid^2.)*nind/((nind-1.)*totwt) )
    if nind eq 1 then rms=!values.f_nan
    rmsresid_frame[i] = rms
    print,i+1,uframe[i],medresid_frame[i],rms,format='(I4,A10,2F10.4)'
  endfor
  print,''

  ; Calculate final RMS
  gd = where(allarr.rejected eq 0,npts)
  npts = n_elements(tarr.mag)
  totwt = total(allarr[gd].weight)
  rms = sqrt( total(allarr[gd].weight*allarr[gd].resid^2.)*npts/((npts-1.)*totwt) )

  ; Write to file

  ; Print out the transformation equations
  ;---------------------------------------
  fmt2 = '(F7.4)'
  undefine,lines
  push,lines,'#'
  push,lines,'# FINAL TRANSFORMATION COEFFICIENTS:'
  push,lines,'# Final RMS      = '+string(rms,format='(F9.6)')
  push,lines,'#'
  ;for i=0,nnights-1 do begin
  ;  push,lines,'Night '+strtrim(unights[i],2)+' zpoint = '+string(fpar[i],format=fmt2)+$
  ;             ' ('+string(sigpar[i],format=fmt2)+')'
  ;end
  for i=0,nnights-1 do begin
    push,lines,'# Night '+strtrim(unights[i],2)
    push,lines,'# Airmass term   = '+string(amtermarr[i],format=fmt2)+$
               '  ('+string(amerrarr[i],format=fmt2)+')'
    push,lines,'# ZP term        = '+string(zptermarr[i],format=fmt2)+$
               '  ('+string(zperrarr[i],format=fmt2)+')'
  endfor
  push,lines,'# Am*Color term  = '+string(alltrans[0].amcolterm,format=fmt2)+$
       '  ('+string(alltrans[0].amcolerr,format=fmt2)+')'
  push,lines,'# ColorSq term   = '+string(alltrans[0].colsqterm,format=fmt2)+$
       '  ('+string(alltrans[0].colsqerr,format=fmt2)+')'
  if not keyword_set(silent) then printline,lines

  ; Figure out what the magnitude and color names are
  tags = TAG_NAMES(arr)
  magname = strlowcase(tags[6])  ; lowercase, ugriz
  colname = strlowcase(tags[8])  ; lowercase, ugriz
  colname = REPSTR(colname,'_','-')


  ; Print transformation equations to file
  ;---------------------------------------
  ;    M    M-T  -0.9990    0.1402     -0.1345    0.0000   0.0000
  ;              1.094E-02  5.037E-03  2.010E-03  0.0000   0.0000

  ; Final transformation equations
  if not keyword_set(nooutput) then $
    WRITELINE,magname+'.trans',lines

  ; Transformation equations for each night and chip
  for i=0,nnights-1 do begin
    print,'Night ',strtrim(unights[i],2)
    undefine,lines1,all
    tlines = strarr(nchips,2)
    for j=0,nchips-1 do begin
      tind = where(alltrans.night eq unights[i] and alltrans.chip eq chips[j],ntind)
      itrans = alltrans[tind[0]]

      undefine,lines1
      add=string(itrans.chip,format='(I3)')+'  '+magname+'  '+colname+'  '
      nadd = strlen(add)
      magline = add+string(itrans.zpterm,format=fmt2)+'   '+string(itrans.amterm,format=fmt2)+$
           '   '+string(itrans.colterm,format=fmt2)+'   '+string(itrans.amcolterm,format=fmt2)+'   '+string(itrans.colsqterm,format=fmt2)
      push,lines1,magline
      spaces = string(byte(lonarr(nadd)+32))
      errline = spaces+string(itrans.zperr,format=fmt2)+'   '+string(itrans.amerr,format=fmt2)+$
                '   '+string(itrans.colerr>0.0001,format=fmt2)+'   '
      if itrans.amcolterm ne 0.0 then errline+=string(itrans.amcolerr>0.0001,format=fmt2)+'   ' else $
         errline+=string(itrans.amcolerr,format=fmt2)+'   '
      if itrans.colsqterm ne 0.0 then errline+=string(itrans.colsqerr>0.0001,format=fmt2) else $
         errline+=string(itrans.colsqerr,format=fmt2)
      push,lines1,errline

      ; Add to trans lines array
      tlines[i,*] = lines1

      push,all,' '
      push,all,lines1
    endfor
    printline,all
    print,''
    WRITELINE,'n'+strtrim(unights[i],2)+magname+'.trans',all

    ; Add to MAG.trans file
    WRITELINE,magname+'.trans',['','Night '+strtrim(unights[i],2)],/append
    WRITELINE,magname+'.trans',all,/append
  endfor

  ; Output the stars used for deriving the transformation
  ; with the REJECTED tag
  arr = allarr

  ; USE MJD NIGHT VALUES!!!!
  

  
stop
  
  return

endif  ; each chip separately


; Add CMAG, COL to arr
ADD_TAG,arr,'CMAG',0.0,arr
arr.cmag = arr.(6)
ADD_TAG,arr,'COL',0.0,arr
arr.col = arr.(8)

; Add a rejection tag
ADD_TAG,arr,'rejected',0,arr
arr.rejected = 0

; Select only stars with good instrumental and calibrated photometry
gdphot = where(arr.mag gt 0 and arr.mag lt 50 and $
               arr.cmag gt 0 and arr.cmag lt 50 and $
               arr.col gt -5 and arr.col lt 20,ngdphot)
if ngdphot eq 0 then stop,'NO stars'
arr = arr[gdphot]


; Error limit
if n_elements(errlim) gt 0 then if errlim gt 0.0001 then begin
  if not keyword_set(silent) then $
    print,'Imposing Error limit <= ',strtrim(errlim,2)
  ;gderr = where(arr.err le errlim,ngderr)
  gderr = where(arr.err le errlim and arr.(7) le errlim,ngderr) ; instr and cal errors
  if ngderr eq 0 then gderr = where(arr.err le 2*errlim,ngderr)
  if ngderr eq 0 then stop,'no stars'
  if not keyword_set(silent) then $
    print,strtrim(ngderr,2),'/',strtrim(narr,2),' observations retained'
  arr = arr[gderr] 
endif


; Figure out what the magnitude and color names are
tags = TAG_NAMES(arr)
magname = tags[6]
colname = tags[8]
colname = REPSTR(colname,'_','-')


; LO and HI indices for the unique stars
si = sort(arr.id)
arr = arr[si]
idlo = where(shift(arr.id,1) ne arr.id)
idhi = where(shift(arr.id,-1) ne arr.id)


;; Get a number for each unique "real" star
; Use the calibrated color and magnitude to do this.
ADD_TAG,arr,'REALSTAR',0L,arr
ui = uniq(arr.id,sort(arr.id))
uniqid = arr[ui].id
nuniqid = n_elements(ui)
for i=0,nuniqid-1 do begin
  gd = where(arr.id eq uniqid[i],ngd)
  arr[gd].realstar=i+1
end


; Nights
; mapnight[night] gives the parameter zero-point
; index for the "night"
; zeronight = par[mapnight[night]]
ui = uniq(arr.night,sort(arr.night))
unights = arr[ui].night
nnights = n_elements(unights)
mapnight = lonarr(max(unights)+1)-1
mapnight[unights] = indgen(nnights)

; Chips
; mapchip[chip] gives the parameter zero-point
; index for the "chip"
; zerochip = par[mapchip[chip]]
ui = uniq(arr.chip,sort(arr.chip))
uchips = arr[ui].chip
nchips = n_elements(uchips)
mapchip = lonarr(max(uchips)+1)-1
mapchip[uchips] = indgen(nchips)

; Outlier rejection loop
flag=0
count=0
nnewrej=0
nrej=0
;nnights = max(arr.night)-min(arr.night)+1
;par = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
;par = replicate(0.0,nnights+2)
;par = replicate(0.0,nnights+3)
par = replicate(0.0,nnights+4)
if not keyword_set(silent) then begin
  print,'-------------------------------------'
  print,' ITER  NPTS   RMS       SIG     NREJ  '
  print,'====================================='
endif
WHILE (flag ne 1) do begin

  ; Get non-rejected stars
  gd = where(arr.rejected eq 0,ngd)
  tarr = arr[gd]

  ; Fitting with MPFIT.PRO
  func = 'std_devtrans'
  ;fa = {y:arr.mag, mag:arr.cmag, col:arr.col, am:arr.airmass, night:arr.night, weight:arr.weight}
  if tag_exist(arr,'CHIP') then chip=tarr.chip else chip=lonarr(ngd)+1
  fa = {y:tarr.cmag, mag:tarr.mag, col:tarr.col, am:tarr.airmass,$
        night:tarr.night, chip:chip, mapnight:mapnight, weight:tarr.weight}

  npar = n_elements(par)
  parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},npar)
  ; zeropoint
  parinfo[npar-5].limited=1 & parinfo[npar-5].limits=[-5,5]
  ; airmass
  parinfo[npar-4].limited=1 & parinfo[npar-4].limits=[-5,5]
  if n_elements(fixam) gt 0 then begin
    par[npar-4] = fixam
    parinfo[npar-4].fixed = 1
  endif
  ; color
  parinfo[npar-3].limited=1 & parinfo[npar-3].limits=[-5,5]
  if n_elements(fixcolr) gt 0 then begin
    par[npar-3] = fixcolr
    parinfo[npar-3].fixed = 1
  endif
  ; am*col
  parinfo[npar-2].fixed = 1                           ; fix am*col by default
  if keyword_set(fitac) then parinfo[npar-2].fixed=0  ; let am*col float
  if n_elements(fixac) gt 0 then begin
    par[npar-2] = fixac
    parinfo[npar-2].fixed = 1
  endif  
  ; col*col
  parinfo[npar-1].fixed = 1                              ; fix col*col by default
  if keyword_set(fitcolsq) then parinfo[npar-1].fixed=0  ; let col*col float
  if n_elements(fixcolsq) gt 0 then begin
    par[npar-1] = fixcolsq
    parinfo[npar-1].fixed = 1
  endif  
  
  par[0:nnights-1] = median(tarr.mag-tarr.cmag)

  fpar = MPFIT(func,par,functargs=fa,status=status,perror=perror,bestnorm=chisq,$
               parinfo=parinfo,dof=dof,autoderivative=1,ftol=ftol,/quiet)

  sigpar = perror * sqrt(chisq/dof)

  ; total resid
  model = std_transfunc(fpar,mag=tarr.mag,col=tarr.col,am=tarr.airmass,night=tarr.night,mapnight=mapnight)
  ;resid = arr.mag-model
  resid = model-tarr.cmag
  ;rms = sqrt(total(resid^2)/(npts-n_elements(fpar)-1))
  npts = n_elements(tarr.mag)
  totwt = total(tarr.weight)
  rms = sqrt( total(tarr.weight*resid^2.)*npts/((npts-1.)*totwt) )
  sig = mad(resid)
  
  model2 = std_transfunc(fpar,mag=arr.mag,col=arr.col,am=arr.airmass,night=arr.night,mapnight=mapnight)
  resid2 = model2-arr.cmag


  ; New Rejecting stars with bad resids
  thresh = 3.0*sig > 0.015 
  bad1 = where(abs(resid2-median(resid2,/even)) gt thresh and arr.rejected eq 0,nbad1)
  if nbad1 gt 0 then begin
    arr[bad1].rejected = 1
    ;remove,bd,arr
    par = fpar
  endif

  ; Chuck any stars that are REALLY off
  nbad2 = 0
  com=''
  nstar_rej = 0
  for i=0,nuniqid-1 do begin
    nind = idhi[i]-idlo[i]+1
    ind = lindgen(nind)+idlo[i]
    ;ind = where(arr.id eq uniqid[i],nind)
    medresid = median(resid2[ind],/even)
    if (abs(medresid) gt 2.0*sig and min(arr[ind].rejected) eq 0 and nind gt 2) then begin
      newbad = where(arr[ind].rejected eq 0,nbad2)
      arr[ind].rejected = 1
      
      if com eq '' then com='Star(s) ' else com=com+', '
      ;com = com+strtrim(uniqid[i],2)

      nstar_rej++
    endif
  end
  if nstar_rej gt 0 then com=strtrim(nstar_rej,2)+' stars rejected'
  ;if com ne '' then com=com+' rejected'

  count++


  ; How many new rejects
  nnewrej = nbad1+nbad2
  ; Ending?
  if nnewrej eq 0 then flag=1

  fmt = '(I3,I6,F10.5,F10.5,I5,A5,A-45)'
  ;print,format=fmt,count,npts,rms,sig,nbd
  if not keyword_set(silent) then $
    print,format=fmt,count,npts,rms,sig,nnewrej,'',com

  ;stop

endwhile
if not keyword_set(silent) then $
  print,'-------------------------------------'

;nnights = max(arr.night)-min(arr.night)+1


; Put the residual information into the structure
ADD_TAG,arr,'resid',0.0,arr
model = std_transfunc(fpar,mag=arr.mag,col=arr.col,am=arr.airmass,night=arr.night,mapnight=mapnight)
resid = model-arr.cmag
arr.resid = resid


; Print out the number of "good" observations per night
;------------------------------------------------------
if not keyword_set(silent) then begin
  print,''
  print,'Number of "good" observations per night'
  print,'---------------------------------------'
endif
for i=0,nnights-1 do begin
  inight = unights[i]
  dum = where(arr.night eq inight and arr.rejected eq 0,ngd_inight)
  if not keyword_set(silent) then $
    print,'Night '+strtrim(inight,2)+' Nobs = '+strtrim(ngd_inight,2)
end

; Print out RMS per night
;------------------------
if not keyword_set(silent) then begin
  print,''
  print,'RMS per night'
  print,'-------------'
endif
for i=0,nnights-1 do begin
  inight = unights[i]
  igd = where(arr.night eq inight and arr.rejected eq 0,ngd_inight)
  irms = sqrt(median([arr[igd].resid]^2.0))
  if ngd_inight eq 1 then irms = 0.0
  if not keyword_set(silent) then $
    print,'Night '+strtrim(inight,2)+' RMS = '+string(irms,format='(F10.5)')
end
if not keyword_set(silent) then print,''

; Plot the residuals
;---------------------------------
if keyword_set(plotresid) then begin

  ;window,15,xsize=600,ysize=400
  !p.multi = [0,3,2]
  charsize = 2.0
  bd = where(arr.rejected eq 1,nbd)

  dy = range(arr.resid)*0.1
  if keyword_set(yrange) then if n_elements(yrange) eq 2 then yr=yrange
  if n_elements(yr) eq 0 then yr = [min(arr.resid)-dy,max(arr.resid)+dy]

  ; Resid vs. Airmass
  ;------------------
  dx = range(arr.airmass)*0.1
  xr = [min(arr.airmass)-dx,max(arr.airmass)+dx]
  plot,arr.airmass,arr.resid,ps=1,xtit='Airmass',ytit='Residuals',tit='Airmass',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize
  if nbd gt 0 then oplot,[arr[bd].airmass],[arr[bd].resid],ps=1,co=250

  ; Resid vs. Color
  ;----------------
  dx = range(arr.col)*0.1
  xr = [min(arr.col)-dx,max(arr.col)+dx]
  plot,arr.col,arr.resid,ps=1,xtit='Color',ytit='Residuals',tit='Color',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize
  if nbd gt 0 then oplot,[arr[bd].col],[arr[bd].resid],ps=1,co=250

  ; Resid vs. Magnitude
  ;--------------------
  dx = range(arr.cmag)*0.1
  xr = [min(arr.cmag)-dx,max(arr.cmag)+dx]
  plot,arr.cmag,arr.resid,ps=1,xtit='Magnitude',ytit='Residuals',tit='Magnitude',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize
  if nbd gt 0 then oplot,[arr[bd].cmag],[arr[bd].resid],ps=1,co=250

  ; Resid vs. night
  ;----------------
  xr = [min(arr.night)-1,max(arr.night)+1]
  plot,arr.night,arr.resid,ps=1,xtit='Night',ytit='Residuals',tit='Night',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize,xminor=1
  if nbd gt 0 then oplot,[arr[bd].night],[arr[bd].resid],ps=1,co=250

  ; Resid vs. frame
  ;----------------
  ;frame = strarr(narr)
  ;for i=0,narr-1 do frame[i] = first_el(strsplit(arr[i].id,'-',/extract))
  ui = uniq(arr.frame,sort(arr.frame))
  uframes = arr[ui].frame
  nframes = n_elements(uframes)

  if nframes le 60 then begin

    if not keyword_set(silent) then begin
      print,' NUM      FRAME'
      for i=0,nframes-1 do print,format='(I3,A12)',i+1,uframes[i]
    endif

    xtickv = indgen(nframes+1)
    xtickname = [' ',strtrim(indgen(nframes)+1,2),' ']

    xr = [0,nframes+1]
    plot,[0],[0],/nodata,ps=1,xtit='Frame',ytit='Residuals',tit='Frame',$
         xr=xr,yr=yr,xs=1,ys=1,charsize=charsize,xtickinterval=1,$
         xtickv=xtickv,xtickname=xtickname,xminor=1
    for i=0,nframes-1 do begin
      frameind = where(arr.frame eq uframes[i],nframeind)
      oplot,[lonarr(nframeind)+i+1],[arr[frameind].resid],ps=1
      bdframeind = where(arr.frame eq uframes[i] and arr.rejected eq 1,nbdframeind)
      if nbdframeind gt 0 then oplot,[lonarr(nbdframeind)+i+1],[arr[bdframeind].resid],ps=1,co=250
    end

  endif else begin
    print,'Too many frames to plot RESID vs. FRAME'
  endelse

  ; Resid vs. realstar
  ;-------------------
  xr = [0,max(arr.realstar)+1]
  plot,arr.realstar,arr.resid,ps=1,xtit='Star',ytit='Residuals',tit='Star',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize
  if nbd gt 0 then oplot,[arr[bd].realstar],[arr[bd].resid],ps=1,co=250

  ui = uniq(arr.realstar,sort(arr.realstar))
  urealstar = arr[ui].realstar
  nrealstar = n_elements(urealstar)

  if not keyword_set(silent) then begin
    print,''
    print,' NUM  STAR ID'
  endif
  for i=0,nrealstar-1 do begin
    starind = where(arr.realstar eq urealstar[i],nstarind)
    if not keyword_set(silent) then $
      print,format='(I3,'+strtrim(nstarind,2)+'A12)',i+1,arr[starind[0]].id
  end

  !p.multi=0

endif



; Print out the transformation equations
;---------------------------------------
fmt2 = '(F7.4)'
undefine,lines
push,lines,'# '
push,lines,'# FINAL TRANSFORMATION COEFFICIENTS:'
push,lines,'# Final RMS      = '+string(rms,format='(F9.6)')
;push,lines,'Final RMS      = '+string(rms,format=fmt2)
push,lines,'# '
for i=0,nnights-1 do begin
  ;push,lines,'Night '+strtrim(i+min(arr.night),2)+' zpoint = '+string(fpar[i],format=fmt2)+$
  ;           ' ('+string(sigpar[i],format=fmt2)+')'
  push,lines,'# Night '+strtrim(unights[i],2)+' zpoint = '+string(fpar[i],format=fmt2)+$
             ' ('+string(sigpar[i],format=fmt2)+')'
end
push,lines,'# Airmass term   = '+string(fpar[nnights],format=fmt2)+' ('+string(sigpar[nnights],format=fmt2)+')'
push,lines,'# Color term     = '+string(fpar[nnights+1],format=fmt2)+' ('+string(sigpar[nnights+1],format=fmt2)+')'
push,lines,'# Am*Color term  = '+string(fpar[nnights+2],format=fmt2)+' ('+string(sigpar[nnights+2],format=fmt2)+')'
push,lines,'# ColorSq term   = '+string(fpar[nnights+3],format=fmt2)+' ('+string(sigpar[nnights+3],format=fmt2)+')'
if not keyword_set(silent) then printline,lines



; Print transformation equations to file
;---------------------------------------
;    M    M-T  -0.9990    0.1402     -0.1345    0.0000   0.0000
;              1.094E-02  5.037E-03  2.010E-03  0.0000   0.0000

; Final transformation equations
if not keyword_set(nooutput) then $
  WRITELINE,magname+'.trans',lines


; Start the trans structure
;  one element for each night
transdum = {magname:magname,colname:colname,night:0L,rms:rms,zpterm:0.0,zperr:0.0,$
            amterm:fpar[nnights],amerr:sigpar[nnights],colterm:fpar[nnights+1],$
            colerr:sigpar[nnights+1],amcolterm:fpar[nnights+2],amcolerr:sigpar[nnights+2],$
            colsqterm:fpar[nnights+3],colsqerr:sigpar[nnights+3],magline:'',errline:''}
trans = REPLICATE(transdum,nnights)


; Transformation equations for each night
undefine,lines1,all
tlines = strarr(nnights,2)
for i=0,nnights-1 do begin

  undefine,lines1
  ;push,lines1,' '
  ;push,lines1,'Night '+strtrim(i+1,2)+' Transformation Equations'
  add='  '+magname+'  '+colname+'  '
  nadd = strlen(add)
  magline = add+string(fpar[i],format=fmt2)+'   '+string(fpar[nnights],format=fmt2)+$
       '   '+string(fpar[nnights+1],format=fmt2)+'   '+string(fpar[nnights+2],format=fmt2)+'   '+string(fpar[nnights+3],format=fmt2)
  push,lines1,magline
  spaces = string(byte(lonarr(nadd)+32))
  errline = spaces+string(sigpar[i],format=fmt2)+'   '+string(sigpar[nnights],format=fmt2)+$
       '   '+string(sigpar[nnights+1]>0.0001,format=fmt2)+'   '+$
       string(sigpar[nnights+2]>0.0001,format=fmt2)+'   '+string(sigpar[nnights+3]>0.0001,format=fmt2)
  push,lines1,errline

  ; Writing transformation equations for this night only
  ;WRITELINE,'n'+strtrim(i+1,2)+magname+'.trans',lines1
  if not keyword_set(nooutput) then $
    WRITELINE,'n'+strtrim(unights[i],2)+magname+'.trans',lines1

  ; Add to trans structure
  trans[i].night = unights[i]
  trans[i].zpterm = fpar[i]
  trans[i].zperr = sigpar[i]
  trans[i].magline = magline
  trans[i].errline = errline
  ; Add to trans lines array
  tlines[i,*] = lines1

  push,all,' '
  ;push,all,'Night '+strtrim(i+1,2)+' Transformation Equations'
  push,all,'Night '+strtrim(unights[i],2)+' Transformation Equations'
  push,all,lines1
endfor

; Write everything to a log file
if not keyword_set(nooutput) then begin
  WRITELINE,magname+'.trans',all,/append

  if not keyword_set(silent) then begin
    print,''
    print,'Transformation equations written to: ',magname+'.trans'
    print,''
  endif
endif


if keyword_set(stp) then stop

end
