;+
;
; SMASHRED_SET_ZEROPOINTS
;
; This determines the absolute zeropint of the magnitude offsets
; using photometric and/or anchor data and applies these to the
; chip-level transformation equations (zeropoints) in CHSTR.
;
; INPUTS:
;  field        The field name, e.g. Field100.
;  chstr        The structure with information on each chip including
;                   the PHOTOMETRIC column and transformation equations.
;  ubercalstr   The structure with information on the ubercal-solved
;                   offsets for each chip.
;  allobj       The structure of average SMASH source values.  Only
;                   needed if /usegaia is set.
;  =gaia        The structure of GAIA sources for this field.  This is
;                   loaded if needed and passed back on subsequent
;                   calls to save time.
;  /usegaia     Use GAIA photometry for zero-point calibration if there
;                   is no other option.
;
; OUTPUTS:
;  The values in the ZPTERM column in CHSTR are updated.
;
; USAGE:
;  IDL>smashred_set_zeropoints,chstr,ubercalstr,allobj
;
; By D. Nidever  April 2016
;-

pro smashred_set_zeropoints,field,chstr,ubercalstr,allobj,usegaia=usegaia,$
                            dir=dir,gaia=gaia

; Not enough inputs
if n_elements(field) eq 0 or n_elements(chstr) eq 0 or n_elements(ubercalstr) eq 0 or n_elements(allobj) eq 0 then begin
  print,'Syntax - smashred_set_zeropoints,field,chstr,ubercalstr,allobj,gaia=gaia,usegaia=usegaia'
  return
endif

filter = chstr[0].filter
if n_elements(dir) eq 0 then dir='/data/smash/cp/red/photred/catalogs/final/'
dir09m = '/data/smash/cp/red/photred/0.9m/'
gaiadir = '/data/smash/cp/red/photred/gaia/'

; Restore the field overlap structure
if file_test(dir+'smash_fieldoverlaps.fits') eq 0 then stop,'SMASH_FIELDOVERLAPS.FITS NOT FOUND!!'
fieldoverlapstr = MRDFITS(dir+'smash_fieldoverlaps.fits',1,/silent)
fieldoverlapstr.field = strtrim(fieldoverlapstr.field,2)  ; remove any trailing spaces
fieldoverlapstr.overlapfield = strtrim(fieldoverlapstr.overlapfield,2)  ; remove any trailing spaces

; Figure out what calibration options are availble for this field/filter
;-----------------------------------------------------------------------
calopt = intarr(4)
; OPT 1: Calibration with transformation equations
calchipind = where(chstr.photometric eq 1 and chstr.badsoln eq 0,ncalchipind)
if ncalchipind gt 0 then calop[0]=1
; OPT 2: Calibration with calibrated overlapping fields
foverlapind = where(strtrim(fieldoverlapstr.field,2) eq field,nfoverlapind)
if nfoverlapind gt 0 then begin
  ; See if we can calibrate this filter with the overlaps
  overlapfilters = ['u','g','r','i','z']
  overlapfiltind = where(overlapfilters eq filter)  ; index for this filter
  ; Indices for overlapping fields
  ofieldind = fieldoverlapstr[foverlapind].overlapindex[0:fieldoverlapstr[foverlapind].noverlap-1]
  calibfilt = fieldoverlapstr[ofieldind].calib[overlapfiltind]  ; calibration flag for this filter
  if total(calibfilt) gt 0 then calopt[1]=1  
endif
; OPT 3: Calibration with 0.9m data
if file_test(dir09m+field+'_phot.fits') eq 1 then calopt[2]=1
; OPT 4: Calibration with Gaia color-color relations
if file_test(gaiadir+field+'_gaia.fits') eq 1 and keyword_set(usegaia) then calopt[3]=1

; Use the first non-zero option
;  they are already in priority order
if total(calopt) gt 0 then caltype=first_el(where(calopt eq 1))+1


; Measure the photometric zero-point offset using the chosen option
;-------------------------------------------------------------------
; ZPMAGOFF is the SUBTRACTIVE magnitude offset to correct the
;  photometric zeropoint
CASE caltype of
  ; 1) Calibration with DECam transformation equations
  ;---------------------------------------------------
  1:  begin
    print,' Setting photometric zeropoint using photometric DECam calibrated data'
    ; Set absolute zeropoint using photometry and/or anchor data
    gdphot = where(chstr.photometric eq 1 and chstr.badsoln eq 0,ngdphot)
    ; Calculate median mag offset for photometric data
    ;  Use median, don't remove outliers, basically want mean of all
    ;  photometric data
    zpmagoff = median(ubercalstr[gdphot].magoff)
    sigmagoff = mad(ubercalstr[gdphot].magoff)
    zpmagofferr = sigmagoff/sqrt(ngdphot)   ; std.dev. of mean
    ; Set the calibration type in CHSTR
    chstr.zpcalibflag = caltype
  end
  ; 2) Calibration with calibrated overlapping fields
  ;--------------------------------------------------
  2:  begin
    print,' Setting photometric zeropoint using overlapping DECam calibrated fields'
    foverlapind = where(strtrim(fieldoverlapstr.field,2) eq field,nfoverlapind)
    overlapfilters = ['u','g','r','i','z']
    overlapfiltind = where(overlapfilters eq filter)  ; index for this filter
    ; Indices for overlapping fields
    ofieldindall = fieldoverlapstr[foverlapind].overlapindex[0:fieldoverlapstr[foverlapind].noverlap-1]
    calibfilt = fieldoverlapstr[ofieldindall].calib[overlapfiltind]  ; calibration flag for this filter
    ; Get fields that are calibrated in this band
    gdfoverlap = where(calibfilt eq 1,ngdfoverlap)
    ofieldind = ofieldindall[gdfoverlap]
    nfoverlap = ngdfoverlap
    ; Load the overlapping field data
    ;   there'll be duplicates, but I think that's okay
    print,'  Loading overlapping fields ('+strtrim(nfoverlap,2)+'): ',strjoin(fieldoverlapstr[ofieldind].field,', ')
    undefine,olapallobj
    for i=0,nfoverlap-1 do begin
      ;print,'    ',strtrim(i+1,2),' Overlapping field: ',fieldoverlapstr[ofieldind[i]].field
      ostr = MRDFITS(dir+fieldoverlapstr[ofieldind[i]].field+'_combined_allobj_bright.fits',1,/silent)
      PUSH,olapallobj,ostr
    endfor
    ; Match them up
    SRCMATCH,allobj.ra,allobj.dec,olapallobj.ra,olapallobj.dec,0.5,ind1,ind2,/sph,count=nmatch
    if nmatch eq 0 then stop,'NO MATCHES FROM OVERLAPPING FIELDS. SOMETHING IS VERY WRONG!!!'
    mallobj = allobj[ind1]   ; matched catalogs
    molapallobj = olapallobj[ind2]
    ; Get the indices for the correct magnitude and error columns
    tags = tag_names(allobj)
    magind = where(tags eq strupcase(filter),nmagind)
    errind = where(tags eq strupcase(filter)+'ERR',nerrind)
    otags = tag_names(olapallobj)
    omagind = where(otags eq strupcase(filter),nomagind)
    oerrind = where(otags eq strupcase(filter)+'ERR',noerrind)
    magdiff = molapallobj.(omagind) - mallobj.(magind)
    errdiff = sqrt( molapallobj.(oerrind)^2 + mallobj.(errind)^2 )
    ; Get good sources
    ;  the allobj_bright.fits catalogs already have most of these cuts
    gdphot1 = where(mallobj.(magind) lt 50 and abs(mallobj.sharp) lt 1 and mallobj.chi lt 4 and $
                   mallobj.(errind) lt 0.05 and molapallobj.(oerrind) lt 0.05,ngdphot1)
    ; Remove outliers
    medmagdiff = median(magdiff[gdphot1])
    sigmagdiff = mad(magdiff[gdphot1])
    gdphot = where(mallobj.(magind) lt 50 and abs(mallobj.sharp) lt 1 and mallobj.chi lt 4 and $
                   mallobj.(errind) lt 0.05 and molapallobj.(oerrind) lt 0.05 and $
                   abs(magdiff-medmagdiff) lt 3*(sigmagdiff>0.01),ngdphot)
    ; Calculate the weighted offset with bootstrap errors
    coef = DLN_POLY_FIT(gdphot*0.0,magdiff[gdphot],0,measure_errors=errdiff[gdphot],sigma=sigma,/bootstrap)
    zpmagoff = -coef[0]
    zpmagofferr = sigma[0]
    ; Set the calibration type in CHSTR
    chstr.zpcalibflag = caltype
    plotc,molapallobj[gdphot].g-molapallobj[gdphot].r,magdiff[gdphot],mallobj[gdphot].(magind),ps=3,xr=[-1,3.5]
  end
  ; 3) Calibration with 0.9m data
  ;------------------------------
  3:  begin
    print,' Setting photometric zeropoint using 0.9m data'
    ; Load the 0.9m photometry file
    str = MRDFITS(dir09m+field+'_phot.fits',1,/silent)
    ; Match them up
    SRCMATCH,allobj.ra,allobj.dec,str.ra,str.dec,0.5,ind1,ind2,/sph,count=nmatch
    if nmatch eq 0 then stop,'NO MATCHES FROM OVERLAPPING FIELDS. SOMETHING IS VERY WRONG!!!'
    mallobj = allobj[ind1]   ; matched catalogs
    mstr = str[ind2]
    ; Get the indices for the correct magnitude and error columns
    tags = tag_names(allobj)
    magind = where(tags eq strupcase(filter),nmagind)
    errind = where(tags eq strupcase(filter)+'ERR',nerrind)
    otags = tag_names(str)
    omagind = where(otags eq strupcase(filter),nomagind)
    oerrind = where(otags eq strupcase(filter)+'ERR',noerrind)
    magdiff = mstr.(omagind) - mallobj.(magind)
    errdiff = sqrt( mstr.(oerrind)^2 + mallobj.(errind)^2 )
    ; Get good sources
    ;  the allobj_bright.fits catalogs already have these cuts
    gdphot1 = where(mallobj.(magind) lt 50 and abs(mallobj.sharp) lt 1 and mallobj.chi lt 4 and $
                    mallobj.(errind) lt 0.05,ngdphot1)
    ; Remove outliers
    medmagdiff = median(magdiff[gdphot1])
    sigmagdiff = mad(magdiff[gdphot1])
    gdphot = where(mallobj.(magind) lt 50 and abs(mallobj.sharp) lt 1 and mallobj.chi lt 4 and $
                   mallobj.(errind) lt 0.05 and abs(magdiff-medmagdiff) lt 3*(sigmagdiff>0.01),ngdphot)
    ; Calculate the weighted offset with bootstrap errors
    coef = DLN_POLY_FIT(gdphot*0.0,magdiff[gdphot],0,measure_errors=errdiff[gdphot],sigma=sigma,/bootstrap)
    zpmagoff = -coef[0]
    zpmagofferr = sigma[0]
    ; Set the calibration type in CHSTR
    chstr.zpcalibflag = caltype
  end
  ; 4) Calibration with Gaia color-color relations
  ;-----------------------------------------------
  4:  begin
    print,' Setting photometric zeropoint using SMASH-GAIA color-color relations'
    ; Load the Gaia data
    if n_elements(gaia) eq 0 then begin
      gaiafile = gaiadir+field+'_gaia.fits'
      if file_test(gaiafile) eq 0 then stop,gaiafile,' NOT FOUND'
      print,'Loading GAIA file'
      gaia = MRDFITS(gaiafile,1,/silent)
      print,strtrim(n_elements(gaia),2),' GAIA sources'
    endif
    ; Load the Gaia color-color relations information
    gaiacolfile = gaiadir+'gaiasmash_colorterms.fits'
    if file_test(gaiacolfile) eq 0 then stop,gaiacolfile,' NOT FOUND'
    gaiacolstr = MRDFITS(gaiacolfile,1,/silent)
    ; Find the offset relative to GAIA
    SMASHRED_MEASURE_GAIA_OFFSET,gaia,allobj,filter,gaiacolstr,zpmagoff,zpmagofferr,/silent
    ; Set the calibration type in CHSTR
    chstr.zpcalibflag = caltype
  end
  ; No calibration option, just remove median
  ;------------------------------------------
  else: begin
    print,' Not calibrating the zeropoint.  Using ALL points to set average zeropoint'
    ; Calculate median mag offset for photometric data
    ;  Use median, don't remove outliers, basically want mean of all
    ;  photometric data
    zpmagoff = median(ubercalstr.magoff)
    zpmagofferr = 0.0
    ; Set the calibration type in CHSTR
    chstr.zpcalibflag = 0
  end
ENDCASE
print,' Removing zeropoint mag offset = ',stringize(zpmagoff,ndec=4),' +/- ',stringize(zpmagofferr,ndec=5),' mag'

; SAVE the zero-point offset values in CHSTR
chstr.zpcalib_magoffset -= zpmagoff
chstr.zpcalib_magofferr = zpmagofferr

; SUBTRACT the photometric zero-point offset from the ubercal magoffsets
ubercalstr.magoff -= zpmagoff

; APPLY the new zeropoints to the chip-level transformation equations
;  we want to ADD the ubercal magnitude offsets to the magnitudes
;  zpterm is subtracted from instrumental mags, so subtract mag offset
chstr.zpterm -= ubercalstr.magoff
chstr.ubercal_magoffset -= ubercalstr.magoff

; Set the ubercal flag
chstr.ubercal_flag = ubercalstr.flag

;stop

end
