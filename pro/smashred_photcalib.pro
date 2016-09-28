;+
;
; SMASHRED_PHOTCALIB
;
; This program calibrates photometry for a field using traditional
; photometric transformation equations, ubercal and additional
; photometric "anchor" data (i.e. 0.9m data).
;
; INPUTS:
;  info    The structure with information on the field.
;  fstr    The structure with information for each exposure.
;  chstr   The structure with information for each chip.
;  allsrc  The structure with information for each source detection.
;  allobj  The structure with information for each unique object.
;  =transfile   The file with the photometric transformation equations.
;  =reduxdir    The reduction directory, the default is "/data/smash/cp/red/photred/"
;  /usegaia  Use GAIA photometry for overall photometric calibration.
;
; OUTPUTS:
;  chstr       The photometric transformation information is added and
;                updated via ubercal.  The MJD and CALIBRATED flags
;                are added.
;  allsrc      The CMAG calibrated magnitudes are updated / calibrated.
;  =error      The error message if one occurred.   
;
; USAGE:
;  IDL>smashred_photcalib,fstr,chstr,allsrc,allobj,/usegaia
;
; D.Nidever  March 2016
;-

pro smashred_photcalib,info,fstr,chstr,allsrc,allobj,transfile=transfile,usegaia=usegaia,$
                       reduxdir=reduxdir,error=error

ninfo = n_elements(info)
nfstr = n_elements(fstr)
nchstr = n_elements(chstr)
nallsrc = n_elements(allsrc)
nallobj = n_elements(allobj)

; Not enough inputs
if ninfo eq 0 or nfstr eq 0 or nchstr eq 0 or nallsrc eq 0 or nallobj eq 0 or n_elements(transfile) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_photcalib,info,fstr,chstr,allsrc,allobj,transfile=transfile,usegaia=usegaia'
  return
endif

; also need the list of photometric nights, transformation equations,
; 0.9m data

; HOW DOES THE 0.9M DATA FIT IN????  Can't be done just as an
; afterburner, since this will affect the mean mags which will in turn
; affect the color terms.  So needs to be a 

; Defaults
if n_elements(reduxdir) eq 0 then reduxdir='/data/smash/cp/red/photred/'
if file_test(reduxdir,/directory) eq 0 then begin
  error = reduxdir+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif
if n_elements(outputdir) eq 0 then outputdir=reduxdir+'catalogs/final/'
if file_test(outputdir,/directory) eq 0 then begin
  if not keyword_set(silent) then print,outputdir+' does NOT exist.  Creating it.'
  FILE_MKDIR,outputdir
endif
; Temporary directory
tmpdir = outputdir+'/tmp/'
if file_test(tmpdir,/directory) eq 0 then FILE_MKDIR,tmpdir
; Gaia directory
gaiadir = '/data/smash/cp/red/photred/gaia/'


; Add MJD to CHSTR
if tag_exist(chstr,'MJD') eq 0 then begin
  add_tag,chstr,'MJD',0L,chstr
  dateobs = chstr.utdate+'T'+chstr.uttime
  for i=0,nchstr-1 do chstr[i].mjd=PHOTRED_GETMJD('','ctio',dateobs=dateobs[i])
endif
; Add CALIBRATED flag to CHSTR
add_tag,chstr,'calibrated',0B,chstr
; Add UBERCAL colums to CHSTR 
add_tag,chstr,'ubercal_magoffset',0.0,chstr
add_tag,chstr,'ubercal_flag',-1L,chstr

; Load the transformation equations
trans_fitstr = MRDFITS(transfile,1,/silent)
trans_chipstr = MRDFITS(transfile,2,/silent)
trans_ntstr = MRDFITS(transfile,3,/silent)

; --- Add the transformation equation information to CHSTR ---

; Add the photometric transformation equation columns to CHSTR
if tag_exist(chstr,'photometric') eq 0 then begin
  newtags = ['photometric','badsoln','band','colband','colsign','zpterm','zptermsig','amterm',$
             'amtermsig','colterm','coltermsig','amcolterm','amcoltermsig','colsqterm','colsqtermsig']
  temp = chstr & undefine,chstr
  add_tags,temp,newtags,['0B','0B','""','""','-1',replicate('0.0d0',10)],chstr
  undefine,temp
endif
chstr.band = chstr.filter

; Transfer over the photometric transformation equations
chstr_mjdchipfilt = strtrim(chstr.mjd,2)+'-'+strtrim(chstr.chip,2)+'-'+strtrim(chstr.filter,2)
trans_mjdchipfilt = strtrim(trans_fitstr.mjd,2)+'-'+strtrim(trans_fitstr.chip,2)+'-'+strtrim(trans_fitstr.filter,2)
for i=0,nchstr-1 do begin
  MATCH,chstr_mjdchipfilt[i],trans_mjdchipfilt,ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    chstr[i].photometric = trans_fitstr[ind2].photometric
    chstr[i].badsoln = trans_fitstr[ind2].badsoln
    chstr[i].band = trans_fitstr[ind2].filter
    chstr[i].colband = trans_fitstr[ind2].colband
    chstr[i].colsign = trans_fitstr[ind2].colsign
    chstr[i].zpterm = trans_fitstr[ind2].zpterm
    chstr[i].zptermsig = trans_fitstr[ind2].zptermerr
    chstr[i].amterm = trans_fitstr[ind2].amterm
    chstr[i].amtermsig = trans_fitstr[ind2].amtermerr
    chstr[i].colterm = trans_fitstr[ind2].colterm
    chstr[i].coltermsig = trans_fitstr[ind2].coltermerr
    chstr[i].amcolterm = trans_fitstr[ind2].colamterm
    chstr[i].amcoltermsig = trans_fitstr[ind2].colamtermerr
    ;chstr[i].colsqterm = trans_fitstr[ind2].colsqterm
    ;chstr[i].colsqtermsig = trans_fitstr[ind2].colsqtermerr
  endif else begin  ; no match
    chstr[i].badsoln = 1  ; no solution
  endelse

  ; BADSOLN or NON-PHOTOMETRIC, zero-out all trans values except for COLTERM
  if chstr[i].badsoln eq 1 or chstr[i].photometric eq 0 then begin
    chstr[i].zpterm = 0.0
    chstr[i].zptermsig = 0.0
    chstr[i].amterm = 0.0
    chstr[i].amtermsig = 0.0
    chstr[i].amcolterm = 0.0
    chstr[i].amcoltermsig = 0.0
    chstr[i].colsqterm = 0.0
    chstr[i].colsqtermsig = 0.0
  endif
endfor


; Double-check that every CHSTR element has the essentials
;  COLBAND, COLSIGN, COLTERM COLTERMSIG
bdchstr = where(chstr.colband eq '',nbdchstr)
chstr_chipfilt = strtrim(chstr.chip,2)+'-'+strtrim(chstr.filter,2)
trans_chipfilt = strtrim(trans_chipstr.chip,2)+'-'+strtrim(trans_chipstr.filter,2)
for i=0,nbdchstr-1 do begin
  ; Get this information from the chip-specific
  ; transformation equation structure
  MATCH,chstr_chipfilt[bdchstr[i]],trans_chipfilt,ind1,ind2,/sort,count=nmatch
  chstr[bdchstr[i]].colband = trans_chipstr[ind2[0]].colband
  chstr[bdchstr[i]].colsign = trans_chipstr[ind2[0]].colsign
  chstr[bdchstr[i]].colterm = trans_chipstr[ind2[0]].colterm
  chstr[bdchstr[i]].coltermsig = trans_chipstr[ind2[0]].coltermerr
endfor

; Copy PHOTOMETRIC and BADSOLN to FSTR
add_tag,fstr,'photometric',0,fstr
add_tag,fstr,'badsoln',0,fstr
for i=0,n_elements(fstr)-1 do begin
  MATCH,fstr[i].expnum,chstr.expnum,ind1,ind2,/sort,count=nmatch
  fstr[i].photometric = chstr[ind2[0]].photometric
  fstr[i].badsoln = chstr[ind2[0]].badsoln
endfor

; Give statistics on how many exposures have calibration information
; (per band) and photometric data
allfilters = ['u','g','r','i','z']
nocalib = intarr(5)
print,'--------------------------'
print,'FILTER   NEXP   NGOODEXP'
print,'=========================='
for i=0,n_elements(allfilters)-1 do begin
  MATCH,allfilters[i],fstr.filter,ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    gdexp = where(fstr[ind2].photometric eq 1 and fstr[ind2].badsoln eq 0,ngdexp)
  endif else ngdexp=0
  if nmatch gt 0 and ngdexp eq 0 then nocalib[i]=1    ; can't calibrate!
  print,allfilters[i],nmatch,ngdexp,format='(A4,I8,I8)'
endfor
print,'--------------------------'

; Some band CANNOT be calibrated, give average zeropoint
bdnocalib = where(nocalib eq 1,nbdnocalib)
if nbdnocalib gt 0 then begin
  print,'!!!! CANNOT CALIBRATE ',strtrim(nbdnocalib,2),' FILTER(S): ',strjoin(allfilters[bdnocalib],' '),' !!!!'
  print,'Giving them average zeropoints for that band'
  for i=0,nbdnocalib-1 do begin
    MATCH,trans_ntstr.filter,allfilters[bdnocalib[i]],ind1,ind2,/sort,count=nmatch
    zpterm = median([trans_ntstr[ind1].zpterm])
    zptermsig = median([trans_ntstr[ind1].zptermsig])
    MATCH,chstr.filter,allfilters[bdnocalib[i]],ind1,ind2,/sort,count=nmatch
    chstr[ind1].zpterm = zpterm
    chstr[ind1].zptermsig = zptermsig
  endfor
endif

; SHOULD I ALSO COMPUTE EXPOSURE-LEVEL APERTURE CORRECTIONS??
; MAYBE FIT A WEAK SPATIAL POLYNOMIAL??

; ----  LOAD the GAIA catalog for this field ----
gaiafile = gaiadir+info[0].field+'_gaia.fits'
if file_test(gaiafile) eq 0 then stop,gaiafile,' NOT FOUND'
print,'Loading GAIA file'
gaia = MRDFITS(gaiafile,1,/silent)
print,strtrim(n_elements(gaia),2),' GAIA sources loaded'

; ---- LOAD the GAIA-SMASH color terms information ----
gaiacolfile = gaiadir+'gaiasmash_colorterms.fits'
if file_test(gaiacolfile) eq 0 then stop,gaiacolfile,' NOT FOUND'
gaiacolstr = MRDFITS(gaiacolfile,1,/silent)

; Get unique filters
ui = uniq(fstr.filter,sort(fstr.filter))
ufilter = fstr[ui].filter
nufilter = n_elements(ufilter)

; Iterate until we converge
;--------------------------
lastcmag = allsrc.cmag
doneflag = 0
niter = 0
last_zpterm = chstr.zpterm
last_maxdiffcmag = 999.9
last_rmsdiffcmag = 999.9
WHILE (doneflag eq 0) do begin

  print & print,'Global field photometric calibration.  Iteration = ',strtrim(niter+1,2) & print

  ;----------------------------------------------------------
  ; Step 1. Regular photometric calibration with trans eqns.
  ;==========================================================
  ; "regular" trans calib, use allobj for the other band
  ; needed to construct the color, but use the magnitude
  ; from that exposure for the other band to construct
  ; the color (just like photcalib.pro does), otherwise
  ; I'd have to run averagemag many times
  ; This also averages mags for allobj to get the color terms
  print,'--- Step 1. Regular photometric calibration with trans eqns. ---'

  ; Do regular calibration with transformation equations
  ;  for non-photometric data only color-terms are applied
  ; this updates CMAG/CERR in ALLSRC
  SMASHRED_APPLY_PHOTTRANSEQN,fstr,chstr,allsrc,allobj

  ;----------------------------------------------------------
  ; Step 2. Ubercal measurement and offsets to individual
  ;           exposure zeropoints
  ;==========================================================
  ; ubercal measurement and application, with filter loop
  ; add these offsets to the zpterm terms for each chip
  ; in the TRANS structure so they'll improve the photometry
  ; calculated in step 1 the next time around
  ; use trans.photometric to help anchor the offsets
  print,'--- Step 2. Ubercal measurement and application ---'

  ; Filter loop
  undefine,overlapstr,ubercalstr
  For f=0,nufilter-1 do begin
    ifilter = ufilter[f]
    chfiltind = where(chstr.filter eq ifilter,nchfiltind)
    chfiltstr = chstr[chfiltind]
    print,'- ',strtrim(f+1,2),' FILTER = ',ifilter,' ',strtrim(nchfiltind,2),' nchips'

    ; Step 2a. Measure photometric offsets between chip pairs
    print,'2a. Measuring relative magnitude offsets between chip pairs'
    SMASHRED_MEASURE_MAGOFFSET,chfiltstr,allsrc,overlapstr,/usecalib

    ; Step 2b. Determine relative magnitude offsets per chip using ubercal 
    print,'2b. Determine relative magnitude offsets per chip using ubercal'
    MATCH,trans_chipstr.filter,ifilter,ind1,ind2,/sort
    chiptrans = trans_chipstr[ind1]
    SMASHRED_SOLVE_UBERCAL,overlapstr,ubercalstr,chiptrans=chiptrans

    ; Step 2c. Set absolute zeropoint of magnitute offsets with
    ;            photometric data and/or anchor points
    print,'2c. Set absolute zeropoint with photometric data and achor points'
    SMASHRED_SET_ZEROPOINTS,chfiltstr,ubercalstr,gaia,allobj,gaiacolstr,usegaia=usegaia

    ; Step 2d. Set CALIBRATED bit
    print,'2d. Set CALIBRATED bit in CHSTR'
    SMASHRED_SET_CALIBRATED_BIT,chfiltstr,overlapstr

    ; Stuff it back in
    chstr[chfiltind] = chfiltstr

    ; Save the overlapstr and ubercalstr structure
    ;  the ubercalstr mag offset are relative
    ;  see chstr.ubercal_magoffset for absolute ones
    outfile = tmpdir+info[0].field+'_'+ufilter[f]+'_photcalib.dat'
    SAVE,chfiltstr,overlapstr,ubercalstr,file=outfile
 
  Endfor  ; filter loop

  ; BRIGHTER-FATTER CORRECTION??
  ; maybe make it a (instrumental) magnitude-dependent correction, bfcorr
  ; compare psf to aperture photometry to figure it out
  ; BUT CORRECT IT TO WHAT???  what's the fiducial value.
  ; it will likely depend on the seeing, so we might need to figure it
  ; out on an exposure-by-exposure level

  ; Check for convergence
  ;  check changes in ALLSRC CMAG values, changes in CHSTR ZPTERM and iterations
  maxiter = 10   ;20  ;50
  maxdiff_thresh = 0.001    ;0.0001
  maxzpterm_thresh = 0.001  ;0.0001
  diffcmag = abs(allsrc.cmag-lastcmag)
  diffzpterm = abs(chstr.zpterm-last_zpterm)
  maxdiffcmag = max(diffcmag)
  maxdiffzpterm = max(diffzpterm)
  rmsdiffcmag = stddev(diffcmag)
  if (niter ge 4 and (maxdiffcmag lt maxdiff_thresh) and (maxdiffzpterm lt maxzpterm_thresh)) or $       ; small differences 
     (niter ge 4 and (maxdiffcmag gt last_maxdiffcmag) and (rmsdiffcmag gt last_rmsdiffcmag)) or $       ; diffs are getting bigger, plateau
     (niter ge maxiter) then doneflag=1                                                                  ; max iterations
;if doneflag eq 1 then stop,'finished'
  last_zpterm = chstr.zpterm  ; save for next time
  last_maxdiffcmag = maxdiffcmag
  last_rmsdiffcmag = rmsdiffcmag
  lastcmag = allsrc.cmag


  print,''
  print,'CMAG diff  ','max=',max(diffcmag),'med=',median(diffcmag),'rms=',stddev(diffcmag),format='(A12,A5,F11.6,A5,F11.6,A5,F11.6)'
  print,'ZPTERM diff','max=',max(diffzpterm),'med=',median(diffzpterm),'rms=',stddev(diffzpterm),format='(A12,A5,F11.6,A5,F11.6,A5,F11.6)'

  niter++

ENDWHILE

; Apply the transformation equations and do average mags one last time
print,'Applying transformation equations one last time'
SMASHRED_APPLY_PHOTTRANSEQN,fstr,chstr,allsrc,allobj


end
