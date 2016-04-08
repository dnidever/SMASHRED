;+
;
; SMASHRED_MEASURE_MAGOFFSET
;
; This program measures photometric offsets between pairs of
; overlapping chip exposures.
;
; INPUTS:
;  chstr       Structure with information on each individual chip.
;  allsrc      Structure with information on each source detection.
;  /usecalib   Use the calibrated photometry (CMAG/CERR).  The default
;                is to use the instrumental photometry (MAG/ERR).
; /silent      Don't print anything to the screen.
;
; OUTPUTS:
;  overlapstr  Structure with overlap and magnitude offset information.
;  =error      The error message if one occurred.
;
; USAGE:
;  IDL>smashred_measure_magoffset,chstr,allsrc,overlapstr
;
; By D.Nidever  March 2016
;-

pro smashred_measure_magoffset,chstr,allsrc,overlapstr,usecalib=usecalib,silent=silent,error=error

undefine,overlapstr,error

nchips = n_elements(chstr)
nallsrc = n_elements(allsrc)

; Not enough inputs
if nchips eq 0 or nallsrc eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_measure_magoffset,chstr,allsrc,overlapstr,usecalib=usecalib,silent=silent,error=error'
  return
endif

; Use calibrated magnitudes, check that we have them
if keyword_set(usecalib) then begin
  srctags = tag_names(allsrc)
  dum = where(srctags eq 'CMAG',ncmag)
  dum = where(srctags eq 'CERR',ncerr)
  if ncmag eq 0 or ncerr eq 0 then begin
    error = '/USECALIB set but NO CMAG/CERR columsn in ALLSRC'
    if not keyword_set(silent) then print,error
    return
  endif
endif

; Add medoff tag
;add_tag,chstr,'magoffset',0.0,chstr

; Calculate all overlaps and magnitude offsets
overlapstr = replicate({expnum1:'',chip1:0L,expnum2:'',chip2:0L,overlap:-1,magoff:99.99,magofferr:9.99,nmatch:-1L},nchips,nchips)
for j=0,nchips-1 do begin
  overlapstr[j,*].expnum1 = chstr[j].expnum
  overlapstr[j,*].chip1 = chstr[j].chip
  overlapstr[*,j].expnum2 = chstr[j].expnum
  overlapstr[*,j].chip2 = chstr[j].chip
endfor

; Outer chip loop
for j=0,nchips-1 do begin
  ; Inner chip loop
  for k=j+1,nchips-1 do begin
    ; Must be different exposures to overlap
    if chstr[j].expnum ne chstr[k].expnum then begin
      ; use vertices to check for overlap
      ;   use code from printVisitOverlap.py
      overlap = DOPOLYGONSOVERLAP(chstr[j].vertices_ra, chstr[j].vertices_dec, chstr[k].vertices_ra, chstr[k].vertices_dec)
      overlapstr[j,k].overlap = overlap
      overlapstr[k,j].overlap = overlap
      ; Measure mag offsets
      if overlap eq 1 then begin
        ; Get the overlapping sources
        id1 = allsrc[chstr[j].allsrcindx:chstr[j].allsrcindx+chstr[j].nsrc-1].id
        id2 = allsrc[chstr[k].allsrcindx:chstr[k].allsrcindx+chstr[k].nsrc-1].id
        MATCH,id1,id2,mind1,mind2,count=nmatch,/sort

        ; We have matching sources
        undefine,magoff,magoffsig,magofferr,ngdmag
        if nmatch gt 1 then begin
          ; Get information for chip 1
          allind1 = lindgen(chstr[j].nsrc)+chstr[j].allsrcindx
          str1 = allsrc[allind1[mind1]]
          ; Get information for chip 2
          allind2 = lindgen(chstr[k].nsrc)+chstr[k].allsrcindx
          str2 = allsrc[allind2[mind2]]

          ; Make the measurement
          if keyword_set(usecalib) then begin  ; calibrated photometry
            mag1=str1.cmag & err1=str1.cerr
            mag2=str2.cmag & err2=str2.cerr
          endif else begin                     ; instrumental photometry
            mag1=str1.mag & err1=str1.err
            mag2=str2.mag & err2=str2.err
          endelse
          magdiff = mag1 - mag2
          magerr = sqrt( err1^2.0 + err1^2.0 )
          gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.07,ngdmag)
          if ngdmag lt 10. then $   ; not enough points, lower error threshold
            gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.1,ngdmag)
          if ngdmag lt 10. then $   ; not enough points, lower error threshold
            gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.2,ngdmag)
          if ngdmag lt 10. then $   ; not enough points, lower error threshold
            gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.5,ngdmag)

          ; Some sources with decent photometry
          if ngdmag gt 0 then begin
            ROBUST_MEAN,magdiff[gdmag],magoffset,magofferr,sig=magerr[gdmag],numrej=numrej
            magoffsig = magofferr/sqrt(ngdmag-numrej)  ; stdev of mean
          endif
        endif else ngdmag=0                    ; enough matches to make a measurement

        if nmatch gt 1 and ngdmag gt 1 then begin
          ;magofferr = magoffsig/sqrt(nmatch)  ; error in the mean
          if magoffsig lt 1e-8 then magofferr=10.0  ; if it's zero that normally means it's pretty bad
          ;magoffsig = magoffsig > 1e-6  ; let lower limit, can sometimes be zero
          overlapstr[j,k].magoff = magoffset
          overlapstr[j,k].magofferr = magofferr  ; eror in the mean
          overlapstr[j,k].nmatch = nmatch
          if not keyword_set(silent) then $
            print,overlapstr[j,k].expnum1,overlapstr[j,k].chip1,overlapstr[j,k].expnum2,overlapstr[j,k].chip2,overlapstr[j,k].nmatch,$
                  overlapstr[j,k].magoff,overlapstr[j,k].magofferr,format='(A10,I5,A10,I5,I7,F10.4,F10.5)'
          ; fill in the reverse situation
          overlapstr[k,j].magoff = -magoffset
          overlapstr[k,j].magofferr = magofferr
          overlapstr[k,j].nmatch = nmatch
        endif ; some matches
      endif ; overlap
    endif  ; different exposures
  endfor ; inner chip loop
endfor ; outer chip loop

;stop

end

