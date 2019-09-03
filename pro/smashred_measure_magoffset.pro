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
;  /silent      Don't print anything to the screen.
;  /verbose     Print lots of information to the screen.
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

pro smashred_measure_magoffset,chstr,allsrc,overlapstr,usecalib=usecalib,silent=silent,verbose=verbose,error=error

undefine,overlapstr,error

nchips = n_elements(chstr)
nallsrc = n_elements(allsrc)

; Not enough inputs
if nchips eq 0 or nallsrc eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_measure_magoffset,chstr,allsrc,overlapstr,usecalib=usecalib,silent=silent,verbose=verbose,error=error'
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

; Initialize overlapstr structure
overlapdata_schema = {expnum1:'',chip1:0L,expnum2:'',chip2:0L,overlap:-1,magoff:99.99,magofferr:9.99,nmatch:-1L,primary:0B}
overlapdata = replicate(overlapdata_schema,nchips*3)
noverlapdata = n_elements(overlapdata)

; Calculate all overlaps and magnitude offsets
;---------------------------------------------
count = 0LL
; Outer chip loop
for j=0,nchips-1 do begin
  ; Inner chip loop
  for k=j+1,nchips-1 do begin
    ; Must be different exposures to overlap
    if chstr[j].expnum ne chstr[k].expnum then begin
      ; use vertices to check for overlap
      ;   use code from printVisitOverlap.py
      overlap = DOPOLYGONSOVERLAP(chstr[j].vertices_ra, chstr[j].vertices_dec, chstr[k].vertices_ra, chstr[k].vertices_dec)
      ; Measure mag offsets
      if overlap eq 1 then begin
        ; Add more elements 
        if noverlapdata lt count+2 then begin
          old = overlapdata
          undefine,overlapdata
          overlapdata = replicate(overlapdata_schema,noverlapdata+(1000L>nchips))
          struct_assign,old,overlapdata,/nozero
          noverlapdata = n_elements(overlapdata)
          undefine,old
        endif
        ; Fill in the information
        overlapdata[count].expnum1 = chstr[j].expnum
        overlapdata[count].chip1 = chstr[j].chip
        overlapdata[count].expnum2 = chstr[k].expnum
        overlapdata[count].chip2 = chstr[k].chip
        overlapdata[count].overlap = overlap
        overlapdata[count].primary = 1
        ; Reverse situation
        overlapdata[count+1].expnum1 = chstr[k].expnum
        overlapdata[count+1].chip1 = chstr[k].chip
        overlapdata[count+1].expnum2 = chstr[j].expnum
        overlapdata[count+1].chip2 = chstr[j].chip
        overlapdata[count+1].overlap = overlap

        ; Get the overlapping sources
        if chstr[j].nsrc gt 0 and chstr[k].nsrc gt 0 then begin
          allind1 = lindgen(chstr[j].nsrc)+chstr[j].allsrcindx
          allind2 = lindgen(chstr[k].nsrc)+chstr[k].allsrcindx
          id1 = allsrc[allind1].fid
          id2 = allsrc[allind2].fid
          MATCH,id1,id2,mind1,mind2,count=nmatch,/sort

        ; One of both of the chips have no sources
        endif else begin
          nmatch = 0
        endelse

        ; We have matching sources
        undefine,magoff,magoffsig,magofferr,ngdmag
        if nmatch gt 1 then begin
          ; Get information for chip 1
          str1 = allsrc[allind1[mind1]]
          ; Get information for chip 2
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
          if ngdmag lt 10 then $   ; not enough points, lower error threshold
            gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.1,ngdmag)
          if ngdmag lt 10 then $   ; not enough points, lower error threshold
            gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.2,ngdmag)
          if ngdmag lt 10 then $   ; not enough points, lower error threshold
            gdmag = where(mag1 lt 50. and mag2 lt 50. and magerr lt 0.5,ngdmag)

          ; Some sources with decent photometry
          if ngdmag gt 0 then begin
            ROBUST_MEAN,magdiff[gdmag],magoffset,magoffsig,sig=magerr[gdmag],numrej=numrej,/usemad
            magofferr = magoffsig/sqrt(ngdmag-numrej)  ; stdev of mean
          endif
        endif else ngdmag=0                    ; enough matches to make a measurement

        ; Save the information in the OVERLAPSTR structure
        if nmatch gt 1 and ngdmag gt 1 then begin
          ;magofferr = magoffsig/sqrt(nmatch)  ; error in the mean
          if magoffsig lt 1e-8 then magofferr=10.0  ; if it's zero that normally means it's pretty bad
          ;magoffsig = magoffsig > 1e-6  ; let lower limit, can sometimes be zero
          overlapdata[count].magoff = magoffset
          overlapdata[count].magofferr = magofferr  ; error in the mean
          overlapdata[count].nmatch = nmatch
          ; Reverse situation
          overlapdata[count+1].magoff = -magoffset
          overlapdata[count+1].magofferr = magofferr  ; error in the mean
          overlapdata[count+1].nmatch = nmatch
          ; Print info
          if keyword_set(verbose) then $
            print,overlapdata[count].expnum1,overlapdata[count].chip1,overlapdata[count].expnum2,overlapdata[count].chip2,overlapdata[count].nmatch,$
                  overlapdata[count].magoff,overlapdata[count].magofferr,format='(A10,I5,A10,I5,I7,F10.4,F10.5)'
          ; Increment by two
          count += 2
        endif ; some matches
      endif ; overlap
    endif  ; different exposures
  endfor ; inner chip loop
endfor ; outer chip loop

; Trim excess elements
overlapdata = overlapdata[0:count-1]
noverlapdata = n_elements(overlapdata)

; Initialize the final output structure
overlapstr = {expnum:strtrim(chstr.expnum,2),chip:long(chstr.chip),noverlap:lonarr(nchips),$
              ind0:lonarr(nchips)-1,ind1:lonarr(nchips)-1,$
              index:lonarr(noverlapdata),revindex:lonarr(noverlapdata),data:overlapdata}
; ind0 and ind1 are the index into "index" and "revindex" which in
; turn are indices into the overlapdata structure,
; e.g. data[index[ind0[i]:ind[i]]] gives the elements of the overlap
; data structure for the ith chip.


; Deal with single exposure situation
uiexp = uniq(chstr.expnum,sort(chstr.expnum))
if n_elements(uiexp) eq 1 then begin
  print,'Only ONE exposure.  No overlaps.'
  overlapstr.data.magoff = 0.0
  overlapstr.data.magofferr = 0.0
  return
endif

; Get reverse indices
oid1 = strtrim(overlapdata.expnum1,2)+'-'+strtrim(overlapdata.chip1,2)
oid2 = strtrim(overlapdata.expnum2,2)+'-'+strtrim(overlapdata.chip2,2)
index_cnt = 0LL
for i=0,nchips-1 do begin
  MATCH,oid1,strtrim(chstr[i].expnum,2)+'-'+strtrim(chstr[i].chip,2),ind1a,ind1b,count=nmatch1,/sort
  ; reverse ones
  MATCH,oid2,strtrim(chstr[i].expnum,2)+'-'+strtrim(chstr[i].chip,2),ind2a,ind2b,count=nmatch2,/sort
  overlapstr.noverlap[i] = nmatch1
  if nmatch1 ne nmatch2 then stop,'DIFFERENT MATCHES!!'
  if nmatch1 gt 0 then begin
    overlapstr.ind0[i] = index_cnt
    overlapstr.ind1[i] = index_cnt+nmatch1-1
    overlapstr.index[index_cnt:index_cnt+nmatch1-1] = ind1a
    overlapstr.revindex[index_cnt:index_cnt+nmatch1-1] = ind2a
    index_cnt += nmatch1
  endif
endfor

; Print out some statistics of the offsets
gd = where(overlapstr.data.magoff lt 50 and overlapstr.data.overlap eq 1 and overlapstr.data.primary eq 1,ngd)
if not keyword_set(silent) then begin
  print,ngd,' mag offsets','med=',median([overlapstr.data[gd].magoff]),'rms=',$
        mad([overlapstr.data[gd].magoff]),format='(I5,A12,A5,F11.6,A5,F11.6)'
  print,'',' errors','med=',median([overlapstr.data[gd].magofferr]),'rms=',$
        mad([overlapstr.data[gd].magofferr]),format='(A5,A12,A5,F11.6,A5,F11.6)'
endif

;stop

end

