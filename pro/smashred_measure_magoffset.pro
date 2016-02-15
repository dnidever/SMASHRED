pro smashred_measure_magoffset,chstr,overlapstr,silent=silent

; Measure magnitude offset between sets of exposures to be used later
; to perform ubercal calibration.

; Add medoff tag
add_tag,chstr,'magoffset',0.0,chstr

sz = size(chstr)
nchips = sz[1]

; Calculate all overlaps and magnitude offsets
overlapstr = replicate({expnum1:'',chip1:0L,expnum2:'',chip2:0L,overlap:-1,magoff:99.99,magofferr:9.99,nmatch:-1L},nchips,nchips)
for j=0,nchips-1 do begin
  overlapstr[j,*].expnum1 = chstr[j].expnum
  overlapstr[j,*].chip1 = chstr[j].chip
  overlapstr[*,j].expnum2 = chstr[j].expnum
  overlapstr[*,j].chip2 = chstr[j].chip
endfor
for j=0,nchips-1 do begin
  for k=j+1,nchips-1 do begin
    ; must be different exposures to overlap
    if chstr[j].expnum ne chstr[k].expnum then begin
      ; use vertices to check for overlap
      ;   use code from printVisitOverlap.py
      overlap = DOPOLYGONSOVERLAP(chstr[j].vertices_ra, chstr[j].vertices_dec, chstr[k].vertices_ra, chstr[k].vertices_dec)
      overlapstr[j,k].overlap = overlap
      overlapstr[k,j].overlap = overlap
      ; measure mag offsets
      if overlap eq 1 then begin
        str1 = *chstr[j].data
        str2 = *chstr[k].data
        dcr = 0.5
        ;  Match the two catalogs and get photometric offsets
        PHOTMATCH,str1,str2,ind1,ind2,dcr=dcr,magoffset=magoffset,magoffsig=magoffsig,astoffset=astoffset,count=nmatch,/silent
        ; too few matches, increase matching radius
        ;if nmatch lt 5 then $
        ;  photmatch,str1,str2,ind1,ind2,dcr=1.0,magoffset=magoffset,magoffsig=magoffsig,astoffset=astoffset,count=nmatch,/silent
        ;if nmatch lt 5 then $
        ;  photmatch,str1,str2,ind1,ind2,dcr=1.5,magoffset=magoffset,magoffsig=magoffsig,astoffset=astoffset,count=nmatch,/silent
        if nmatch gt 1 and n_elements(magoffset) gt 0 then begin
          magofferr = magoffsig/sqrt(nmatch)  ; error in the mean
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

