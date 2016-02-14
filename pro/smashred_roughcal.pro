pro smashred_roughcal,cat

; Perform rough calibration using SLR and APASS
; input the "final" catalog of merged sources

rar = minmax(cat.ra)
decr = minmax(cat.dec)

; Get APASS data for this field 
;-------------------------------
; get from disk or VizieR?
; /data/smash/apass, in 30 deg blocks of RA
apasslo = lindgen(12)*30
apasshi = lindgen(12)*30+30
lo = first_el(where(apasslo le rar[0]),/last)
hi = first_el(where(apasshi ge rar[1]))
n = hi-lo+1
ind = indgen(n)+lo
undefine,apass
for i=0,n-1 do begin
  apass1 = mrdfits('/data/smash/apass/apassdr7_ra'+strtrim(apasslo[ind[i]],2)+'-'+strtrim(apasshi[ind[i]],2)+'.fits.gz',1)
  push,apass,apass1
endfor
undefine,apass1

; Find the matches
dcr = 0.5
srcmatch,apass.ra,apass.dec,cat.ra,cat.dec,dcr,ind1,ind2,count=nmatch,/sph,domains=4000

; Measure photometric offset in g-band

stop

; Apply color correction (to other bands) using SLR
;---------------------------------------------------

stop

end
