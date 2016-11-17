pro smashred_2cdcuts,str

; Use color-color cuts for star-galaxy separation

; Use b-spline fits to low-err sources and
; then use simple differences in Y to remove "bad" sources
; maybe only use the cuts for regions where the slope is "low".

nstr = n_elements(str)

; this does all of the bspline fitting
make_stellar_locus,str,tstr

; u-g, u-r, u-i vs. g-i  good separation
; g-z, r-z, i-z vs. g-i is okay but still quite a bit of overlap with galaxies
; u-z, g-r, r-i vs. g-i not very useful

color = ['u-g','u-r','u-i','g-z','r-z','i-z']
thresh = [0.3,0.3,0.3,0.13,0.13,0.13]
;color = ['u-g','u-r','u-i']
;thresh = [0.3,0.3,0.3]
ncolor = n_elements(color)

stags = tag_names(str)
ttags = tag_names(tstr)
mask = bytarr(nstr)+1
for i=0,ncolor-1 do begin

  band1 = first_el(strsplit(color[i],'-',/extract))
  band2 = first_el(strsplit(color[i],'-',/extract),/last)

  if band1 eq 'i' or band2 eq 'i' then begin
    if band1 eq 'i' then begin
      tcolind = where(ttags eq strupcase(band2+'i'),ntcolind)
      tcolbin = -tstr.(tcolind)
    endif else begin
      tcolind = where(ttags eq strupcase(band1+'i'),ntcolind)
      tcolbin = tstr.(tcolind)
    endelse
  endif else begin
    tcolind1 = where(ttags eq strupcase(band1+'i'),ntcolind1)
    tcolind2 = where(ttags eq strupcase(band2+'i'),ntcolind2)
    tcolbin = tstr.(tcolind1) - tstr.(tcolind2)
  endelse

  magind1 = where(stags eq strupcase(band1),nmagind1)
  magind2 = where(stags eq strupcase(band2),nmagind2)
  col = str.(magind1) - str.(magind2)

  ; Interpolate the mode color for each source  
  gd = where(str.g lt 50 and str.i lt 50 and str.(magind1) lt 50 and str.(magind2) lt 50 and $
             str.g-str.i ge min(tstr.gibin) and str.g-str.i le max(tstr.gibin),ngd)
  INTERP,tstr.gibin,tcolbin,str[gd].g-str[gd].i,tcol

  coldiff = fltarr(nstr)+99.99
  coldiff[gd] = col[gd] - tcol

  ; Divide by observation errors
  errind1 = where(stags eq strupcase(band1)+'ERR',nerrind1)
  errind2 = where(stags eq strupcase(band2)+'ERR',nerrind2)
  obserr = sqrt( str.(errind1)^2 + str.(errind2)^2 )
  reldiff = coldiff
  reldiff[gd] /= obserr[gd]

  ; Mask
  ;mask = mask and (abs(coldiff) lt thresh[i] or coldiff gt 90)
  mask = mask and (abs(reldiff) lt 3 or reldiff gt 90)

  hess,str[gd].g-str[gd].i,coldiff[gd],dx=0.05,dy=0.02,yr=[-2,2],/log,xtit='g-i',ytit=color[i],tit=color[i]+' '+strtrim(i+1,2)
  oplot,[-5,5],[0,0]  
  oplot,[-5,5],[1,1]*thresh[i],linestyle=2,co=250
  oplot,[-5,5],[-1,-1]*thresh[i],linestyle=2,co=250

  ; I think I need a better way to set the threshold level

  ;stop

endfor

;gd = where(mask eq 1 and str.g lt 50 and str.i lt 50 and str.u lt 50,ngd)
gd = where(mask eq 1 and str.g lt 50 and str.i lt 50 and str.z lt 50 and $
           str.ndet gt 5 and abs(str.sharp) lt 1 and str.chi lt 3,ngd)
hess,str[gd].g-str[gd].i,str[gd].g,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[25,16],/log


stop

end
