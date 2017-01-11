;+
;
; SMASHRED_2CDCUTS
;
; Perform star/galaxy separation with the stellar locus in multiple
; color-color planes.
;
; INPUTS:
;  str      A structure of SMASH photometry.
;  =gmax    The maximum g-value to use when defining the stellar
;             locus.  gmas=23.0 by default.
;  /pl      Make some plots.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  ind      The indices for elements that passed the cuts.
;  tstr     The stellar locus structure.
;
; USAGE:
;  IDL>smashred_2cdcuts,str,ind
;
; By D.Nidever  Jan 2017
;-

pro smashred_2cdcuts,str,ind,tstr,gmax=gmax,pl=pl,stp=stp

; Use color-color cuts for star-galaxy separation

; Use b-spline fits to low-err sources and
; then use simple differences in Y to remove "bad" sources
; maybe only use the cuts for regions where the slope is "low".

undefine,ind
undefine,tstr
  
; Not enough inputs
nstr = n_elements(str)
if nstr eq 0 then begin
  print,'Syntax - smashred_2cdcuts,str,ind,tstr,gmax=gmax,pl=pl,stp=stp'
  return
endif

; de-extincting
;magext = [4.239, 3.303, 2.285, 1.698, 1.263]
;u0=str.u-magext[0]*str.ebv
;g0=str.g-magext[1]*str.ebv
;r0=str.r-magext[2]*str.ebv
;i0=str.i-magext[3]*str.ebv 
;z0=str.z-magext[4]*str.ebv

; this does all of the bspline fitting
if n_elements(gmax) eq 0 then gmax=23.0
MAKE_STELLAR_LOCUS,str,tstr,gmax=gmax

; u-g, u-r, u-i vs. g-i  good separation
; g-z, r-z, i-z vs. g-i is okay but still quite a bit of overlap with galaxies
; u-z, g-r, r-i vs. g-i not very useful

color = ['g-z','r-z','i-z']
thresh = [0.13,0.13,0.13]
athresh = [0.2,0.2,0.2]
rthresh = [2.5,2.5,2.5]
;color = ['u-g','u-r','u-i','g-z','r-z','i-z']
;thresh = [0.3,0.3,0.3,0.13,0.13,0.13]
;color = ['u-g','u-r','u-i']
;thresh = [0.3,0.3,0.3]
ncolor = n_elements(color)

; looking at other color-color plots in ~/decam/plots/FieldB/FieldB_2cdcuts.pdf
; g-z vs. g-r, okay
; r-z vs. g-r, okay
; i-z vs. g-r, pretty good
; g-z vs. g-i, okay
; r-z vs. g-i, pretty good
; i-z vs. g-i, quite good
; r-i vs. g-z, fine
; r-z vs. g-z, okay
; i-z vs. g-z, okay
; r-z vs. r-z, fine
; i-z vs. r-i, okay
; i-z vs. r-z, okay

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

  ; Interpolate the model color for each source  
  gd = where(str.g lt 50 and str.i lt 50 and str.(magind1) lt 50 and str.(magind2) lt 50 and $
             str.g-str.i ge min(tstr.gibin) and str.g-str.i le max(tstr.gibin),ngd)
  INTERP,tstr.gibin,tcolbin,str[gd].g-str[gd].i,tcol

  coldiff = fltarr(nstr)+99.99
  coldiff[gd] = col[gd] - tcol

  ; Divide by observation errors
  errind1 = where(stags eq strupcase(band1)+'ERR',nerrind1)
  errind2 = where(stags eq strupcase(band2)+'ERR',nerrind2)
  obserr = sqrt( str.(errind1)^2 + str.(errind2)^2 ) > 0.03  ; lower threshold
  reldiff = coldiff
  reldiff[gd] /= obserr[gd]

  ; Mask
  ;mask = mask and (abs(coldiff) lt thresh[i] or coldiff gt 90)
  ;mask = mask and (abs(reldiff) lt 3 or reldiff gt 90)
  ;mask = mask and ((abs(reldiff) lt 3 and abs(coldiff) lt 0.5) or reldiff gt 90)
  ;mask = mask and (abs(reldiff) lt rthresh[i] or reldiff gt 90)
  mask = mask and ((abs(reldiff) lt rthresh[i] and abs(coldiff) lt athresh[i]) or reldiff gt 90)

  ;; Plotting
  if keyword_set(pl) then begin
    ;hess,str[gd].g-str[gd].i,coldiff[gd],dx=0.05,dy=0.03,yr=[-2,2],/log,xtit='g-i',ytit='Delta '+color[i],$
    ;     tit=color[i]+' '+strtrim(i+1,2),position=[0,0.5,1,1]
    ;oplot,[-5,5],[0,0]  
    ;oplot,[-5,5],[1,1]*thresh[i],linestyle=2,co=250
    ;oplot,[-5,5],[-1,-1]*thresh[i],linestyle=2,co=250
    hess,str[gd].g-str[gd].i,reldiff[gd],dx=0.05,dy=0.1,yr=[-10,10],/log,xtit='g-i',ytit='Relative delta '+color[i],$
         tit=color[i]+' '+strtrim(i+1,2),position=[0,0.5,1,1]
    oplot,[-5,5],[0,0]  
    oplot,[-5,5],[1,1]*rthresh[i],linestyle=2,co=250
    oplot,[-5,5],[-1,-1]*rthresh[i],linestyle=2,co=250

    ;hess,str[gd].g,coldiff[gd],str[gd].g-str[gd].i,dx=0.1,dy=0.03,xr=[15,26],yr=[-2,2],/log,xtit='g',$
    ;     ytit='Delta '+color[i],tit=color[i]+' '+strtrim(i+1,2),position=[0,0,1,0.5],/noerase
    ;oplot,[10,30],[0,0]  
    ;oplot,[10,30],[1,1]*thresh[i],linestyle=2,co=250
    ;oplot,[10,30],[-1,-1]*thresh[i],linestyle=2,co=250
    hess,str[gd].g,reldiff[gd],str[gd].g-str[gd].i,dx=0.1,dy=0.1,xr=[15,26],yr=[-10,10],/log,xtit='g',$
         ytit='Relative delta '+color[i],tit=color[i]+' '+strtrim(i+1,2),position=[0,0,1,0.5],/noerase
    oplot,[10,30],[0,0]  
    oplot,[10,30],[1,1]*rthresh[i],linestyle=2,co=250
    oplot,[10,30],[-1,-1]*rthresh[i],linestyle=2,co=250

    stop
  endif
  
  ; I think I need a better way to set the threshold level

  ;stop

endfor

ind = where(mask eq 1,ningd)

if keyword_set(pl) then begin
  ;gd = where(mask eq 1 and str.g lt 50 and str.i lt 50 and str.u lt 50,ngd)
  gd = where(mask eq 1 and str.g lt 50 and str.i lt 50 and str.z lt 50 and $
             str.ndet gt 5 and abs(str.sharp) lt 1 and str.chi lt 3,ngd)
  hess,str[gd].g-str[gd].i,str[gd].g,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[25,16],/log
endif
  
; at faint mags the stars depart from the "bright" stellar locus in
; u-X bands.  Don't see this in the other colors.

; The cut sources
;bd=where(mask eq 0 and str.g lt 50 and str.i lt 50 and str.z lt 50)
;hess,str[bd].g-str[bd].i,str[bd].g,dx=0.02,dy=0.05,xr=[-1,3.5],yr=[25,16],/log

if keyword_set(stp) then stop

end
