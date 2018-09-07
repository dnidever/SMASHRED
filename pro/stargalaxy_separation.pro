pro stargalaxy_separation,obj0,ind,astobj=astobj0,astind=astind

;; Make "star" catalogs.

undefine,ind

nobj0 = n_elements(obj0)
if nobj0 eq 0 then begin
  print,'Syntax - stargalaxy_separation,obj,ind'
  return
endif

tags = tag_names(obj0)

; Maybe only look at depthflag=2 or 3
  
;; Deredden
;obj = obj0
;print,'NO DEREDDENING'
;print,'  Dereddening'
;magext = [4.239, 3.303, 2.285, 1.698, 1.263]
;obj.u -= magext[0]*obj.ebv
;obj.g -= magext[1]*obj.ebv
;obj.r -= magext[2]*obj.ebv
;obj.i -= magext[3]*obj.ebv
;obj.z -= magext[4]*obj.ebv

;; Select the stars
depththresh = 1
if max(obj0.depthflag) eq 1 then depththresh=0  ; in case only short exposures
;gdstars = where(abs(obj0.sharp) lt 1 and obj0.chi lt 2 and obj0.prob gt 0.2 and $
;; stars1 had the PROB cut but it was removing lots of bright stars,
;; probably because of the brighter-fatter effect.
gdstars = where(abs(obj0.sharp) lt 1 and obj0.chi lt 2 and $
                obj0.ndet gt 5 and obj0.depthflag gt depththresh,ngdstars)
obj = obj0[gdstars]

; ASTs
nastobj = n_elements(astobj0)
if nastobj gt 0 then begin
  ;ndet = astobj0.ndetu + astobj0.ndetg + astobj0.ndetr + astobj0.ndeti + astobj0.ndetz
  gdast = where(abs(astobj0.sharp) lt 1 and astobj0.chi lt 2 and astobj0.prob gt 0.2 and $
                astobj0.ndet gt 5,ngdast)
  ; ASTs only run on deep exposures, so no DEPTHFLAG column
  ;                astobj0.ndet gt 5 and astobj0.depthflag gt 1,ngdast)
  astobj = astobj0[gdast]
endif

;; Morphology cuts
SMASHRED_MORPHCUTS,obj,ind1,nsig=3,astobj=astobj,astind=astind1
obj1 = obj[ind1]
if nastobj gt 0 then astobj1=astobj[astind1]
  
;; Color-color cuts
SMASHRED_2CDCUTS,obj1,ind2,aststr=astobj1,astind=astind2
obj2 = obj1[ind2]
if nastobj gt 0 then astobj2=astobj1[astind2]
; this cuts out some stars with g>24 on Field56

; Final indices
ind = gdstars[ind1[ind2]]
if nastobj gt 0 then astind = gdast[astind1[astind2]]

;stop

end
