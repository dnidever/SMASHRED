pro synthphot,iso,mass,synthstr

;+
;
; Synthphot
;
; Create a synthetic isochrone stellar population from a single isochrone.
;
; INPUTS:
;  iso         One isochrone.
;  mass        Mass array for each distance and isochrone [Ndist,Niso].
;  /silent     No printing to the screen.
;
; OUTPUTS:
;  synthstr    The synthetic photometry structure.
;
; USAGE:
; IDL>synthphot,iso,mass,synthstr
;
; By D. Nidever  Feb 2018
;    based on cmd_model_smash_synthphot.pro
;-

t0 = systime(1)
deg2rad = !dpi/180d

undefine,synthstr

; Not enough inputs
if n_elements(iso) eq 0 or n_elements(mass) eq 0 then begin
  print,'Syntax - synthphot,iso,mass,synthstr,silent=silent'
  return
endif

; Interpolate to 3x higher sampling of the isochrone
niso = n_elements(iso)
newiso = replicate(iso[0],3*niso-2)
tags = tag_names(iso[0])
ntags = n_tags(iso[0])
for k=0,niso-2 do begin
  ; Linearly interpolate between two points
  iso2 = replicate(iso[k],2)
  for t=0,ntags-1 do begin
    if tags[t] ne 'FEH' and tags[t] ne 'ZMETAL' and tags[t] ne
       'LOGAGE' and tags[t] ne 'PMODE'  and $
       tags[t] ne 'STAGE' then iso2.(t) =
          INTERPOL(iso[k:k+1].(t),[0,1],[1./3.,2./3.])
    endfor
    newiso[k*3:k*3+2] = [iso[k],iso2]
  endfor
  newiso[3*niso-3] = iso[niso-1] ; last one
  oldiso = iso
  iso = newiso & undefine,newiso
  niso = n_elements(iso)

  ; PDF, probability distribution function
  ; int_IMF, which is the integral of the IMF under consideration (as selected in the form, in number of stars,
  ; and normalised to a total mass of 1 Msun) from 0 up to the current M_ini. Differences between 2 values of
  ; int_IMF give the absolute number of stars occupying that isochrone section per unit mass of stellar
  ; population initially born, as expected for the selected IMF.
  ;pdf = [iso[0].int_imf, slope(iso.int_imf) > 1e-9]
  pdf = slope(iso.int_imf) > 1e-9
  pdf = [pdf,first_el(pdf,/last)]
  add_tag,iso,'pdf',0.0,iso
  iso.pdf = pdf
  isoarr[i,j].data = ptr_new(iso)
endfor

stop

; Create the synthetic photometry from the isochrones
;-----------------------------------------------------

str0 = {gc:0,age:0.0,feh:0.0,l:0.0d0,b:0.0d0,dmod:0.0,m_act:0.0,g:99.0,i:99.0,logte:0.0,logg:0.0,gobs:99.0,iobs:99.0}
;str0 = {gc:0,age:0.0,feh:0.0,l:0.0d0,b:0.0d0,dmod:0.0,g:99.0,i:99.0,gobs:99.0,iobs:99.0}
;if n_elements(inpsynthstr) gt 0 then synthstr=inpsynthstr else synthstr=replicate(str0,1e6)
;cnt = 0LL

cnt0 = cnt ; inital count

niso = n_elements(isoarr)
nd = n_elements(dist)
ntotstars = 0LL
nstarsarr = fltarr(nd,niso)

;dt1=0 & dt2=0 & dt3=0 & dt4=0 & dt5=0 & dt6=0

t1 = systime(1) 

; Isochrone loop
For i=0,niso-1 do begin

  iso = *isoarr[i].data
  pdf = iso.pdf
  int_imf = iso.int_imf
  ;pdf = slope(int_imf)
  ;pdf = [pdf,first_el(pdf,/last)]
  ggemax = isoarr[i].ggemax
  glemin = isoarr[i].glemin
  g0 = isoarr[i].g0
  dg = isoarr[i].dg
  ng = isoarr[i].ng
  
  if not keyword_set(silent) then $
    print,'  ',strtrim(i+1,2),' age=',strtrim((10^iso[0].logage)/1e9,2),' Gyr  [Fe/H]=',strtrim(iso[0].feh,2)

  ; Distance loop
  For d=0,nd-1 do begin

;profiler,/reset & profiler & profiler,/system

;t2 = systime(1)
    dmod = dmodarr[d]
    ; Get indices from the GGEMAX and GELEMIN index arrays
    ;  generally the magnitudes decrease with index number
    gmin = minmag-dmod
    gmax = maxmag-dmod
    ;gmin_ind = 0 > floor( (gmin-isoarr[i].g0)/isoarr[i].dg ) < (isoarr[i].ng-1)
    ;gmax_ind = 0 > floor( (gmax-isoarr[i].g0)/isoarr[i].dg ) < (isoarr[i].ng-1)
    ;isohi = isoarr[i].ggemax[gmin_ind]
    ;isolo = isoarr[i].glemin[gmax_ind]
    ;gmin_ind = 0 > floor( (gmin-g0)/dg ) < (ng-1)
    ;gmax_ind = 0 > floor( (gmax-g0)/dg ) < (ng-1)
    ;  using floor/round/long makes this 5x slower
    gmin_ind = 0 > (gmin-g0)/dg < (ng-1)
    gmax_ind = 0 > (gmax-g0)/dg < (ng-1)   
    isohi = ggemax[gmin_ind]
    isolo = glemin[gmax_ind]
    nisoind = isohi-isolo+1
    if isolo lt 0 or isohi lt 0 then goto,BOMB    ; outside range

    ; How many stars do we need for a given TOTAL MASS
    totpdf = iso[isohi].int_imf - iso[isolo].int_imf
    ;totpdf = int_imf[isohi] - int_imf[isolo]
    nstars = round(totpdf*mass[d,i])
    if not keyword_set(silent) then print,'    ',strtrim(dist[d],2),' kpc   ',strtrim(mass[d,i],2),' Msun   Nstars=',strtrim(nstars,2)
    ntotstars += nstars

    ; Add more elements
    if cnt+nstars gt n_elements(synthstr) then begin
      if not keyword_set(silent) then print,'Adding new elements to SYNTHSTR'
      newstr = replicate(str0,1e6)
      synthstr = [temporary(synthstr),newstr]
      undefine,newstr
    endif
    
    ; Create synthetic photometry
    if nstars ge 1 then begin

      ; Create normalized index array (from 0-1)
      indx = findgen(nisoind)/(nisoind-1)
      ; Create cumulative distribution function
      ;cdf = total(pdf[isoind],/cum)
      cdf = int_imf[isolo:isohi]-int_imf[isolo]+pdf[isolo]
      cdf /= max(cdf)

      ; Get the indices in the structure from the CDF
      ;INTERP,cdf,indx,randomu(seed,nstars),newindx
      newindx0 = INTERPOL(indx,cdf,randomu(seed,nstars))

      ; Round to the nearest integer
      ;newindx = min(isoind) > round(newindx0*(nisoind-1)) < max(isoind)
      newindx = isolo > (round(newindx0*(nisoind-1))+isolo) < isohi

      synthstr[cnt:cnt+nstars-1].m_act = iso[newindx].m_act
      synthstr[cnt:cnt+nstars-1].g = iso[newindx].g
      synthstr[cnt:cnt+nstars-1].i = iso[newindx].i
      synthstr[cnt:cnt+nstars-1].logg = iso[newindx].logg
      synthstr[cnt:cnt+nstars-1].logte = iso[newindx].logte
      synthstr[cnt:cnt+nstars-1].age = iso[0].logage
      synthstr[cnt:cnt+nstars-1].feh = iso[0].feh
      ; spread out the stars in distance
      synthstr[cnt:cnt+nstars-1].dmod = dmod + randomu(seed,nstars)*delta_dist[d]
      
      ; Increment the counter
      cnt += nstars
      
     BOMB:
   endif

 Endfor                         ; distance loop
Endfor                          ; isochrone loop

;; No stars
;if cnt eq 0 then begin
;  undefine,synthstr
;  count = 0
;  return,-1
;endif
;synthstr = synthstr[0:cnt-1]    ; trim excess elements
;nsynth = n_elements(synthstr)
;count = nsynth
count = ntotstars
if not keyword_set(silent) then print,strtrim(count,2),' synthetic stars'

; Add observations uncertainties
;gcoef = [-0.984334, 0.182836,  -0.0110896,  0.000223300]
;gerr = poly(synthstr.g+synthstr.dmod,gerr_coef) > 0.017
;icoef = [-3.32317,  0.598515,   -0.0356460,  0.000705355]
;ierr = poly(synthstr.g+synthstr.dmod,ierr_coef) > 0.016
; This block uses 1e-4 sec, very fast
if count gt 0 then begin
  ; spread out the stars in l/b
  synthstr[cnt0:cnt-1].l = glon + (randomu(seed,ntotstars)-0.5)*sqrt(area)/cos(glat/!radeg)
  synthstr[cnt0:cnt-1].b = glat + (randomu(seed,ntotstars)-0.5)*sqrt(area)
  ; observations uncertainties
  gerr_coef = field[0].gerr_coef
  ;gerr = poly(synthstr.g+synthstr.dmod,gerr_coef) > 0.01
  gerr = 10^poly(synthstr[cnt0:cnt-1].g+synthstr[cnt0:cnt-1].dmod,gerr_coef) > 0.001
  ierr_coef = field[0].ierr_coef
  ;ierr = poly(synthstr.g+synthstr.dmod,ierr_coef) > 0.01
  ierr = 10^poly(synthstr[cnt0:cnt-1].g+synthstr[cnt0:cnt-1].dmod,ierr_coef) > 0.001
  synthstr[cnt0:cnt-1].gobs = synthstr[cnt0:cnt-1].g + synthstr[cnt0:cnt-1].dmod + randomn(seed,count)*gerr
  synthstr[cnt0:cnt-1].iobs = synthstr[cnt0:cnt-1].i + synthstr[cnt0:cnt-1].dmod + randomn(seed,count)*ierr
endif

dt = systime(1)-t0
if not keyword_set(silent) then print,'dt = ',dt
;print,'dt = ',dt

;stop

;return,synthstr
;return,-1

end
