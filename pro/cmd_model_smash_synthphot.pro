pro cmd_model_smash_synthphot,field,mass,dist,isoarr,synthstr,cnt,count=count,silent=silent

;+
;
; MDF_MODEL_SYNTHPHOT
;
; Create the synthetic photometry for a single field
;
; INPUTS:
;  field       Structure for a single field with GLON, GLAT, AREA, and
;                GERR_COEF, IERR_COEF, GMIN, and GMAX columns.
;  mass        Mass array for each distance and isochrone [Ndist,Niso].
;  dist        The distance array in kpc.
;  isoarr      Structure array of isochrones to use.
;  /silent     No printing to the screen.
;
; OUTPUTS:
;  synthstr    The synthetic photometry structure.
;  =count      Number of synthetic stars
;
; USAGE:
; IDL>synthstr = mdf_model_synthphot(field,mass,dist,iso)
;
; By D. Nidever  Feb 2015
;-

profiler,/reset & profiler & profiler,/system
  
t0 = systime(1)
deg2rad = !dpi/180d

count = 0

; Not enough inputs
if n_elements(field) eq 0 or n_elements(mass) eq 0 or n_elements(dist) eq 0 or n_elements(isoarr) eq 0 then begin
  print,'Syntax - synthstr =  mdf_model_synthphot(field,mass,dist,isoarr,silent=silent)'
  return
  ;return,-1
endif

glon = field[0].glon
glat = field[0].glat
area = field[0].area
maxmag = field[0].gmax
minmag = field[0].gmin

; Delta distance array
if n_elements(dist) gt 1 then begin
  delta_dist = slope(dist)
  delta_dist = [delta_dist[0],delta_dist]
endif else delta_dist = 0.0
dmodarr = 5*alog10(dist*1e2)  

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
