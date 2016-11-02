pro solve_transphot_jan572014

; Figure out the photometric zeropoints for the Jan 507, 2014
; nights which were photometric but no standard star fields
; were observed.  Use fields that were also observed and
; calibration from other nights.

;    8   56663   20140105        1       no clouds  157, 52, 64
;    9   56664   20140106        1       no clouds  157, 20, 22, 26, 44, 56, car1-4
;   10   56665   20140107        1       no clouds  157, 19, 26, 55, 60
;no calibration data for those nights
;157  griz calib, u not
;52   non-calib, LMC-main-body
;64   non-calib, 0.9m data
;20   non-calib
;22   non-calib
;26   non-calib
;44   non-calib, LMC-main-body
;56   iz calib, ugr not
;19   non-calib
;55   non-calib, 0.9m data, LMC-main-body
;60   non-calib
;
;So I can use 157, 52, 64, 44, 55, 56 (iz) to derive the photometric
;zeropoints for these three nights.

fields = 'Field'+['157','52','64','44','55','56']
nfields = n_elements(fields)
dir = '/data/smash/cp/red/photred/catalogs/final/v3/'
; Load the calibrated catalogs
undefine,allobj,allchstr
for i=0,nfields-1 do begin
  print,strtrim(i+1,2),' ',fields[i]
  ;file = file_search(dir+fields[i]+'_combined_allobj_bright.fits',count=nfile)
  chfile = file_search(dir+fields[i]+'_combined_chips.fits.gz',count=nchfile)
  if nchfile gt 0 then begin
    ;obj = mrdfits(file,1)
    ;push,allobj,obj
    chstr = mrdfits(chfile,1)

    ; remove alftiletype
    if tag_exist(chstr,'alftiletype') then begin
      chstr0 = chstr
      dum = allchstr[0]
      struct_assign,{dum:''},dum
      chstr = replicate(dum,n_elements(chstr0))
      struct_assign,chstr0,chstr
    endif

    push,allchstr,chstr
  endif
endfor

; use the ZPTERM information to derive chip-level
; ZPTERM for each night


;  8   56663   20140105        1       no clouds
;  9   56664   20140106        1       no clouds
; 10   56665   20140107        1       no clouds

mjd = ['56663','56664','56665']
nmjd = n_elements(mjd)
filters = ['u','g','r','i','z']
nfilters = n_elements(filters)
ui = uniq(allchstr.chip,sort(allchstr.chip))
uchips = allchstr[ui].chip
;match,uchips,2,ind1,ind2  ; don't need chip 2
;remove,ind1,uchips
nchips = n_elements(uchips)

; Restore the night-level information to get the extinction terms
mfitstr = mrdfits('/data/smash/cp/red/photred/stdred/smashred_transphot_eqns.fits',3)

; Initialize the trans structure
fitstr = replicate({mjd:0L,chip:0L,nightnum:0,filter:'',color:'',colband:'',colsign:0,nstars:0L,$
                   nrejected:0L,seeing:0.0,zpterm:0.0,zptermerr:0.0,colterm:0.0,coltermerr:0.0,$
                   amterm:0.0,amtermerr:0.0,colamterm:0.0,colamtermerr:0.0,rms:99.99,sig:99.99,$
                   chisq:99.99,medresid:0.0,nbrt:0L,brtrms:99.99,brtsig:99.99,brtchisq:99.99,$
                   badsoln:0,photometric:1,amavgflag:0,namavg:0},nmjd*nfilters*nchips)
cnt = 0LL
undefine,mchstr
for i=0,nmjd-1 do begin
  imjd = mjd[i]
  for j=0,nfilters-1 do begin
    mfitind = where(mfitstr.mjd ge 56662 and mfitstr.mjd le 56666 and mfitstr.filter eq filters[j],nmfitind)
    amterm = mean([mfitstr[mfitind].amterm])
    amtermsig = mean([mfitstr[mfitind].amtermsig])
    amavgflag = mfitstr[mfitind[0]].amavgflag
    namavg = mfitstr[mfitind[0]].namavg

    for k=0,nchips-1 do begin
      fitstr[cnt].mjd = imjd
      fitstr[cnt].chip = uchips[k]
      ;fitstr[cnt].nightnum = 
      fitstr[cnt].filter = filters[j]
      ;fitstr[cnt].color = 
      ;fitstr[cnt].colband 
      ;fitstr[cnt].colsign = 
      ;fitstr[cnt].nstars = 
      ind = where(allchstr.mjd eq imjd and allchstr.filter eq filters[j] and $
                  allchstr.chip eq uchips[k] and allchstr.zpcalibflag gt 0 and allchstr.zpcalibflag lt 4,nind)
      if nind gt 0 then begin
        push,mchstr,allchstr[ind]
        fitstr[cnt].seeing = median([allchstr[ind].fwhm*allchstr[ind].pixscale])
        fitstr[cnt].colband = allchstr[ind[0]].colband
        fitstr[cnt].colsign = allchstr[ind[0]].colsign
        if fitstr[cnt].colsign eq 1 then begin
          fitstr[cnt].color = fitstr[cnt].filter+'-'+allchstr[ind[0]].colband
        endif else begin
          fitstr[cnt].color = allchstr[ind[0]].colband+'-'+fitstr[cnt].filter
        endelse
        ; Stuff in the airterm information
        fitstr[cnt].amterm = amterm
        fitstr[cnt].amtermerr = amtermsig
        fitstr[cnt].amavgflag = amavgflag
        fitstr[cnt].namavg = namavg
        ; Stuff in the colterm information
        fitstr[cnt].colterm = allchstr[ind[0]].colterm
        fitstr[cnt].coltermerr = allchstr[ind[0]].coltermsig

        ; Need to remove the extinction portion of the zpterm
        med = median([allchstr[ind].zpterm-(allchstr[ind].airmass*amterm)])
        ; Has to be this sign, other way the answer is WAY OFF!!
        ;med = median([allchstr[ind].zpterm])
        sig = mad([allchstr[ind].zpterm])
        fitstr[cnt].zpterm = med
        fitstr[cnt].sig = sig
        fitstr[cnt].zptermerr = sig/sqrt(nind)
        print,imjd,'  ',filters[j],'  ',strtrim(uchips[k],2),'  ',strtrim(nind,2),'  ',med,sig
      endif else begin
        print,imjd,'  ',filters[j],'  ',strtrim(uchips[k],2),' none'
      endelse
      cnt++
    endfor ; chips loop
  endfor ; filter loop
endfor  ; mjd loop

; use mean of neighboring chips for chip 2

; No U-BAND for second night!!!
; We could use the calibration for night 1+3 to calibrate Field157
; and then use that to calibrate u-band for night 2.
; OR use new 0.9m calibration from Field22 (new, sdssphot)

stop

end
