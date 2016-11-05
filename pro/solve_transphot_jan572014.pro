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
;22   non-calib, 0.9m data (new)
;26   non-calib
;44   non-calib, LMC-main-body
;56   iz calib, ugr not
;19   non-calib
;55   non-calib, 0.9m data, LMC-main-body
;60   non-calib
;
;So I can use 157, 52, 64, 44, 55, 56 (iz) to derive the photometric
;zeropoints for these three nights.

fields = 'Field'+['157','52','64','44','55','56','22']
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
mfitstr0 = mrdfits('/data/smash/cp/red/photred/stdred/smashred_transphot_eqns.fits',3)
mfitstr_schema = mfitstr0[0]
struct_assign,{dum:''},mfitstr_schema
mfitstr = replicate(mfitstr_schema,3*5)  ; initialize new MJD-level structure

; Initialize the trans structure
fitstr = replicate({mjd:0L,chip:0L,nightnum:0,filter:'',color:'',colband:'',colsign:0,nstars:-1L,$
                   nrejected:-1L,seeing:0.0,zpterm:0.0,zptermerr:0.0,colterm:0.0,coltermerr:0.0,$
                   amterm:0.0,amtermerr:0.0,colamterm:0.0,colamtermerr:0.0,rms:99.99,sig:99.99,$
                   chisq:99.99,medresid:0.0,nbrt:0L,brtrms:99.99,brtsig:99.99,brtchisq:99.99,$
                   badsoln:0,photometric:1,amavgflag:0,namavg:0},nmjd*nfilters*nchips)
cnt = 0LL
undefine,mchstr
observatory,'ctio',obs
for i=0,nmjd-1 do begin
  imjd = mjd[i]
  for j=0,nfilters-1 do begin
    mfitind = where(mfitstr0.mjd ge 56662 and mfitstr0.mjd le 56666 and mfitstr0.filter eq filters[j],nmfitind)
    amterm = mean([mfitstr0[mfitind].amterm])
    amtermsig = mean([mfitstr0[mfitind].amtermsig])
    amavgflag = mfitstr0[mfitind[0]].amavgflag
    namavg = mfitstr0[mfitind[0]].namavg

    mfitstr[i*5+j].mjd = imjd
    ;mfitstr[i*5+j].nightnum  ; filled in before
    ; Convert MJD to YYYYMMDD date
    jd = imjd+2400000.5d0-0.5+obs.tz/24.
    caldat,jd,month,day,year,hour,min,sec
    date = strtrim(year,2)+string(month,format='(I02)')+string(day,format='(I02)')
    mfitstr[i*5+j].date = date
    mfitstr[i*5+j].filter = filters[j]
    ;mfitstr[i*5+j].nstars
    ;mfitstr[i*5+j].brtrms
    ;mfitstr[i*5+j].airmass0
    ;mfitstr[i*5+j].airmass1
    ;mfitstr[i*5+j].amrange
    mfitstr[i*5+j].amterm = amterm
    mfitstr[i*5+j].amtermsig = amtermsig
    ;mfitstr[i*5+j].colamterm
    ;mfitstr[i*5+j].colamtermsig
    mfitstr[i*5+j].badsoln = 0
    mfitstr[i*5+j].photometric = 1
    mfitstr[i*5+j].amavgflag = amavgflag
    mfitstr[i*5+j].namavg = namavg

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

    ; Fill in other night-level information
    gdchip = where(fitstr.mjd eq imjd and fitstr.filter eq filters[j] and fitstr.zpterm ne 0.0,ngdchip)
    mfitstr[i*5+j].nchips = ngdchip
    mfitstr[i*5+j].seeing = median(fitstr[gdchip].seeing)
    mfitstr[i*5+j].rms = median(fitstr[gdchip].sig)
    mfitstr[i*5+j].zpterm = median(fitstr[gdchip].zpterm)
    mfitstr[i*5+j].zptermsig = mad(fitstr[gdchip].zpterm)
    mfitstr[i*5+j].color = fitstr[cnt-1].color
    mfitstr[i*5+j].colband = fitstr[cnt-1].colband
    mfitstr[i*5+j].colsign = fitstr[cnt-1].colsign

  endfor ; filter loop
endfor  ; mjd loop

; chip 2 was already bad at this time.  Just remove chip=2 elements
bd = where(fitstr.chip eq 2,nbd)
REMOVE,bd,fitstr

; Add new elements to FITSTR
fitstr0 = mrdfits('/data/smash/cp/red/photred/stdred/smashred_transphot_eqns.fits',1)
before = where(fitstr0.mjd lt min(mjd),nbefore)
after = where(fitstr0.mjd gt max(mjd),nafter)
newfitstr = [fitstr0[before],fitstr,fitstr0[after]]

; Add new elements to MFITSTR
before = where(mfitstr0.mjd lt min(mjd),nbefore)
after = where(mfitstr0.mjd gt max(mjd),nafter)
newmfitstr = [mfitstr0[before],mfitstr,mfitstr0[after]]

; Renumber NIGHTNUM in fitstr and mfitstr
ui = uniq(newmfitstr.mjd,sort(newmfitstr.mjd))
umjd = newmfitstr[ui].mjd
numjd = n_elements(umjd)
unightnum = lindgen(numjd)+1
for i=0,numjd-1 do begin
  ; Match to FITSTR
  MATCH,newfitstr.mjd,umjd[i],ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then newfitstr[ind1].nightnum = unightnum[i]
  ; Match to MFITSTR
  MATCH,newmfitstr.mjd,umjd[i],ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then newmfitstr[ind1].nightnum = unightnum[i]
endfor

; Load the chip-specific information
chstr0 = mrdfits('/data/smash/cp/red/photred/stdred/smashred_transphot_eqns.fits',2)

; Write the new file
outfile = '/data/smash/cp/red/photred/stdred/smashred_transphot_eqns_new.fits'
print,'Writing new SMASH transformation equation file to: ',outfile
;MWRFITS,newfitstr,outfile,/create
;MWRFITS,chstr0,outfile
;MWRFITS,newmfitstr,outfile

; Add header comments
;head0 = headfits(outfile,exten=0)
;sxaddhist,'SMASH photometric transformation equations for ugriz filters',head0
;sxaddhist,'HDU1: chip+night information (colterm, zeropoint)',head0
;sxaddhist,'HDU2: unique chip-level information (colterm, relative zeropoint)',head0
;sxaddhist,'HDU3: unique night-level information (amterm, nightly zeropoint)',head0
;sxaddpar,head0,'FILTER','giruz'
;sxaddhist,'NOREJECT = 0',head0
;sxaddhist,'NOEXPREJECT = 1',head0
;sxaddhist,'NOSMITH = 1',head0
;sxaddhist,'FITCHIPZPMJDTERM = 0',head0
;sxaddhist,'FIXCHIPZPTERM = 0',head0
;MODFITS,outfile,0,head0,exten=0

MKHDR,head,0
sxaddhist,'SMASH photometric transformation equuations for ugriz filters',head
sxaddhist,'HDU1: chip+night information (colterm, zeropoint)',head
sxaddhist,'HDU2: unique chip-level information (colterm, relative zeropoint)',head
sxaddhist,'HDU3: unique night-level information (amterm, nightly zeropoint)',head
sxaddpar,head,'filter','giruz'
sxaddhist,'NOREJECT = 0',head
sxaddhist,'NOEXPREJECT = 1',head
sxaddhist,'NOSMITH = 1',head
sxaddhist,'FITCHIPZPMJDTERM = 0',head
sxaddhist,'FIXCHIPZPTERM = 0',head
FITS_WRITE,outfile,0,head
MWRFITS,newfitstr,outfile,/silent
MWRFITS,chstr0,outfile,/silent
MWRFITS,newmfitstr,outfile,/silent

;HISTORY SMASH photometric transformation equations for ugriz filters            
;HISTORY HDU1: chip+night information (colterm, zeropoint)                       
;HISTORY HDU2: unique chip-level information (colterm, relative zeropoint)       
;HISTORY HDU3: unique night-level information (amterm, nightly zeropoint)        
;FILTER  = 'giruz   '           /                                                
;HISTORY NOREJECT = 0                                                            
;HISTORY NOEXPREJECT = 1                                                         
;HISTORY NOSMITH = 1                                                             
;HISTORY FITCHIPZPMJDTERM = 0                                                    
;HISTORY FIXCHIPZPTERM = 0   

stop

end
