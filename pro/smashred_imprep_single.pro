pro smashred_imprep_single,file,redo=redo,nodiffmaskflag=nodiffmaskflag

;+
;
; SMASHRED_IMPREP_SINGLE
;
; Get SMASH "calibrated" images ready for PHOTRED
;
; INPUTS:
;  file   The CP ooi (calibrated image) file. c4d_140528_023354_ooi_g_v1.fits.fz
;  /nodiffmaskflag   Turn off the CP difference image masking flag
;                      which sometimes has problems.  Set by default.
;  /redo  Redo/overwrite previous one.
; 
; OUTPUTS:
;  This gets CP-reduced images ready for PHOTRED/DAOPHOT
; 
;+ By D.Nidever  2015

if n_elements(file) eq 0 then begin
  print,'Syntax - smashred_imprep_single,file,redo=redo'
  return
endif
gd = where(stregex(file,'_ooi_',/boolean) eq 1,ngd)
if ngd eq 0 then begin
  print,file,' is NOT an ooi file'
  return
endif

; Get information on all the images
fstr = {file:'',base:'',object:'',ra:'',dec:'',exptime:0.0,filter:'',expnum:'',$
                  shname:'',fname:'',newfile:''}
head = headfits(file,exten=0)
base = strmid(file_basename(file),0,17)

; Get header information
orig_filename = sxpar(head,'DTACQNAM')
; /data_local/images/DTS/2013B-0440/DECam_00317720.fits.fz
expnum = strmid(file_basename(orig_filename),6,8)  

fstr.file = file
fstr.base = base
fstr.object = strtrim(sxpar(head,'object'),2)
fstr.ra = sxpar(head,'ra')
fstr.dec = sxpar(head,'dec')
fstr.exptime = sxpar(head,'exptime')
fstr.filter = strmid(strtrim(sxpar(head,'filter'),2),0,1)
fstr.expnum = expnum

print,file_basename(file),'  ',fstr.object,'  ',fstr.filter,'  ',$
      strtrim(fstr.exptime,2),'  ',fstr.ra,'  ',fstr.dec



; Make the PHOTRED-ready FITS files for the science exposures
;-------------------------------------------------------------
imfile = fstr.file
base = strmid(file_basename(imfile),0,17)
maskfile = repstr(imfile,'ooi','ood')
wtfile = repstr(imfile,'ooi','oow')
; ood  calibrated dqmask
; ooi  calibrated image  (nextend=60)
; oow  calibrated weight

if file_test(maskfile) eq 0 then begin
  print,maskfile,' NOT FOUND'
  return
endif

; Get header information
head0 = headfits(imfile,exten=0)
orig_filename = sxpar(head0,'DTACQNAM')
; /data_local/images/DTS/2013B-0440/DECam_00317720.fits.fz
expnum = strmid(file_basename(orig_filename),6,8)  
filter = sxpar(head0,'filter')
exptime = sxpar(head0,'exptime')
;print,strtrim(i+1,2),'/',strtrim(nfstr,2),' ',base,' ',exptime,' ',filter

; Has this one been done already
if file_test(expnum+'_01.fits') and not keyword_set(redo) then begin
  print,expnum,' already done and not /redo set.  Skipping'
  return
endif

; Uncompress
print,'Uncompressing'
if file_test(file_basename(imfile,'.fz')) eq 0 then spawn,'funpack '+imfile
if file_test(file_basename(maskfile,'.fz')) eq 0 then spawn,'funpack '+maskfile
;if file_test(file_basename(wtfile,'.fz')) eq 0 then spawn,'funpack '+wtfile
imfile = file_basename(imfile,'.fz')
maskfile = file_basename(maskfile,'.fz')
;wtfile = file_basename(wtfile,'.fz')

; How many extensions are there
fits_open,imfile,fcb
next = fcb.nextend
fits_close,imfile  
print,strtrim(next,2),' extensions'

; Chip loop
for j=1,next-1 do begin
  fits_read,imfile,fim,fhead,exten=j
  fits_read,maskfile,mim,mhead,exten=j
  ;fits_read,wtfile,wim,whead,exten=j
  ccdnum = sxpar(fhead,'CCDNUM')

  ; Need to add gain, rdnoise, saturation

  ; Fix the background of the two amps
  ;med1b = median(fim[0:1023,0:1998])
  ;med1t = median(fim[0:1023,1999:*])
  ;med2b = median(fim[1024:*,0:1998])
  ;med2t = median(fim[1024:*,1999:*])
  ;med = median(fim)
  ;print,med
  ;print,[med1t,med2t]-med
  ;print,[med1b,med2b]-med
  
  med1 = median(fim,dim=1)
  ;sig = mad(med1[0:1900])
  med1slp = slope(med1)
  
  ;med1 = median(fim[*,1800:1998])
  ;med2 = median(fim[*,1998:2200])
  ;sig = mad(fim[*,1800:1990])
  ;print,med1,med2,sig

  ;sig = mad(med1slp)
  ;;bd = where(abs(med1slp[1950:2050]) gt 4*sig,nbd)
  ;bd = where(abs(med1slp[1997]) gt 4*sig,nbd)
  ;print,nbd

  newim = fim
  ;if nbd gt 0 then begin
  ;  med1 = median(fim[*,1995:1997])
  ;  med2 = median(fim[*,1998:2000])
  ;  diff = med2-med1
  ;  newim[*,0:1997] += diff*0.5
  ;  newim[*,1998:*] -= diff*0.5
  ;endif

  ; Check for differences in amp background levels
  med1 = median(newim[800:1023,*])
  med2 = median(newim[1024:1200,*])
  err1 = mad(newim[800:1023,*])/sqrt(n_elements(newim[800:1023,*]))
  err2 = mad(newim[1024:1200,*])/sqrt(n_elements(newim[1024:1200,*]))
  err = sqrt(err1^2 + err2^2)
  ;print,med1,med2

  ;if abs(med1-med2) gt 5*err then begin
  ;  newim[0:1023,*] -= (med1-med2)*0.5
  ;  newim[1024:*,*] += (med1-med2)*0.5
  ;  print,'amp offset ',med1,med2,abs(med1-med2)
  ;endif


  ; Set bad pixels to saturation value
  ; --DESDM bit masks (from Gruendl):
  ; BADPIX_BPM 1          /* set in bpm (hot/dead pixel/column)        */
  ; BADPIX_SATURATE 2     /* saturated pixel                           */
  ; BADPIX_INTERP 4
  ;     /* interpolated pixel                        */
  ; BADPIX_LOW     8      /* too little signal- i.e. poor read         */
  ; BADPIX_CRAY   16      /* cosmic ray pixel                          */
  ; BADPIX_STAR   32      /* bright star pixel                         */
  ; BADPIX_TRAIL  64      /* bleed trail pixel                         */
  ; BADPIX_EDGEBLEED 128  /* edge bleed pixel                          */
  ; BADPIX_SSXTALK 256    /* pixel potentially effected by xtalk from super-saturated source */
  ; BADPIX_EDGE   512     /* pixel flagged to exclude CCD glowing edges */
  ; BADPIX_STREAK 1024    /* pixel associated with satellite (airplane/meteor) streak     */
  ; BADPIX_FIX    2048    /* a bad pixel that was fixed                */
  ; --CP bit masks, Pre-V3.5.0 (PLVER)
  ; Bit   DQ Type  PROCTYPE
  ; 1  detector bad pixel          InstCal
  ; 1  detector bad pixel/no data  Resampled
  ; 1  No data                     Stacked
  ; 2  saturated                   InstCal/Resampled
  ; 4  interpolated                InstCal/Resampled
  ; 16  single exposure cosmic ray InstCal/Resampled
  ; 64  bleed trail                InstCal/Resampled
  ; 128  multi-exposure transient  InstCal/Resampled
  ; --CP bit masks, V3.5.0 on (after ~10/28/2014), integer masks
  ;  1 = bad (in static bad pixel mask)
  ;  2 = no value (for stacks)
  ;  3 = saturated
  ;  4 = bleed mask
  ;  5 = cosmic ray
  ;  6 = low weight
  ;  7 = diff detect
  ; You can't have combinations but the precedence as in the order
  ; of the list (which is also the order in which the processing
  ; discovers them).  So a pixel marked as "bad" (1) won't ever be
  ; flagged as "diff detect" (7) later on in the processing.
  ;
  ; "Turn off" the "difference image masking", clear the 8th bit
  ; 128 for Pre-V3.5.0 images and set 7 values to zero for V3.5.0 or later.
  if n_elements(nodiffmaskflag) eq 0 then nodiffmaskflag = 1  ; set by default
  if keyword_set(nodiffmaskflag) then begin
    print,'Turning off the CP difference image masking flags'
    plver = sxpar(fhead,'plver',count=nplver)  ; DESDM doesn't have this
    plver = strtrim(plver,2)
    if nplver gt 0 then begin  ; CP data
      ; V3.5.0 and on, Integer masks
      versnum = long(strsplit(strmid(plver,1),'.',/extract))
      if versnum[0] gt 3 or (versnum[0] eq 3 and versnum[1] ge 5) then begin
        bdpix = where(mim eq 7,nbdpix)
        if nbdpix gt 0 then mim[bdpix]=0

      ; Pre-V3.5.0, Bitmasks
      endif else begin
        bdpix = where( (mim and 2^7) eq 2^7,nbdpix)
        if nbdpix gt 0 then mim[bdpix]-=128   ; clear 128
      endelse
      print,strtrim(nbdpix,2),' pixels cleared of difference image mask flag'
    endif
  endif

  bdpix = where(mim gt 0.0,nbdpix)
  if nbdpix gt 0 then newim[bdpix]=6e4

  ; add gain, rdnoise, saturation
  newhead = fhead
  if strmid(newhead[0],0,5) eq 'XTENS' then newhead[0]='SIMPLE  =                    T / Fits standard'

  ;gain = (arr[ccd-1].gaina+arr[ccd-1].gainb)*0.5
  ;rdnoise = (arr[ccd-1].rdnoisea+arr[ccd-1].rdnoiseb)*0.5
  ;gain = sxpar(fhead,'ARAWGAIN')
  gainA = sxpar(fhead,'GAINA')
  gainB = sxpar(fhead,'GAINB')
  gain = (gainA+gainB)*0.5
  rdnoiseA = sxpar(fhead,'RDNOISEA')
  rdnoiseB = sxpar(fhead,'RDNOISEB')
  rdnoise = (rdnoiseA+rdnoiseB)*0.5
  sxaddpar,newhead,'GAIN',gain
  sxaddpar,newhead,'RDNOISE',rdnoise

  ; Remove second SIMPLE
  ;indsimple = where(stregex(newhead,'^SIMPLE',/boolean) eq 1,nindsimple)
  ;if nindsimple gt 1 then remove,indsimple[1:*],newhead

  ; REMOVE DUPLICATE KEYWORDS!!  They cause lots of annoying errors
  ; EXTVER, CHECKSUM, DATASUM
  bd = where(strmid(newhead,0,6) eq 'EXTVER',nbd)
  if nbd gt 1 then remove,bd[1:*],newhead
  bd = where(strmid(newhead,0,8) eq 'CHECKSUM',nbd)
  if nbd gt 1 then remove,bd[1:*],newhead
  bd = where(strmid(newhead,0,7) eq 'DATASUM',nbd)
  if nbd gt 1 then remove,bd[1:*],newhead

  ; Add "COMMENT " before "BEGIN EXTENSION HEADER ---", it causes problems in daophot
  bd = where(strmid(newhead,0,5) eq 'BEGIN',nbd)
  if nbd gt 0 then newhead[bd]='COMMENT '+newhead[bd]

  ; Write new image
  newfile = expnum+'_'+string(ccdnum,format='(I02)')+'.fits'
  ;;fits_write,newfile,newim,newhead
  print,'Writing ',newfile,' ',gain,' ',rdnoise
  mwrfits,newim,newfile,newhead,/create

  ;stop
endfor

; Delete the uncompressed files
if file_test(imfile) and file_test(imfile+'.fz') then file_delete,imfile,/verbose
if file_test(maskfile) and file_test(maskfile+'.fz') then file_delete,maskfile,/verbose
;if file_test(wtfile) and file_test(wtfile+'.fz') then file_delete,wtfile,/verbose


;stop

end
