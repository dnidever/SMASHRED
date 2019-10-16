pro dr2_cpsymlinks

;; Make symlinks for the SMASH DR2 flat files


catdir = '/dl1/users/dnidever/smash/cp/red/photred/'
outdir = '/dl1/users/dnidever/smash/dr2/'

;; Restore the list of images
str = mrdfits('/dl1/users/dnidever/smash/cp/red/photred/observations/smash_dr2_images.fits.gz',1)
str.proctype = strtrim(str.proctype,2)
str.prodtype = strtrim(str.prodtype,2)
str.date_obs = strtrim(str.date_obs,2)
str.fullpath = strtrim(str.fullpath,2)

goto,here

;; Raw
;;-----
graw = where(stregex(str.proctype,'raw',/boolean,/fold_case) eq 1 and str.prodtype eq 'image',ngraw)
for i=0,ngraw-1 do begin
  ;; Get Night
  dateobs = str[graw[i]].date_obs
  dateobs = strmid(dateobs,0,10)+'T'+strmid(dateobs,11)
  mjd = photred_getmjd('','ctio',dateobs=dateobs)  ; night number
  CALDAT,mjd+2400000.5d0-0.5,month,day,year,hour,minute,second
  night = string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')
  ;; Make subdirectory if necessary
  outdir1 = outdir+'raw/'+strtrim(night,2)
  if file_test(outdir1,/directory) eq 0 then file_mkdir,outdir1
  ;; Make symlink
  ;;  use the new naming convention c4d_140603_060636_ori.fits.fz
  base = file_basename(str[graw[i]].fullpath)
  if strmid(base,0,3) ne 'c4d' then begin   ;; construct new name
    origbase = base
    ;; c4d_140106_052819, 2014-01-06T05:28:19
    base = 'c4d_'+strmid(dateobs,2,2)+strmid(dateobs,5,2)+strmid(dateobs,8,2)+'_'+$
                  strmid(dateobs,11,2)+strmid(dateobs,14,2)+strmid(dateobs,17,2)+'_ori.fits.fz'
  endif
  if file_test(outdir1+'/'+base) eq 0 then $
    FILE_LINK,str[graw[i]].fullpath,outdir1+'/'+base
endfor

stop

here:


;; Instcal
;;---------
ginstcal = where(stregex(str.proctype,'InstCal',/boolean,/fold_case) eq 1,nginstcal)
for i=0,nginstcal-1 do begin
  ;; Get Night
  dateobs = str[ginstcal[i]].date_obs
  dateobs = strmid(dateobs,0,10)+'T'+strmid(dateobs,11)
  mjd = photred_getmjd('','ctio',dateobs=dateobs)  ; night number
  CALDAT,mjd+2400000.5d0-0.5,month,day,year,hour,minute,second
  night = string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')
  ;; Make subdirectory if necessary
  outdir1 = outdir+'instcal/'+strtrim(night,2)
  if file_test(outdir1,/directory) eq 0 then file_mkdir,outdir1
  ;; Make symlink
  ;;  use the new naming convention c4d_140603_060636_ori.fits.fz
  base = file_basename(str[ginstcal[i]].fullpath)
  if strmid(base,0,3) ne 'c4d' then begin   ;; construct new name
    origbase = base
    ;; c4d_161030_013526_ooi_z_v1.fits.fz
    ;; ooi - prodtype = image
    ;; ood - prodtype = dqmask
    ;; oow - prodtype = wtmap
    base = 'c4d_'+strmid(dateobs,2,2)+strmid(dateobs,5,2)+strmid(dateobs,8,2)+'_'+$
                  strmid(dateobs,11,2)+strmid(dateobs,14,2)+strmid(dateobs,17,2)+'_oo'
    case str[ginstcal[i]].prodtype of
    'image': type='i'
    'dqmask': type='d'
    'wtmap': type='w'
    else: stop,'prodtype not understood'
    endcase
    base += type+'_'+strmid(str[ginstcal[i]].filter,0,1)+'_v1.fits.fz'
  endif
  if file_test(outdir1+'/'+base) eq 0 then $
    FILE_LINK,str[ginstcal[i]].fullpath,outdir1+'/'+base
endfor

stop

;; Resampled
;;-----------
gresamp = where(stregex(str.proctype,'Resampled',/boolean,/fold_case) eq 1,ngresamp)
for i=0,ngresamp-1 do begin
  ;; Get Night
  dateobs = str[gresamp[i]].date_obs
  dateobs = strmid(dateobs,0,10)+'T'+strmid(dateobs,11)
  mjd = photred_getmjd('','ctio',dateobs=dateobs)  ; night number
  CALDAT,mjd+2400000.5d0-0.5,month,day,year,hour,minute,second
  night = string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')
  ;; Make subdirectory if necessary
  outdir1 = outdir+'resampled/'+strtrim(night,2)
  if file_test(outdir1,/directory) eq 0 then file_mkdir,outdir1
  ;; Make symlink
  ;;  use the new naming convention c4d_140603_060636_ori.fits.fz
  base = file_basename(str[gresamp[i]].fullpath)
  ;if strmid(base,0,3) ne 'c4d' then begin   ;; construct new name
  ;  origbase = base
  ;  ;; c4d_140106_052819, 2014-01-06T05:28:19
  ;  base = 'c4d_'+strmid(dateobs,2,2)+strmid(dateobs,5,2)+strmid(dateobs,8,2)+'_'+$
  ;                strmid(dateobs,11,2)+strmid(dateobs,14,2)+strmid(dateobs,17,2)+'_ori.fits.fz'
  ;endif
  if file_test(outdir1+'/'+base) eq 0 then $
    FILE_LINK,str[gresamp[i]].fullpath,outdir1+'/'+base
endfor

;; Stacked
;;---------
gstack = where(stregex(str.proctype,'stack',/boolean,/fold_case) eq 1,ngstack)
for i=0,ngstack-1 do begin
  ;; Get Night
  dateobs = str[gstack[i]].date_obs
  dateobs = strmid(dateobs,0,10)+'T'+strmid(dateobs,11)
  mjd = photred_getmjd('','ctio',dateobs=dateobs)  ; night number
  CALDAT,mjd+2400000.5d0-0.5,month,day,year,hour,minute,second
  night = string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')
  ;; Make subdirectory if necessary
  outdir1 = outdir+'stacked/'+strtrim(night,2)
  if file_test(outdir1,/directory) eq 0 then file_mkdir,outdir1
  ;; Make symlink
  ;;  use the new naming convention c4d_140603_060636_ori.fits.fz
  base = file_basename(str[gstack[i]].fullpath)
  ;if strmid(base,0,3) ne 'c4d' then begin   ;; construct new name
  ;  origbase = base
  ;  ;; c4d_140106_052819, 2014-01-06T05:28:19
  ;  base = 'c4d_'+strmid(dateobs,2,2)+strmid(dateobs,5,2)+strmid(dateobs,8,2)+'_'+$
  ;                strmid(dateobs,11,2)+strmid(dateobs,14,2)+strmid(dateobs,17,2)+'_ori.fits.fz'
  ;endif
  if file_test(outdir1+'/'+base) eq 0 then $
    FILE_LINK,str[gstack[i]].fullpath,outdir1+'/'+base
endfor


stop

end
