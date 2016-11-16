pro qacheck_wcs,night

; Double-check the WCS fitting

reduxdir = '/data/smash/cp/red/photred/'
if n_elements(night) gt 0 then begin
  dirs = file_search(reduxdir+night,/test_directory,count=ndirs)
endif else begin
  dirs = file_search(reduxdir+'20??????',/test_directory,count=ndirs)
endelse
undefine,str
for i=0,ndirs-1 do begin
  files = file_search(dirs[i]+'/F*/F*-*_??.fits',count=nfiles)
  print,strtrim(i+1,2),'/',strtrim(ndirs,2),' ',dirs[i],' ',strtrim(nfiles,2)
  str1 = replicate({dir:'',file:'',rms:0.0,nmatch:0L,refname:''},nfiles)
  for j=0,nfiles-1 do begin
    head = headfits(files[j])
    ;readline,files[j],gaiahead
    ; Get GAIA RMS and NMATCH
    ;  HISTORY WCSFIT: RMS=0.027 arcsec on  Mon Sep 26 03:09:50 2016
    rmsind = first_el(where(stregex(head,'WCSFIT: RMS',/boolean) eq 1,nrmsind),/last)  ; last one                                                                              
    wcsline = head[rmsind[0]]
    lo = strpos(wcsline,'RMS=')
    tmp = strmid(wcsline,lo+4)
    rms = float( first_el(strsplit(tmp,' ',/extract)) )
    ;  HISTORY WCSFIT: NMATCH=269
    nmatchind = first_el(where(stregex(head,'WCSFIT: NMATCH',/boolean) eq 1,n_nmatchind),/last)
    wcsline = head[nmatchind[0]]
    lo = strpos(wcsline,'NMATCH=')
    tmp = strmid(wcsline,lo+7)
    nmatch = long( first_el(strsplit(tmp,' ',/extract)) )
    ;  HISTORY WCSFIT: Reference catalog=GAIA/GAIA
    refind = first_el(where(stregex(head,'WCSFIT: Reference',/boolean) eq 1,n_nrefind),/last)
    refline = head[refind[0]]
    lo = strpos(refline,'catalog=')
    tmp = strmid(refline,lo+8)
    refname = first_el(strsplit(tmp,' ',/extract))
    str1[j].dir = dirs[i]
    str1[j].file = files[j]
    str1[j].rms = rms
    str1[j].nmatch = nmatch
    str1[j].refname = refname
  endfor
  push,str,str1
endfor

; mwrfits,str,'/data/smash/cp/red/photred/gaia/qacheck_wcs.fits',/create

FINDOUTLIERS:

; Find outliers chips compare to the rest in their exposure
allbase = file_basename(str.file,'.fits')
dum = strsplitter(allbase,'_',/extract)
allexpnum = reform(dum[0,*])
allchip = reform(dum[1,*])
ui = uniq(allexpnum,sort(allexpnum))
expnum = allexpnum[ui]
nexpnum = n_elements(expnum)
bad = intarr(n_elements(str))
for i=0,nexpnum-1 do begin
  MATCH,allexpnum,expnum[i],ind1,ind2,/sort,count=nmatch
  rms = str[ind1].rms
  medrms = median(rms)
  bd1 = where(rms gt medrms*1.5,nbd1)
  bad[ind1[bd1]] = 1
endfor

;mwrfits,bad,'/data/smash/cp/red/photred/gaia/qacheck_wcs_outliers.fits',/create

stop

end
