pro qacheck_gaiawcs

; Double-check the GAIA WCS fitting

;gaiastr = mrdfits('/data/smash/cp/red/photred/gaia/qacheck_gaiawcs.fits',1)
;goto,findoutliers

reduxdir = '/data/smash/cp/red/photred/'
dirs = file_search(reduxdir+'20??????',/test_directory,count=ndirs)
undefine,gaiastr
for i=0,ndirs-1 do begin
  files = file_search(dirs[i]+'/F*/F*-*_??.gaiawcs.head',count=nfiles)
  print,strtrim(i+1,2),'/',strtrim(ndirs,2),' ',dirs[i],' ',strtrim(nfiles,2)
  gaiastr1 = replicate({dir:'',file:'',wcsfile:'',rms:0.0,nmatch:0L,refname:''},nfiles)
  for j=0,nfiles-1 do begin
    readline,files[j],gaiahead
    ; Get GAIA RMS and NMATCH
    ;  HISTORY WCSFIT: RMS=0.027 arcsec on  Mon Sep 26 03:09:50 2016
    rmsind = first_el(where(stregex(gaiahead,'WCSFIT: RMS',/boolean) eq 1,nrmsind),/last)  ; last one                                                                              
    wcsline = gaiahead[rmsind[0]]
    lo = strpos(wcsline,'RMS=')
    tmp = strmid(wcsline,lo+4)
    rms = float( first_el(strsplit(tmp,' ',/extract)) )
    ;  HISTORY WCSFIT: NMATCH=269
    nmatchind = first_el(where(stregex(gaiahead,'WCSFIT: NMATCH',/boolean) eq 1,n_nmatchind),/last)
    wcsline = gaiahead[nmatchind[0]]
    lo = strpos(wcsline,'NMATCH=')
    tmp = strmid(wcsline,lo+7)
    nmatch = long( first_el(strsplit(tmp,' ',/extract)) )
    ;  HISTORY WCSFIT: Reference catalog=GAIA/GAIA
    refind = first_el(where(stregex(gaiahead,'WCSFIT: Reference',/boolean) eq 1,n_nrefind),/last)
    refline = gaiahead[refind[0]]
    lo = strpos(refline,'catalog=')
    tmp = strmid(refline,lo+8)
    refname = first_el(strsplit(tmp,' ',/extract))
    gaiastr1[j].dir = dirs[i]
    gaiastr1[j].file = repstr(files[j],'.gaiawcs.head','.fits')
    gaiastr1[j].wcsfile = files[j]
    gaiastr1[j].rms = rms
    gaiastr1[j].nmatch = nmatch
    gaiastr1[j].refname = refname
  endfor
  push,gaiastr,gaiastr1
endfor

; mwrfits,gaiastr,'/data/smash/cp/red/photred/gaia/qacheck_gaiawcs.fits',/create

; 255 >0.1
; 3335 >0.05

FINDOUTLIERS:

; Find outliers chips compare to the rest in their exposure
allbase = file_basename(gaiastr.file,'.gaiawcs.head')
dum = strsplitter(allbase,'_',/extract)
allexpnum = reform(dum[0,*])
allchip = reform(dum[1,*])
ui = uniq(allexpnum,sort(allexpnum))
expnum = allexpnum[ui]
nexpnum = n_elements(expnum)
bad = intarr(n_elements(gaiastr))
for i=0,nexpnum-1 do begin
  MATCH,allexpnum,expnum[i],ind1,ind2,/sort,count=nmatch
  rms = gaiastr[ind1].rms
  medrms = median(rms)
  bd1 = where(rms gt medrms*1.5,nbd1)
  bad[ind1[bd1]] = 1
endfor

;mwrfits,bad,'/data/smash/cp/red/photred/gaia/qacheck_gaiawcs_outliers.fits',/create

stop

end
