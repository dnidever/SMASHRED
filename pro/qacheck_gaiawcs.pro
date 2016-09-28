pro qacheck_gaiawcs

; Double-check the GAIA WCS fitting

reduxdir = '/data/smash/cp/red/photred/'
dirs = file_search(reduxdir+'20??????',/test_directory,count=ndirs)
undefine,gaiastr
for i=0,ndirs-1 do begin
  files = file_search(dirs[i]+'/F*/F*-*_??.gaiawcs.head',count=nfiles)
  print,strtrim(i+1,2),'/',strtrim(ndirs,2),' ',dirs[i],' ',strtrim(nfiles,2)
  gaiastr1 = replicate({dir:'',file:'',rms:0.0,nmatch:0L,refname:''},nfiles)
  for j=0,nfiles-1 do begin
    readline,files[i],gaiahead
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
    gaiastr1[j].file = files[i]
    gaiastr1[j].rms = rms
    gaiastr1[j].nmatch = nmatch
    gaiastr1[j].refname = refname
  endfor
  push,gaiastr,gaiastr1
endfor

; mwrfits,gaiastr,'/data/smash/cp/red/photred/gaia/qacheck_gaiawcs.fits',/create

stop

end
