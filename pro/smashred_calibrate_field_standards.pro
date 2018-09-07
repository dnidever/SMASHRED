;+
;
; SMASHRED_CALIBRATE_FIELD_STANDARDS
;
; This program calibrates all of the photometry for a single SMASH
; standard star field using transformation equations and ubercal techniques.
;
; INPUTS:
;  field      The field name, e.g. "Field24"
;  =version   The version of the final catalogs.  This is only used
;               if OUTPUTDIR is not input.
;  =sumfiles  Input list of summary files to use.  Otherwise smashred_getredinfo
;               uses all summary files for FIELD in dir+'/20??????/*summary.fits'
;  =transfile The file with the photometric transformation equations.
;  =reduxdir  The base directory for the PHOTRED reductions.  The
;               default is "/data/smash/cp/red/photred/"
;  =outputdir The output directory for the catalogs.
;  /usegaia   Use GAIA photometry to set the photometric zeropoint.
;               This is the default
;  /compress  Gzip compress the output FITS files.  This is the default.
;  /redo      Redo a field that was already calibrated.
;  /silent    Don't print anything to the screen.
;
; OUTPUTS:
;  Five binary FITS files are created in the output directory:
;   FIELD_combined_exposures.fits - information on each unique exposure
;   FIELD_combined_chips.fits - information on each unique chip
;   FIELD_combined_allsrc.fits - information on each unique source detection
;   FIELD_combined_allobj.fits - information on each unique object including
;                                     average magnitudes.
;   FIELD_combined_allobj_bright.fits - bright stars in allobj used for
;                                       cross-matching between fields.
;   =error    The error message, if one occurred.
;
; USAGE:
;  IDL>smashred_calibrate_field_standards,'Field24'
;
; By D.Nidever March 2016
;-

pro smashred_calibrate_field_standards,field,all=all
;  ,version=version,sumfiles=sumfiles,transfile=transfile,reduxdir=reduxdir,outputdir=outputdir,$
;                             usegaia=usegaia,redo=redo,compress=compress,silent=silent,error=error

undefine,error
t0 = systime(1)

;field = 'SDSSJ0000-0000_sdss'

; Not enough inputs
if n_elements(field) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_calibrate_field_standards,field,version=version,sumfiles=sumfiles,reduxdir=reduxdir,outputdir=outputdir,redo=redo,'
  print,'                                  usegaia=usegaia,compress=compress,error=error'
  return
endif

; Defaults
;if n_elements(reduxdir) eq 0 then reduxdir=SMASHRED_ROOTDIR()+'cp/red/photred/'
;;if n_elements(reduxdir) eq 0 then reduxdir='/data/smash/cp/red/photred/'
;if file_test(reduxdir,/directory) eq 0 then begin
;  error = reduxdir+' NOT FOUND'
;  if not keyword_set(silent) then print,error
;  return
;endif
;if n_elements(outputdir) eq 0 then begin
;  outputdir = reduxdir+'catalogs/final/'
;  if n_elements(version) gt 0 then outputdir+=version+'/'
;endif

transfile = '/Users/nidever/projects/SMASHRED/data/smashred_transphot_eqns.fits'
reduxdir = '/Users/nidever/smash/reduction/stdred/'
outputdir = '/Users/nidever/smash/reduction/stdred/standards/'

;if file_test(outputdir,/directory) eq 0 then begin
;  if not keyword_set(silent) then print,outputdir+' does NOT exist.  Creating it.'
;  FILE_MKDIR,outputdir
;endif
;if n_elements(transfile) eq 0 then transfile=SMASHRED_ROOTDIR()+'cp/red/photred/stdred/smashred_transphot_eqns.fits'
;;if n_elements(transfile) eq 0 then transfile='/data/smash/cp/red/photred/stdred/smashred_transphot_eqns.fits'
;; Temporary directory
;tmpdir = outputdir+'/tmp/'
;if file_test(tmpdir,/directory) eq 0 then FILE_MKDIR,tmpdir
;; Compression
;if n_elements(compress) eq 0 then compress=1
;; GAIA photometry
;if n_elements(usegaia) eq 0 then usegaia=0

;; Load standards ALLSRC file
;;  created by hand.  see notes for 5/11/2018
if n_elements(all) eq 0 then begin
  print,'Loading standards ALLSRC file'
  all = mrdfits(outputdir+'stdred_allsrc.fits.gz',1)
  all.id = strtrim(all.id,2)
  all.stdfield = strtrim(all.stdfield,2)
  all.frame = strtrim(all.frame,2)
endif

; Check if the final files already exists
outfile = outputdir+field+'_combined'
;testfiles = outfile+['_exposures','_chips','_allsrc','_allobj']+'.fits'
testfiles = outfile+['_allsrc','_allobj']+'.fits'
ntestfiles = n_elements(testfiles)
if (total(file_test(testfiles)) eq ntestfiles or total(file_test(testfiles+'.gz')) eq ntestfiles) and $
  not keyword_set(redo) then begin
  print,'Final output files already exist for ',field,' and /redo NOT set.'
  return
endif

;; 39 separate fields, only use ones that start with SA or SDSS

;; Create INFO structure
info = {field:field}

;; Get source measurements for this field
gfield = where(all.stdfield eq field,ngfield)
allsrc = all[gfield]

;; THE STANDARD STAR PHOTOMETRY WAS ALREADY CORRECTED FOR EXPTIME
;; IN STDRED_COMBINECAT.PRO.  ADD IT BACK IN!!
;stdred_combinecat.pro:      temp.mag = cat.mag + 2.5*alog10(exptime)   ; CORRECT for exposure time!!!
allsrc.mag -= 2.5*alog10(allsrc.exptime)

;; Construct a CHIP structure
cindex = create_index(allsrc.frame)
uic = cindex.index[cindex.lo]
;uic = uniq(allsrc.frame,sort(allsrc.frame))
chschema = {field:'',file:'',expnum:'',frame:'',chip:0L,filter:'',exptime:0.0,utdate:'',uttime:'',airmass:0.0,apcor:0.0,night:'',mjd:0L,$
            nsrc:0L,allsrcindx:0LL,vertices_ra:dblarr(4),vertices_dec:dblarr(4)}
nchstr = n_elements(cindex.value)
chstr = replicate(chschema,nchstr)
chstr.field = allsrc[uic].stdfield
chstr.airmass = allsrc[uic].airmass
chstr.chip = allsrc[uic].chip
chstr.exptime = allsrc[uic].exptime
chstr.uttime = allsrc[uic].ut
chstr.mjd = allsrc[uic].mjd
chstr.expnum = allsrc[uic].expnum
chstr.frame = allsrc[uic].frame
chstr.filter = allsrc[uic].filter
chstr.nsrc = cindex.num
;; APCOR is 0.0 for standards since we have total photometry


print,'Running SMASHRED_CALIBRATE_FIELD_STANDARDS on ',field

;; COPIED ALL THIS FROM SMASHRED_PHOTCALIB.PRO
;-----------------------------------------------

; Field name
field = info[0].field

; Add MJD to CHSTR
if tag_exist(chstr,'MJD') eq 0 then begin
  add_tag,chstr,'MJD',0L,chstr
  dateobs = chstr.utdate+'T'+chstr.uttime
  for i=0,nchstr-1 do chstr[i].mjd=PHOTRED_GETMJD('','ctio',dateobs=dateobs[i])
endif
; Add CALIBRATED flag to CHSTR
add_tag,chstr,'calibrated',0B,chstr
; Add UBERCAL columns to CHSTR 
add_tag,chstr,'ubercal_magoffset',0.0,chstr
add_tag,chstr,'ubercal_flag',-1L,chstr
;; Add GAIACAL columns to CHSTR
;add_tag,chstr,'gaiacal_magoffset',0.0,chstr
;add_tag,chstr,'gaiacal_magofferr',99.99,chstr
; Add ZEROPOINT offset columns to CHSTR
add_tag,chstr,'zpcalib_magoffset',0.0,chstr
add_tag,chstr,'zpcalib_magofferr',99.99,chstr
add_tag,chstr,'zpcalibflag',0,chstr

; Load the transformation equations
trans_fitstr = MRDFITS(transfile,1,/silent)
trans_chipstr = MRDFITS(transfile,2,/silent)
trans_ntstr = MRDFITS(transfile,3,/silent)

; --- Add the transformation equation information to CHSTR ---

; Add the photometric transformation equation columns to CHSTR
if tag_exist(chstr,'photometric') eq 0 then begin
  newtags = ['photometric','badsoln','band','colband','colsign','zpterm','zptermsig','amterm',$
             'amtermsig','colterm','coltermsig','amcolterm','amcoltermsig','colsqterm','colsqtermsig']
  temp = chstr & undefine,chstr
  add_tags,temp,newtags,['0B','0B','""','""','-1',replicate('0.0d0',10)],chstr
  undefine,temp
endif
chstr.band = chstr.filter

; Transfer over the photometric transformation equations
chstr_mjdchipfilt = strtrim(chstr.mjd,2)+'-'+strtrim(chstr.chip,2)+'-'+strtrim(chstr.filter,2)
trans_mjdchipfilt = strtrim(trans_fitstr.mjd,2)+'-'+strtrim(trans_fitstr.chip,2)+'-'+strtrim(trans_fitstr.filter,2)
for i=0,nchstr-1 do begin
  MATCH,chstr_mjdchipfilt[i],trans_mjdchipfilt,ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    chstr[i].photometric = trans_fitstr[ind2].photometric
    chstr[i].badsoln = trans_fitstr[ind2].badsoln
    chstr[i].band = trans_fitstr[ind2].filter
    chstr[i].colband = trans_fitstr[ind2].colband
    chstr[i].colsign = trans_fitstr[ind2].colsign
    chstr[i].zpterm = trans_fitstr[ind2].zpterm
    chstr[i].zptermsig = trans_fitstr[ind2].zptermerr
    chstr[i].amterm = trans_fitstr[ind2].amterm
    chstr[i].amtermsig = trans_fitstr[ind2].amtermerr
    chstr[i].colterm = trans_fitstr[ind2].colterm
    chstr[i].coltermsig = trans_fitstr[ind2].coltermerr
    chstr[i].amcolterm = trans_fitstr[ind2].colamterm
    chstr[i].amcoltermsig = trans_fitstr[ind2].colamtermerr
    ;chstr[i].colsqterm = trans_fitstr[ind2].colsqterm
    ;chstr[i].colsqtermsig = trans_fitstr[ind2].colsqtermerr
  endif else begin  ; no match
    chstr[i].badsoln = 1  ; no solution
  endelse

  ; BADSOLN or NON-PHOTOMETRIC, zero-out all trans values except for COLTERM
  if chstr[i].badsoln eq 1 or chstr[i].photometric eq 0 then begin
    chstr[i].zpterm = 0.0
    chstr[i].zptermsig = 0.0
    chstr[i].amterm = 0.0
    chstr[i].amtermsig = 0.0
    chstr[i].amcolterm = 0.0
    chstr[i].amcoltermsig = 0.0
    chstr[i].colsqterm = 0.0
    chstr[i].colsqtermsig = 0.0
  endif
endfor

; Double-check that every CHSTR element has the essentials
;  COLBAND, COLSIGN, COLTERM COLTERMSIG
bdchstr = where(chstr.colband eq '',nbdchstr)
chstr_chipfilt = strtrim(chstr.chip,2)+'-'+strtrim(chstr.filter,2)
trans_chipfilt = strtrim(trans_chipstr.chip,2)+'-'+strtrim(trans_chipstr.filter,2)
for i=0,nbdchstr-1 do begin
  ; Get this information from the chip-specific
  ; transformation equation structure
  MATCH,chstr_chipfilt[bdchstr[i]],trans_chipfilt,ind1,ind2,/sort,count=nmatch
  chstr[bdchstr[i]].colband = trans_chipstr[ind2[0]].colband
  chstr[bdchstr[i]].colsign = trans_chipstr[ind2[0]].colsign
  chstr[bdchstr[i]].colterm = trans_chipstr[ind2[0]].colterm
  chstr[bdchstr[i]].coltermsig = trans_chipstr[ind2[0]].coltermerr
endfor

;; KEEP ONLY *PHOTOMETRIC* DATA
gdch = where(chstr.photometric eq 1,ngdch)
print,'Keeping ',strtrim(ngdch,2),' photometric chip observations'
chstr0 = chstr
chstr = chstr[gdch]
nchstr = n_elements(chstr)

;; Construct EXPOSURE structure
uif = uniq(chstr.expnum,sort(chstr.expnum))
fschema = {field:'',file:'',expnum:'',filter:'',exptime:0.0,utdate:'',uttime:'',airmass:0.0,apcor:0.0,night:'',mjd:0L}
fstr = replicate(fschema,n_elements(uif))
fstr.field = chstr[uif].field
fstr.airmass = chstr[uif].airmass
fstr.exptime = chstr[uif].exptime
fstr.uttime = chstr[uif].uttime
fstr.mjd = chstr[uif].mjd
fstr.expnum = chstr[uif].expnum
fstr.filter = chstr[uif].filter

; Copy PHOTOMETRIC and BADSOLN to FSTR
add_tag,fstr,'photometric',0,fstr
add_tag,fstr,'badsoln',0,fstr
for i=0,n_elements(fstr)-1 do begin
  MATCH,fstr[i].expnum,chstr.expnum,ind1,ind2,/sort,count=nmatch
  fstr[i].photometric = chstr[ind2[0]].photometric
  fstr[i].badsoln = chstr[ind2[0]].badsoln
endfor

; Give statistics on how many exposures have calibration information
; (per band) and photometric data
allfilters = ['u','g','r','i','z']
nocalib = intarr(5)
print,'--------------------------'
print,'FILTER   NEXP   NGOODEXP'
print,'=========================='
for i=0,n_elements(allfilters)-1 do begin
  MATCH,allfilters[i],fstr.filter,ind1,ind2,/sort,count=nmatch
  if nmatch gt 0 then begin
    gdexp = where(fstr[ind2].photometric eq 1 and fstr[ind2].badsoln eq 0,ngdexp)
  endif else ngdexp=0
  if nmatch gt 0 and ngdexp eq 0 then nocalib[i]=1    ; can't calibrate!
  print,allfilters[i],nmatch,ngdexp,format='(A4,I8,I8)'
endfor
print,'--------------------------'

; Some band CANNOT be calibrated, give average zeropoint
;  Still do this even if we are doing Gaia calibration
;  it can help to have a good first guess
bdnocalib = where(nocalib eq 1,nbdnocalib)
if nbdnocalib gt 0 then begin
  print,'!!!! CANNOT CALIBRATE ',strtrim(nbdnocalib,2),' FILTER(S): ',strjoin(allfilters[bdnocalib],' '),' !!!!'
  ;print,'Giving them average zeropoints for that band'
  for i=0,nbdnocalib-1 do begin
    MATCH,trans_ntstr.filter,allfilters[bdnocalib[i]],ind1,ind2,/sort,count=nmatch
    zpterm = median([trans_ntstr[ind1].zpterm])
    zptermsig = median([trans_ntstr[ind1].zptermsig])
    MATCH,chstr.filter,allfilters[bdnocalib[i]],ind1,ind2,/sort,count=nmatch
    chstr[ind1].zpterm = zpterm
    chstr[ind1].zptermsig = zptermsig
    ; Set COLOR TERMS to 0.0 for these chips
    ;  we can still use these if we have Gaia calibration
    if not keyword_set(usegaia) then begin
      chstr[ind1].colterm = 0.0
      chstr[ind1].coltermsig = 0.0
    endif
  endfor
endif


;; Sort ALLSRC by the chip/frame
;; ONLY KEEP PHOTOMETRIC ONES
MATCH,chstr.frame,cindex.value,chind1,chind2,/sort
cnt = 0LL
newind = lonarr(total(chstr.nsrc))
for i=0,n_elements(chstr)-1 do begin
  ind1 = cindex.index[cindex.lo[chind2[i]]:cindex.hi[chind2[i]]]
  ind1 = ind1[sort(ind1)]
  nind1 = n_elements(ind1)
  newind[cnt:cnt+nind1-1] = ind1
  chstr[i].allsrcindx = cnt
  cnt += nind1
  ;; Calculate the vertices from the detections
  rar = minmax(all[ind1].ra)
  decr = minmax(all[ind1].dec)  
  chstr[i].vertices_ra = [rar[0],rar[1],rar[1],rar[0]]
  chstr[i].vertices_dec = [decr[0],decr[0],decr[1],decr[1]]
endfor
allsrc = allsrc[newind]



;; Create the ALLOBJ structure
indx = create_index(allsrc.id)
nuexp = n_elements(fstr)
nan = !values.f_nan
dnan = !values.d_nan
allobj_schema = {id:'',ra:0.0d0,dec:0.0d0,raerr:99.99,decerr:99.99,rascatter:99.99,decscatter:99.99,ndet:0,depthflag:0B,$
                 srcindx:lonarr(nuexp)-1,srcfindx:lonarr(nuexp)-1,u:99.99,uerr:9.99,uscatter:99.99,ndetu:0,g:99.99,gerr:9.99,$
                 gscatter:99.99,ndetg:0,r:99.99,rerr:9.99,rscatter:99.9,ndetr:0,i:99.99,ierr:9.99,iscatter:99.99,ndeti:0,$
                 z:99.99,zerr:9.99,zscatter:99.99,ndetz:0,chi:nan,sharp:nan,flag:-1,prob:nan,ebv:99.99}
allobj = replicate(allobj_schema,n_elements(indx.num))
allobj.id = indx.value
for i=0,n_elements(indx.value)-1 do begin
  ;; SRCINDX has all detections pushed to the front
  ind1 = indx.index[ lindgen(indx.num[i])+indx.lo[i] ]
  ind1 = ind1[sort(ind1)]
  allobj[i].srcindx = ind1
  ;; SRCFINDX have the detections in their appropriate exposure column
  MATCH,fstr.expnum,allsrc[ind1].expnum,index1,index2,/sort
  allobj[i].srcfindx[index1] = ind1[index2]
  ;; CMBINDX
  allsrc[ind1].cmbindx = i
endfor
allobj.ra = all[indx.index[indx.lo]].ra
allobj.dec = all[indx.index[indx.lo]].dec
allobj.ndet =  indx.num





;; NOW DO THE CALIBRATION AND AVERAGING
;---------------------------------------
SMASHRED_AVERAGEMORPHCOORD,fstr,chstr,allsrc,allobj
;; Calibrate
SMASHRED_APPLY_PHOTTRANSEQN,fstr,chstr,allsrc,allobj

; Compute average morphology and coordinate values
print,'Calculating average morphology and coordinate parameters'
SMASHRED_AVERAGEMORPHCOORD,fstr,chstr,allsrc,allobj

;; Compute exposure map
;print,'Computing the exposure maps'
;SMASHRED_COMPUTE_EXPMAP,field,chstr,redo=redo,outputdir=outputdir

;; Set non-detections based on the exposure map
;print,'Setting non-detections based on the exposure maps'
;SMASHRED_SET_NONDETECTIONS,field,allobj,dir=outputdir

; Calculate extinction
print,'Getting SFD E(B-V) extinctions'
GLACTC,allobj.ra,allobj.dec,2000.0,lon,lat,1,/deg
ebv = DUST_GETVAL(lon,lat,/interp,/noloop)
; add EBV, G0 and I0 to the catalogs
allobj.ebv = ebv


;; COLOR-COLOR DIAGRAM LOOKS GOOD
;str=mrdfits('~/smash/reduction/catalogs/final/v6/stars1/Field246_allobj_deep_stars.fits.gz',1)  
;loadcol,3
;hess,str.u-str.g,str.g-str.r,dx=0.02,dy=0.02,xr=[-1,4],yr=[-1,2],/log
;loadct,39
;oplot,allobj.u-allobj.g,allobj.g-allobj.r,ps=1,co=250


; Write out the final files
print,'Writing combined file to ',outfile
MWRFITS,fstr,outfile+'_exposures.fits',/create
MWRFITS,chstr,outfile+'_chips.fits',/create
MWRFITS,allsrc,outfile+'_allsrc.fits',/create
MWRFITS,allobj,outfile+'_allobj.fits',/create
; Compress
if keyword_set(compress) then begin
  print,'Compressing output files'
  spawn,['gzip','-f',outfile+'_exposures.fits'],out,errout,/noshell
  spawn,['gzip','-f',outfile+'_chips.fits'],out,errout,/noshell
  spawn,['gzip','-f',outfile+'_allsrc.fits'],out,errout,/noshell
  spawn,['gzip','-f',outfile+'_allobj.fits'],out,errout,/noshell
  spawn,['gzip','-f',outfile+'_expmap.fits'],out,errout,/noshell
endif

; Make bright allobj catalog
;SMASHRED_MAKE_BRIGHTCAT,field,redo=redo,dir=outputdir

; Print processing time
dt = systime(1)-t0
print,'' & print,'Processing time = ',strtrim(string(dt/60.0,format='(F20.3)'),2),' min.'

;stop

end
