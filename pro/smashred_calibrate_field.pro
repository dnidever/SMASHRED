;+
;
; SMASHRED_CALIBRATE_FIELD
;
; This program calibrates all of the photometry for a single SMASH field
; using transformation equations and ubercal techniques.
;
; INPUTS:
;  field      The field name, e.g. "Field24"
;  =version   The version of the final catalogs.  This is only used
;               if OUTPUTDIR is not input.
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
;  IDL>smashred_calibrate_field,'Field24'
;
; By D.Nidever March 2016
;-

pro smashred_calibrate_field,field,version=version,transfile=transfile,reduxdir=reduxdir,outputdir=outputdir,$
                             usegaia=usegaia,redo=redo,compress=compress,silent=silent,error=error

undefine,error
t0 = systime(1)

; Not enough inputs
if n_elements(field) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_calibrate_field,field,version=version,reduxdir=reduxdir,outputdir=outputdir,redo=redo,'
  print,'                                  usegaia=usegaia,compress=compress,error=error'
  return
endif

; Defaults
if n_elements(reduxdir) eq 0 then reduxdir='/data/smash/cp/red/photred/'
if file_test(reduxdir,/directory) eq 0 then begin
  error = reduxdir+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif
if n_elements(outputdir) eq 0 then begin
  outputdir = reduxdir+'catalogs/final/'
  if n_elements(version) gt 0 then outputdir+=version+'/'
endif
if file_test(outputdir,/directory) eq 0 then begin
  if not keyword_set(silent) then print,outputdir+' does NOT exist.  Creating it.'
  FILE_MKDIR,outputdir
endif
if n_elements(transfile) eq 0 then transfile='/data/smash/cp/red/photred/stdred/smashred_transphot_eqns.fits'
; Temporary directory
tmpdir = outputdir+'/tmp/'
if file_test(tmpdir,/directory) eq 0 then FILE_MKDIR,tmpdir
; Compression
if n_elements(compress) eq 0 then compress=1
; GAIA photometry
if n_elements(usegaia) eq 0 then usegaia=0

; Check if the final files already exists
outfile = outputdir+field+'_combined'
testfiles = outfile+['_exposures','_chips','_allsrc','_allobj']+'.fits'
ntestfiles = n_elements(testfiles)
if (total(file_test(testfiles)) eq ntestfiles or total(file_test(testfiles+'.gz')) eq ntestfiles) and $
  not keyword_set(redo) then begin
  print,'Final output files already exist for ',field,' and /redo NOT set.'
  return
endif

; Get reduction info
print,'Getting reduction information'
SMASHRED_GETREDINFO,allinfo,/silent
gdinfo = where(allinfo.field eq field,ngdinfo)
if ngdinfo eq 0 then begin
  error = 'No reduced exposures for '+field
  if not keyword_set(silent) then print,error
  return
endif
info = allinfo[gdinfo]
ninfo = ngdinfo

; MAYBE GET THE CATALOGS FROM THE PHOTRED DIRECTOREIS DIRECTLY!!!!???
; HOW DO WE KNOW IF THE CATALOGS/PHOT FILES ARE CALIBRATED OR INST??
;  Want instrumental photometry.  Maybe only do instrumental mags with PHOTRED
;  Maybe just us the AST files!!!??  We are going to calibrate the photometry
;  anyway so why not also correct for exptime and apcor as well??!!
;  Another option would be to always include the instrumental phot in
;  the PHOT files and use those instead.

; The various steps for this calibration:
; 1.) Load all of the chip data from .phot files
; 2.) Crossmatch all of the sources and build ALLSRC and ALLOBJ structures
; 3.) Calibrate:  smashred_photcalib.pro
;     a) measure photometric offsets between chip pairs (using calibrated photometry)
;     b) determine relative magnitude offsets per chip using ubercal and
;     c) set zeropoint per band with "photometric" data or 0.9m data
;     d) calculate average photometry per object per filter
;     e) calibrate instrumental photometry with trans eqns. and ubercal mag offsets
;
; HOW DOES THE 0.9M DATA FIT IN????  Can't be done just as an afterburner,
; since this will affect the mean mags which will in turn affect the color
; terms.  So needs to be a 
;
; 4.) Output the data and compress

; Load the transformation equations for all chips, nights and filters
; Load the list of photometric/non-photometric nights
; Load 0.9m data

print,'Running SMASHRED_CALIBRATE_FIELD on ',field
print,strtrim(ninfo,2),' PHOTRED catalogs found'
print,info.file


; Check for temporarily saved crossmatch file
crossmatchfile = tmpdir+field+'_crossmatch.dat'
if file_test(crossmatchfile) eq 1 and not keyword_set(redo) then goto,crossmatch  ; skip loading


; Loop through the catalogs
print,'--------------------------------------------------'
print,'--- STEP 1. Load all of the PHOTRED photometry ---'
print,'=================================================='
ninfo = n_elements(info)
undefine,allfstr,allchstr,allsrc
For c=0,ninfo-1 do begin
  info1 = info[c]
  print,strtrim(c+1,2),'/',strtrim(ninfo,2),' ',info1.file
  undefine,chstr1,allsrc1
  SMASHRED_LOAD_CATPHOT,info1,chstr1,allsrc1,/useast,reduxdir=reduxdir,redo=redo,outputdir=outputdir

  ; Some good data to add
  if n_elements(allsrc1) gt 0 then begin
    ; offset indices
    chstr1.allsrcindx += n_elements(allsrc)
    allsrc1.chipindx += n_elements(allchstr)

    ; Combine structures
    push,allfstr,*info1.fstr
    push,allchstr,chstr1
    push,allsrc,allsrc1
  endif
Endfor ; catalog loop
; Kludge!  Remove duplicate exposures in short/deep for Field130
if field eq 'Field130' then begin
  print,'KLUDGE! Removing duplicate exposures in short/deep for Field130'
  bd = where((allfstr.expnum eq '00423440' and allfstr.alf_nsources lt 0) or $
             (allfstr.expnum eq '00426607' and allfstr.alf_nsources lt 0),nbd)
  REMOVE,bd,allfstr
endif
; Rename
fstr = allfstr
chstr = allchstr
undefine,allfstr,allchstr

print,'---------------------------------------------------------------------'
print,'--- STEP 2. Crossmatch all of the sources and build ALLSRC/ALLOBJ ---'
print,'====================================================================='
crossmatch:
SMASHRED_CROSSMATCH,field,fstr,chstr,allsrc,allobj,reduxdir=reduxdir,redo=redo,outputdir=outputdir


print,'-----------------------------------------------'
print,'--- STEP 3. Calibrate all of the photometry ---'
print,'==============================================='
SMASHRED_PHOTCALIB,info,fstr,chstr,allsrc,allobj,transfile=transfile,usegaia=usegaia,reduxdir=reduxdir,$
                   outputdir=outputdir


; Compute average morphology and coordinate values
print,'Calculating average morphology and coordinate parameters'
SMASHRED_AVERAGEMORPHCOORD,fstr,chstr,allsrc,allobj

; Compute exposure map
print,'Computing the exposure maps'
SMASHRED_COMPUTE_EXPMAP,field,chstr,redo=redo,outputdir=outputdir

; Set non-detections based on the exposure map
print,'Setting non-detections based on the exposure maps'
SMASHRED_SET_NONDETECTIONS,field,allobj,dir=outputdir

; Calculate extinction
print,'Getting SFD E(B-V) extinctions'
GLACTC,allobj.ra,allobj.dec,2000.0,lon,lat,1,/deg
ebv = DUST_GETVAL(lon,lat,/interp,/noloop)
; add EBV, G0 and I0 to the catalogs
allobj.ebv = ebv

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
SMASHRED_MAKE_BRIGHTCAT,field,redo=redo,dir=outputdir

; Print processing time
dt = systime(1)-t0
print,'' & print,'Processing time = ',strtrim(string(dt/60.0,format='(F20.3)'),2),' min.'

;stop

end
