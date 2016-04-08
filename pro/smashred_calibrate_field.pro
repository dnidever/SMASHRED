;+
;
; SMASHRED_CALIBRATE_FIELD
;
; This program calibrates all of the photometry for a single SMASH field
; using transformation equations and ubercal techniques.
;
; INPUTS:
;  field      The field name, e.g. "Field24"
;  =reduxdir  The base directory for the PHOTRED reductions.  The
;               default is "/data/smash/cp/red/photred/"
;  =outputdir The output directory for the catalogs.
;  /compress  Gzip compress the output FITS files.  This is the default.
;  /redo      Redo a field that was already calibrated.
;  /silent    Don't print anything to the screen.
;
; OUTPUTS:
;  Four binary FITS files are created in the output directory:
;   FIELD_exposures.fits - information on each unique exposure
;   FIELD_chips.fits - information on each unique chip
;   FIELD_allsrc.fits - information on each unique source detection
;   FIELD_allobj.fits - information on each unique object including
;                       average magnitudes.
;   =error    The error message, if one occurred.
;
; USAGE:
;  IDL>smashred_calibrate,'Field24'
;
; By D.Nidever March 2016
;-

pro smashred_calibrate_field,field,reduxdir=reduxdir,outputdir=outputdir,redo=redo,$
                             compress=compress,silent=silent,error=error

undefine,error

; Not enough inputs
if n_elements(field) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - smashred_calibrate_field,field,reduxdir=reduxdir,outputdir=outputdir,redo=redo,'
  print,'                                  compress=compress,error=error'
  return
endif

; Defaults
if n_elements(reduxdir) eq 0 then reduxdir='/data/smash/cp/red/photred/'
if file_test(reduxdir,/directory) eq 0 then begin
  error = reduxdir+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif
if n_elements(outputdir) eq 0 then outputdir=reduxdir+'catalogs/final/'
if file_test(outputdir,/directory) eq 0 then begin
  if not keyword_set(silent) then print,outputdir+' does NOT exist.  Creating it.'
  FILE_MKDIR,outputdir
endif
; Temporary directory
tmpdir = outputdir+'/tmp/'
if file_test(tmpdir,/directory) eq 0 then FILE_MKDIR,tmpdir
; Compression
if n_elements(compress) eq 0 then compress=1

; Check if the final files already exists
outfile = outputdir+field+'_combined'
testfiles = outfile+['_exposures','_chips','_allsrc','_allobj']+'.fits'
ntestfiles = n_elements(testfiles)
if (total(file_test(testfiles)) eq ntestfiles or total(file_test(testfiles+'.gz')) eq ntestfiles) and $
  not keyword_set(redo) then begin
  print,'Final output files arleady exist for ',field,' and /redo NOT set.'
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

goto,crossmatch

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
  SMASHRED_LOAD_CATPHOT,info1,chstr1,allsrc1,/useast,reduxdir=reduxdir,redo=redo

  ; offset indices
  chstr1.allsrcindx += n_elements(allsrc)
  allsrc1.chipindx += n_elements(allchstr)

  ; Combine structures
  push,allfstr,*info1.fstr
  push,allchstr,chstr1
  push,allsrc,allsrc1
Endfor ; catalog loop

; Rename
fstr = allfstr
chstr = allchstr
undefine,allfstr,allchstr

print,'---------------------------------------------------------------------'
print,'--- STEP 2. Crossmatch all of the sources and build ALLSRC/ALLOBJ ---'
print,'====================================================================='
crossmatch:
;SMASHRED_CROSSMATCH,field,fstr,chstr,allsrc,allobj
;save,fstr,chstr,allsrc,allobj,file=tmpdir+field+'_crossmatch.dat'
print,'restoring temporary allsrc/allobj file'
restore,tmpdir+field+'_crossmatch.dat'


print,'-----------------------------------------------'
print,'--- STEP 3. Calibrate all of the photometry ---'
print,'==============================================='
SMASHRED_PHOTCALIB,fstr,chstr,allsrc,allobj


; Compute average morphology values
print,'Calculating average morphology parameters'
SMASHRED_AVERAGEMORPH,fstr,chstr,allsrc,allobj

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
  spawn,['gzip',outfile+'_exposures.fits'],out,errout,/noshell
  spawn,['gzip',outfile+'_chips.fits'],out,errout,/noshell
  spawn,['gzip',outfile+'_allsrc.fits'],out,errout,/noshell
  spawn,['gzip',outfile+'_allobj.fits'],out,errout,/noshell
endif


;stop

end
